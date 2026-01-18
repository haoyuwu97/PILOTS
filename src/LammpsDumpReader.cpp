#include "pilots/io/LammpsDumpReader.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <chrono>
#include <thread>

#include "pilots/util/Parse.hpp"

namespace pilots {

namespace {

inline std::runtime_error die(const std::string& msg) {
  return std::runtime_error("LammpsDumpReader: " + msg);
}

struct NeedMoreData {};

} // namespace

LammpsDumpReader::LammpsDumpReader(const std::string& path)
: ifs_(path) {
  if (!ifs_) {
    throw die("failed to open file: " + path);
  }
  ifs_.exceptions(std::ios::badbit);
}

void LammpsDumpReader::set_follow(bool enable, int poll_ms) {
  follow_ = enable;
  poll_ms_ = (poll_ms > 0) ? poll_ms : 200;
}

void LammpsDumpReader::set_allow_header_change(bool enable) {
  allow_header_change_ = enable;
}

void LammpsDumpReader::set_stop_flag(const volatile std::sig_atomic_t* flag) {
  stop_flag_ = flag;
}

void LammpsDumpReader::seek_offset(std::uint64_t off) {
  ifs_.clear();
  ifs_.seekg(static_cast<std::streamoff>(off), std::ios::beg);
  if (!ifs_) {
    throw die("seek_offset failed (offset out of range?)");
  }
  last_good_off_ = off;
}

bool LammpsDumpReader::next(Frame& frame) {
  while (true) {
    const std::streampos pos = ifs_.tellg();
    try {

        if (!read_header_(frame)) {
          if (follow_) throw NeedMoreData{};
          return false;
        }

        // ATOMS header
        std::string line;
        if (!std::getline(ifs_, line)) {
          if (follow_ && ifs_.eof()) throw NeedMoreData{};
          throw die("unexpected EOF while reading ITEM: ATOMS header");
        }
        if (line.rfind("ITEM: ATOMS", 0) != 0) {
          // In follow mode, a partially-written header line at EOF should be treated as incomplete data.
          if (follow_ && ifs_.eof()) throw NeedMoreData{};
          throw die("expected 'ITEM: ATOMS ...' but got: " + line);
        }
        parse_atoms_header_(line);

        // Update per-frame presence flags (based on the validated header schema).
        frame.has_type = (col_.type >= 0);
        frame.has_mol  = (col_.mol >= 0);

        // Optional core arrays are allocated on-demand (Frame::resize() runs before ATOMS header
        // parsing on the first frame).
        if (first_frame_) {
          frame.ensure_type_storage();
          frame.ensure_mol_storage();
        }

        // Bind the reader-owned extra schema to the frame and mark present fields for this frame.
        frame.bind_extra_schema(&extra_schema_);
        frame.clear_extra_presence();
        std::vector<double*> extra_ptrs;
        extra_ptrs.resize(extra_specs_.size(), nullptr);
        for (std::size_t e = 0; e < extra_specs_.size(); ++e) {
          const auto fid = extra_specs_[e].fid;
          frame.mark_extra_present(fid);
          extra_ptrs[e] = frame.ensure_extra_storage(fid).data();
        }

        if (!col_.have_unwrapped() && !col_.have_wrapped_with_images()) {
          throw die("dump must contain either (xu,yu,zu) or (x,y,z) with (ix,iy,iz)");
        }

        if (frame.natoms != natoms_) {
          throw die("natoms mismatch between frames");
        }

        // Read atoms (N lines)
        std::string atom_line;
        std::vector<std::int64_t> ids_tmp;
        if (first_frame_) {
          ids_tmp.resize(natoms_);
        }

        // Local storage for extra field values for one atom line.
        // Each extra field corresponds to a concrete column present in the ATOMS header.
        std::vector<double> extra_values(extra_specs_.size(), 0.0);

        const int expected_cols = static_cast<int>(extra_slot_by_col_.size());
        for (std::size_t row = 0; row < natoms_; ++row) {
          if (!std::getline(ifs_, atom_line)) {
      if (follow_ && ifs_.eof()) throw NeedMoreData{};
      throw die("unexpected EOF while reading atom lines");
    }

          const char* p = atom_line.data();
          const char* end = p + atom_line.size();

          std::string_view tok;
          int col_idx = 0;

          std::int64_t id = -1;
          int type = 0;
          std::int64_t mol = 0;

          double xu = 0.0, yu = 0.0, zu = 0.0;
          double x = 0.0, y = 0.0, z = 0.0;
          std::int64_t ix = 0, iy = 0, iz = 0;

          bool have_id = false;
          // We only store type/mol from the first frame (canonical ordering). For subsequent frames,
          // type/mol may still appear in the dump, but we do not require parsing them.
          bool have_type = (!first_frame_) || (col_.type < 0);
          bool have_mol  = (!first_frame_) || (col_.mol  < 0);
          bool have_xu = !col_.have_unwrapped();
          bool have_yu = !col_.have_unwrapped();
          bool have_zu = !col_.have_unwrapped();
          bool have_x  = col_.have_unwrapped();
          bool have_y  = col_.have_unwrapped();
          bool have_z  = col_.have_unwrapped();
          bool have_ix = col_.have_unwrapped();
          bool have_iy = col_.have_unwrapped();
          bool have_iz = col_.have_unwrapped();

          while (next_token(p, end, tok)) {
            if (col_idx == col_.id) {
              if (!parse_int(tok, id)) {
                if (follow_ && ifs_.eof()) throw NeedMoreData{};
                throw die("failed to parse id");
              }
              have_id = true;
            }
            // Extra scalar fields (auto-captured from header)
            if (!extra_specs_.empty() && col_idx >= 0 && col_idx < static_cast<int>(extra_slot_by_col_.size())) {
              const int slot = extra_slot_by_col_[col_idx];
              if (slot >= 0) {
                double v = 0.0;
                if (!parse_double(tok, v)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse extra field '" + extra_specs_[static_cast<std::size_t>(slot)].name + "'");
                }
                extra_values[static_cast<std::size_t>(slot)] = v;
              }
            }
            if (first_frame_) {
              if (col_.type >= 0 && col_idx == col_.type) {
                if (!parse_int(tok, type)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse type");
                }
                have_type = true;
              } else if (col_.mol >= 0 && col_idx == col_.mol) {
                if (!parse_int(tok, mol)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse mol");
                }
                have_mol = true;
              }
            }

            if (col_.have_unwrapped()) {
              if (col_idx == col_.xu) {
                if (!parse_double(tok, xu)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse xu");
                }
                have_xu = true;
              } else if (col_idx == col_.yu) {
                if (!parse_double(tok, yu)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse yu");
                }
                have_yu = true;
              } else if (col_idx == col_.zu) {
                if (!parse_double(tok, zu)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse zu");
                }
                have_zu = true;
              }
            } else {
              if (col_idx == col_.x) {
                if (!parse_double(tok, x)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse x");
                }
                have_x = true;
              } else if (col_idx == col_.y) {
                if (!parse_double(tok, y)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse y");
                }
                have_y = true;
              } else if (col_idx == col_.z) {
                if (!parse_double(tok, z)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse z");
                }
                have_z = true;
              } else if (col_idx == col_.ix) {
                if (!parse_int(tok, ix)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse ix");
                }
                have_ix = true;
              } else if (col_idx == col_.iy) {
                if (!parse_int(tok, iy)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse iy");
                }
                have_iy = true;
              } else if (col_idx == col_.iz) {
                if (!parse_int(tok, iz)) {
                  if (follow_ && ifs_.eof()) throw NeedMoreData{};
                  throw die("failed to parse iz");
                }
                have_iz = true;
              }
            }

            ++col_idx;
          }

          if (!have_id || !have_type || !have_mol || !have_xu || !have_yu || !have_zu || !have_x || !have_y || !have_z || !have_ix || !have_iy || !have_iz) {
            if (follow_ && ifs_.eof()) throw NeedMoreData{};
            throw die("atom line missing required fields (incomplete line)");
          }

          if (expected_cols > 0 && col_idx != expected_cols) {
            // If the line has too few/many columns compared to the ATOMS header, treat as incomplete in follow mode at EOF.
            if (follow_ && ifs_.eof()) throw NeedMoreData{};
            throw die("atom line column count mismatch vs ATOMS header");
          }

          double ux, uy, uz;
          if (col_.have_unwrapped()) {
            ux = xu; uy = yu; uz = zu;
          } else {
            auto u = frame.box.unwrap(x, y, z, ix, iy, iz);
            ux = u[0]; uy = u[1]; uz = u[2];
          }

          if (first_frame_) {
            // Canonical order = row order in the first frame
            ids_tmp[row] = id;
            frame.id[row] = id;
            if (frame.has_type) frame.type[row] = type;
            if (frame.has_mol)  frame.mol[row] = mol;
            frame.xu[row] = ux;
            frame.yu[row] = uy;
            frame.zu[row] = uz;

            // Store extras (canonical row order on first frame)
            for (std::size_t e = 0; e < extra_specs_.size(); ++e) {
              extra_ptrs[e][row] = extra_values[e];
              extra_values[e] = 0.0;
            }
          } else {
            const std::size_t idx = lookup_idx_(id);
            frame.xu[idx] = ux;
            frame.yu[idx] = uy;
            frame.zu[idx] = uz;

            for (std::size_t e = 0; e < extra_specs_.size(); ++e) {
              extra_ptrs[e][idx] = extra_values[e];
              extra_values[e] = 0.0;
            }
          }
        }

        if (first_frame_) {
          ids_first_ = std::move(ids_tmp);
          build_id_map_();
          first_frame_ = false;
        }

        // Record checkpoint offset after a fully-read frame.
        // Note: tellg() may return -1 if eofbit is set; clear/restore to obtain a stable offset.
        std::streamoff off = static_cast<std::streamoff>(ifs_.tellg());
        if (off < 0) {
          const auto st = ifs_.rdstate();
          ifs_.clear();
          off = static_cast<std::streamoff>(ifs_.tellg());
          ifs_.setstate(st);
        }
        if (off >= 0) last_good_off_ = static_cast<std::uint64_t>(off);

        return true;
    } catch (const NeedMoreData&) {
      if (!follow_) return false;
      if (stop_flag_ && *stop_flag_) return false;
      ifs_.clear();
      if (pos != std::streampos(-1)) {
        ifs_.seekg(pos);
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(poll_ms_));
      continue;
    }
  }
}


bool LammpsDumpReader::read_header_(Frame& frame) {
  std::string line;

  // Strict ITEM sequence: the next non-empty line must be ITEM: TIMESTEP.
  while (std::getline(ifs_, line)) {
    if (line.empty()) continue;
    if (line.rfind("ITEM: TIMESTEP", 0) == 0) break;
    // In follow mode, an incomplete header line at EOF should be treated as incomplete data.
    if (follow_ && ifs_.eof()) throw NeedMoreData{};
    throw die("expected 'ITEM: TIMESTEP' but got: " + line);
  }

  if (!ifs_) {
    return false; // EOF
  }

  // timestep value
  if (!std::getline(ifs_, line)) {
    if (follow_ && ifs_.eof()) throw NeedMoreData{};
    throw die("unexpected EOF while reading timestep");
  }
  {
    std::vector<std::string_view> toks;
    split_ws(line, toks);
    if (toks.empty()) throw die("empty timestep line");
    std::int64_t ts = 0;
    if (!parse_int(toks[0], ts)) {
      if (follow_ && ifs_.eof()) throw NeedMoreData{};
      throw die("failed to parse timestep");
    }
    frame.timestep = ts;
  }

  // ITEM: NUMBER OF ATOMS
  if (!std::getline(ifs_, line)) {
    if (follow_ && ifs_.eof()) throw NeedMoreData{};
    throw die("unexpected EOF while reading NUMBER OF ATOMS header");
  }
  if (line.rfind("ITEM: NUMBER OF ATOMS", 0) != 0) {
    throw die("expected 'ITEM: NUMBER OF ATOMS' but got: " + line);
  }
  if (!std::getline(ifs_, line)) {
    if (follow_ && ifs_.eof()) throw NeedMoreData{};
    throw die("unexpected EOF while reading natoms");
  }
  {
    std::vector<std::string_view> toks;
    split_ws(line, toks);
    if (toks.empty()) throw die("empty natoms line");
    std::int64_t n = 0;
    if (!parse_int(toks[0], n) || n <= 0) {
      if (follow_ && ifs_.eof()) throw NeedMoreData{};
      throw die("failed to parse natoms");
    }
    if (first_frame_) {
      natoms_ = static_cast<std::size_t>(n);
      frame.resize(natoms_);
    } else {
      if (static_cast<std::size_t>(n) != natoms_) {
        throw die("natoms changed between frames; unsupported");
      }
      frame.natoms = natoms_;
    }
  }

  // ITEM: BOX BOUNDS ...
  if (!std::getline(ifs_, line)) {
    if (follow_ && ifs_.eof()) throw NeedMoreData{};
    throw die("unexpected EOF while reading BOX BOUNDS header");
  }
  if (line.rfind("ITEM: BOX BOUNDS", 0) != 0) {
    throw die("expected 'ITEM: BOX BOUNDS' but got: " + line);
  }

  if (!read_box_(frame, line)) {
    throw die("failed to parse BOX BOUNDS");
  }

  return true;
}

bool LammpsDumpReader::read_box_(Frame& frame, const std::string& header_line) {
  std::vector<std::string_view> toks;
  split_ws(header_line, toks);

  bool has_xy = false, has_xz = false, has_yz = false;
  for (auto& t : toks) {
    if (t == "xy") has_xy = true;
    if (t == "xz") has_xz = true;
    if (t == "yz") has_yz = true;
  }
  const bool triclinic = has_xy || has_xz || has_yz;
  frame.box.triclinic = triclinic;

  std::string line;
  if (!std::getline(ifs_, line)) { if (follow_ && ifs_.eof()) throw NeedMoreData{}; return false; }
  {
    std::vector<std::string_view> v;
    split_ws(line, v);
    if ((!triclinic && v.size() < 2) || (triclinic && v.size() < 3)) {
      if (follow_ && ifs_.eof()) throw NeedMoreData{};
      throw die("BOX BOUNDS line 1 has wrong number of fields");
    }
    if (!parse_double(v[0], frame.box.xlo) || !parse_double(v[1], frame.box.xhi)) {
      if (follow_ && ifs_.eof()) throw NeedMoreData{};
      throw die("failed to parse xlo/xhi");
    }
    if (triclinic) {
      if (!parse_double(v[2], frame.box.xy)) {
        if (follow_ && ifs_.eof()) throw NeedMoreData{};
        throw die("failed to parse xy");
      }
    } else {
      frame.box.xy = 0.0;
    }
  }

  if (!std::getline(ifs_, line)) { if (follow_ && ifs_.eof()) throw NeedMoreData{}; return false; }
  {
    std::vector<std::string_view> v;
    split_ws(line, v);
    if ((!triclinic && v.size() < 2) || (triclinic && v.size() < 3)) {
      if (follow_ && ifs_.eof()) throw NeedMoreData{};
      throw die("BOX BOUNDS line 2 has wrong number of fields");
    }
    if (!parse_double(v[0], frame.box.ylo) || !parse_double(v[1], frame.box.yhi)) {
      if (follow_ && ifs_.eof()) throw NeedMoreData{};
      throw die("failed to parse ylo/yhi");
    }
    if (triclinic) {
      if (!parse_double(v[2], frame.box.xz)) {
        if (follow_ && ifs_.eof()) throw NeedMoreData{};
        throw die("failed to parse xz");
      }
    } else {
      frame.box.xz = 0.0;
    }
  }

  if (!std::getline(ifs_, line)) { if (follow_ && ifs_.eof()) throw NeedMoreData{}; return false; }
  {
    std::vector<std::string_view> v;
    split_ws(line, v);
    if ((!triclinic && v.size() < 2) || (triclinic && v.size() < 3)) {
      if (follow_ && ifs_.eof()) throw NeedMoreData{};
      throw die("BOX BOUNDS line 3 has wrong number of fields");
    }
    if (!parse_double(v[0], frame.box.zlo) || !parse_double(v[1], frame.box.zhi)) {
      if (follow_ && ifs_.eof()) throw NeedMoreData{};
      throw die("failed to parse zlo/zhi");
    }
    if (triclinic) {
      if (!parse_double(v[2], frame.box.yz)) {
        if (follow_ && ifs_.eof()) throw NeedMoreData{};
        throw die("failed to parse yz");
      }
    } else {
      frame.box.yz = 0.0;
    }
  }

  // For triclinic boxes, LAMMPS writes BOX BOUNDS as the *bounding box* values (xlo_bound/xhi_bound, etc.)
  // along with tilt factors (xy,xz,yz). To build a consistent H-matrix representation for
  // fractional wrapping/min-image, we reconstruct the true box origin (xlo,ylo,zlo) and
  // lengths (Lx,Ly,Lz) from the bounds and tilts.
  // See LAMMPS triclinic box documentation: xlo_bound = xlo + min(0,xy,xz,xy+xz), etc.
  if (triclinic) {
    const double xlo_bound = frame.box.xlo;
    const double xhi_bound = frame.box.xhi;
    const double ylo_bound = frame.box.ylo;
    const double yhi_bound = frame.box.yhi;
    const double zlo_bound = frame.box.zlo;
    const double zhi_bound = frame.box.zhi;
    const double xy = frame.box.xy;
    const double xz = frame.box.xz;
    const double yz = frame.box.yz;

    const double min_x = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
    const double max_x = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
    const double min_y = std::min(0.0, yz);
    const double max_y = std::max(0.0, yz);

    frame.box.xlo = xlo_bound - min_x;
    frame.box.xhi = xhi_bound - max_x;
    frame.box.ylo = ylo_bound - min_y;
    frame.box.yhi = yhi_bound - max_y;
    frame.box.zlo = zlo_bound;
    frame.box.zhi = zhi_bound;
  }

  return true;
}

void LammpsDumpReader::parse_atoms_header_(const std::string& atoms_header_line) {
  std::vector<std::string_view> toks;
  split_ws(atoms_header_line, toks);
  if (toks.size() < 3) {
    if (follow_ && ifs_.eof()) throw NeedMoreData{};
    throw die("ATOMS header too short");
  }

  // Schema validation: enforce stable header across frames.
  // We compare the full field list (in order) excluding the leading "ITEM: ATOMS" tokens.
  {
    std::vector<std::string> fields;
    fields.reserve(toks.size() - 2);
    for (std::size_t i = 2; i < toks.size(); ++i) fields.emplace_back(std::string(toks[i]));
    if (!atoms_fields_set_) {
      atoms_fields_first_ = fields;
      atoms_fields_set_ = true;
    } else {
      bool same = (fields.size() == atoms_fields_first_.size());
      if (same) {
        for (std::size_t i = 0; i < fields.size(); ++i) {
          if (fields[i] != atoms_fields_first_[i]) { same = false; break; }
        }
      }
      if (!same) {
        if (!allow_header_change_) {
          if (follow_ && ifs_.eof()) throw NeedMoreData{};
          throw die("ATOMS header changed between frames (field name/order mismatch)");
        }
        // Accept header changes: rebuild column mapping for this frame.
        // Newly seen columns will extend the reader-owned FieldSchema.
        atoms_fields_first_ = fields;
      }
    }
  }

  col_ = ColSpec{};

  // Fields start at toks[2]
  for (std::size_t i = 2; i < toks.size(); ++i) {
    const std::string_view name = toks[i];
    const int idx = static_cast<int>(i - 2);

    if (name == "id") col_.id = idx;
    else if (name == "type") col_.type = idx;
    else if (name == "mol" || name == "molecule") col_.mol = idx;
    else if (name == "xu") col_.xu = idx;
    else if (name == "yu") col_.yu = idx;
    else if (name == "zu") col_.zu = idx;
    else if (name == "x") col_.x = idx;
    else if (name == "y") col_.y = idx;
    else if (name == "z") col_.z = idx;
    else if (name == "ix") col_.ix = idx;
    else if (name == "iy") col_.iy = idx;
    else if (name == "iz") col_.iz = idx;
  }

  if (col_.id < 0) {
    throw die("ATOMS header missing required field 'id'");
  }

  // Core field presence must remain consistent with the first frame.
  // This is enforced even when allow_header_change is enabled, because optional
  // core arrays (type/mol) are only stored in the canonical (first-frame) order.
  // Allowing them to appear/disappear mid-run would introduce ambiguous semantics
  // and can lead to crashes when measures expect stable availability.
  const bool cur_has_type = (col_.type >= 0);
  const bool cur_has_mol  = (col_.mol  >= 0);
  if (first_frame_) {
    first_has_type_ = cur_has_type;
    first_has_mol_  = cur_has_mol;
  } else {
    if (cur_has_type != first_has_type_) {
      throw die(std::string("ATOMS header core field presence changed: 'type' is ") +
                (cur_has_type ? "present" : "absent") +
                " but in the first frame it was " + (first_has_type_ ? "present" : "absent") +
                ". This is not allowed; keep core fields stable across frames.");
    }
    if (cur_has_mol != first_has_mol_) {
      throw die(std::string("ATOMS header core field presence changed: 'mol'/'molecule' is ") +
                (cur_has_mol ? "present" : "absent") +
                " but in the first frame it was " + (first_has_mol_ ? "present" : "absent") +
                ". This is not allowed; keep core fields stable across frames.");
    }
  }

  // Build extra field column mapping.
  // We automatically capture ALL non-core numeric columns present in the ATOMS header.
  // Only columns explicitly present in the dump header will be available via Frame.extra.
  extra_specs_.clear();
  extra_slot_by_col_.clear();

  const int ncols = static_cast<int>(toks.size() - 2);
  extra_slot_by_col_.assign(static_cast<std::size_t>(ncols), -1);

  auto is_core = [&](std::string_view name) {
    return (name == "id" || name == "type" || name == "mol" || name == "molecule" ||
            name == "x"  || name == "y"  || name == "z" ||
            name == "xu" || name == "yu" || name == "zu" ||
            name == "ix" || name == "iy" || name == "iz");
  };

  extra_specs_.reserve(toks.size() > 2 ? (toks.size() - 2) : 0);

  // Insert all non-core columns as extra scalar fields in the reader-owned schema.
  for (std::size_t i = 2; i < toks.size(); ++i) {
    const std::string_view name = toks[i];
    if (is_core(name)) continue;
    const int c = static_cast<int>(i - 2);
    const std::string sname(name);
    const std::size_t fid = extra_schema_.ensure(sname, FieldSchema::DType::F64, 1);
    extra_specs_.push_back(ExtraSpec{sname, fid, c});
  }

  // Build quick lookup from column index -> extra slot index.
  for (std::size_t slot = 0; slot < extra_specs_.size(); ++slot) {
    const int c = extra_specs_[slot].col;
    if (c < 0) {
      throw die("internal error: negative extra field column for '" + extra_specs_[slot].name + "'");
    }
    if (c >= ncols) {
      throw die("internal error: extra field column index out of range for '" + extra_specs_[slot].name + "'");
    }
    if (extra_slot_by_col_[static_cast<std::size_t>(c)] != -1) {
      throw die("duplicate column mapping for extra field at column " + std::to_string(c));
    }
    extra_slot_by_col_[static_cast<std::size_t>(c)] = static_cast<int>(slot);
  }
}

void LammpsDumpReader::build_id_map_() {
  if (ids_first_.empty()) {
    throw die("internal error: cannot build id map from empty first frame");
  }

  std::int64_t max_id = 0;
  for (auto id : ids_first_) {
    if (id < 0) throw die("negative atom id in first frame");
    if (id > max_id) max_id = id;
  }

  const double sparsity = static_cast<double>(max_id + 1) / static_cast<double>(natoms_);
  use_dense_map_ = (sparsity <= 10.0) && (max_id <= static_cast<std::int64_t>(std::numeric_limits<int>::max() - 1));

  if (use_dense_map_) {
    id2idx_.assign(static_cast<std::size_t>(max_id) + 1, -1);
    for (std::size_t i = 0; i < ids_first_.size(); ++i) {
      const std::int64_t id = ids_first_[i];
      if (id2idx_[static_cast<std::size_t>(id)] != -1) {
        throw die("duplicate atom id in first frame: " + std::to_string(id));
      }
      id2idx_[static_cast<std::size_t>(id)] = static_cast<int>(i);
    }
  } else {
    id2idx_hash_.reserve(ids_first_.size() * 2);
    for (std::size_t i = 0; i < ids_first_.size(); ++i) {
      const std::int64_t id = ids_first_[i];
      auto [it, ok] = id2idx_hash_.emplace(id, i);
      if (!ok) {
        throw die("duplicate atom id in first frame (hash): " + std::to_string(id));
      }
    }
  }
}

std::size_t LammpsDumpReader::lookup_idx_(std::int64_t id) const {
  if (use_dense_map_) {
    if (id < 0 || static_cast<std::size_t>(id) >= id2idx_.size()) {
      throw die("atom id out of dense map range: " + std::to_string(id));
    }
    const int idx = id2idx_[static_cast<std::size_t>(id)];
    if (idx < 0) {
      throw die("atom id not found in first frame: " + std::to_string(id));
    }
    return static_cast<std::size_t>(idx);
  }

  auto it = id2idx_hash_.find(id);
  if (it == id2idx_hash_.end()) {
    throw die("atom id not found in first frame (hash): " + std::to_string(id));
  }
  return it->second;
}

} // namespace pilots
