#pragma once

#include <cstddef>
#include <cstdint>
#include <csignal>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "pilots/core/Frame.hpp"

namespace pilots {

class LammpsDumpReader {
public:
  // Reader automatically captures ALL per-atom columns present in the dump "ITEM: ATOMS" header.
  // - Core fields (id/type/mol and coordinates) are stored in dedicated Frame arrays.
  // - Any additional numeric scalar columns are stored as Frame extra fields using a reader-owned
  //   FieldSchema (name -> field_id), and written into Frame::extra_values[field_id] (SoA).
  // NOTE: Only columns explicitly present in the dump header are available.
  // Measures must explicitly require the fields they need; missing fields should trigger an error.
  explicit LammpsDumpReader(const std::string& path);

  // Follow mode: if enabled, next(frame) will block and wait for additional data when reaching EOF.
  // Useful for reading a dump file while LAMMPS is still writing it.
  void set_follow(bool enable, int poll_ms = 200);

  // By default, PILOTS requires a stable "ITEM: ATOMS" header across frames.
  // If allow_header_change is enabled, the reader will accept header changes by
  // rebuilding the per-frame column mapping (and expanding the internal FieldSchema
  // for new columns). Missing fields required by measures will still fail-fast.
  // NOTE: even with allow_header_change enabled, the presence of core identity
  // fields (type/mol) must remain consistent with the first frame. This avoids
  // ambiguous semantics and prevents crashes when optional core arrays would
  // otherwise appear/disappear mid-run.
  void set_allow_header_change(bool enable);

  // Optional: stop flag for follow mode. If provided and becomes non-zero, next() will return false.
  void set_stop_flag(const volatile std::sig_atomic_t* flag);

  // Read next frame into `frame`. On the first frame, also populates:
  // frame.id, frame.type, frame.mol (canonical ordering = first frame file order)
  // Returns false on EOF.
  bool next(Frame& frame);

  std::size_t natoms() const { return natoms_; }

  // P0: checkpoint/resume support.
  // Returns the byte offset (from file begin) immediately after the last successfully read full frame.
  std::uint64_t tell_offset() const { return last_good_off_; }

  // Seek to a previously stored offset (from tell_offset). This should land on a frame boundary.
  void seek_offset(std::uint64_t off);

private:
  struct ColSpec {
    int id = -1;
    int type = -1;
    int mol = -1;

    int xu = -1, yu = -1, zu = -1;
    int x = -1, y = -1, z = -1;
    int ix = -1, iy = -1, iz = -1;

    bool have_unwrapped() const { return xu >= 0 && yu >= 0 && zu >= 0; }
    bool have_wrapped_with_images() const { return x >= 0 && y >= 0 && z >= 0 && ix >= 0 && iy >= 0 && iz >= 0; }
  };

  std::ifstream ifs_;
  bool follow_ = false;
  int poll_ms_ = 200;
  const volatile std::sig_atomic_t* stop_flag_ = nullptr;
  bool allow_header_change_ = false;
  bool first_frame_ = true;
  std::size_t natoms_ = 0;

  // Byte offset after the last successfully read full frame.
  std::uint64_t last_good_off_ = 0;

  // id -> idx mapping (Level-2 optimization). Dense vector if ids are reasonably dense, else hash.
  bool use_dense_map_ = true;
  std::vector<int> id2idx_;
  std::unordered_map<std::int64_t, std::size_t> id2idx_hash_;

  std::vector<std::int64_t> ids_first_;
  ColSpec col_;

  // Extra scalar fields: resolved to a stable field_id in the reader-owned schema.
  struct ExtraSpec { std::string name; std::size_t fid = 0; int col = -1; };
  std::vector<ExtraSpec> extra_specs_;
  std::vector<int> extra_slot_by_col_; // size = ncols, value = slot index or -1

  // Reader-owned schema for extra per-atom fields.
  FieldSchema extra_schema_;

  // Schema validation: we assume a stable ATOMS header across frames.
  // If it changes (column names/order), we fail fast to avoid silent data corruption.
  bool atoms_fields_set_ = false;
  std::vector<std::string> atoms_fields_first_; // full field list in order, excluding "ITEM: ATOMS"

  // Core field presence in the FIRST frame's ATOMS header.
  // Even when allow_header_change_ is true, these must remain consistent.
  bool first_has_type_ = false;
  bool first_has_mol_ = false;

  bool read_header_(Frame& frame);
  bool read_box_(Frame& frame, const std::string& header_line);
  void parse_atoms_header_(const std::string& atoms_header_line);
  void build_id_map_();

  std::size_t lookup_idx_(std::int64_t id) const;
};

} // namespace pilots
