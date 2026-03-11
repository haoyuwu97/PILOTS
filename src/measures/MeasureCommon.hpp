#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <optional>
#include <numbers>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "pilots/alg/index/ChainIndex.hpp"
#include "pilots/config/IniConfig.hpp"
#include "pilots/core/Frame.hpp"
#include "pilots/core/SystemContext.hpp"
#include "pilots/correlate/CorrelatorSpec.hpp"
#include "pilots/correlate/LagAxis.hpp"
#include "pilots/select/SelectionProvider.hpp"
#include "pilots/select/SelectionView.hpp"
#include "pilots/topology/Topology.hpp"

namespace pilots::measure_ext {
namespace fs = std::filesystem;

inline fs::path resolve_path(const fs::path& base_dir, const std::string& p) {
  fs::path path(p);
  if (path.is_absolute()) return path;
  return (base_dir / path).lexically_normal();
}

template <class EnvLike>
inline fs::path resolve_measure_output_dir(const IniConfig& cfg,
                                           const std::string& section,
                                           const EnvLike& env) {
  const std::string outdir_s = cfg.get_string(section, "output_dir", std::optional<std::string>(""));
  const fs::path output_dir = outdir_s.empty() ? env.output_dir_general : resolve_path(env.cfg_dir, outdir_s);
  if (!env.dry_run) fs::create_directories(output_dir);
  return output_dir;
}

template <class EnvLike>
inline fs::path resolve_measure_output_path(const IniConfig& cfg,
                                            const std::string& section,
                                            const EnvLike& env,
                                            const std::string& key,
                                            const std::string& default_name) {
  return (resolve_measure_output_dir(cfg, section, env)
          / cfg.get_string(section, key, std::optional<std::string>(default_name))).lexically_normal();
}

inline std::string dstr(double v) {
  std::ostringstream oss;
  oss.precision(17);
  oss << v;
  return oss.str();
}

inline std::string x_unit_for_axis(LagAxis a) {
  switch (a) {
    case LagAxis::Frame: return "frames";
    case LagAxis::Timestep: return "timesteps";
    case LagAxis::TimeBin: return "simulation_time";
  }
  return "frames";
}

inline int parse_diag_mask(const std::string& s_raw) {
  std::string s = s_raw;
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  if (s == "xxyyzz" || s == "all" || s == "xyz") return 7;
  if (s == "xx" || s == "x") return 1;
  if (s == "yy" || s == "y") return 2;
  if (s == "zz" || s == "z") return 4;
  if (s == "xxyy" || s == "xy") return 3;
  if (s == "xxzz" || s == "xz") return 5;
  if (s == "yyzz" || s == "yz") return 6;
  throw std::runtime_error(
      "invalid components string: '" + s_raw + "' "
      "(use xxyyzz, xx, yy, zz, xxyy, xxzz, yyzz)");
}

inline int count_diag_dims(int diag_mask) {
  int d = 0;
  if ((diag_mask & 1) != 0) ++d;
  if ((diag_mask & 2) != 0) ++d;
  if ((diag_mask & 4) != 0) ++d;
  return d;
}

inline std::size_t count_lammps_frames(const std::string& path) {
  std::ifstream ifs(path);
  if (!ifs) {
    throw std::runtime_error("failed to open input for frame count: " + path);
  }
  std::string line;
  std::size_t n = 0;
  while (std::getline(ifs, line)) {
    if (line.rfind("ITEM: TIMESTEP", 0) == 0) ++n;
  }
  return n;
}

inline std::int64_t resolve_exact_frame_end(const CorrelatorSpec& corr,
                                            bool follow,
                                            std::int64_t frame_start,
                                            std::int64_t frame_end,
                                            const fs::path& input_path) {
  if (corr.type != "exact") return frame_end;
  if (follow) {
    throw std::runtime_error(
        "correlator=exact is not supported in follow mode (frame_end is unbounded). "
        "Use correlator=multitau.");
  }
  if (frame_end >= 0) return frame_end;
  const std::size_t total_frames = count_lammps_frames(input_path.string());
  if (total_frames == 0) throw std::runtime_error("input contains 0 frames");
  const std::int64_t eff = static_cast<std::int64_t>(total_frames - 1);
  if (eff < frame_start) {
    throw std::runtime_error("effective frame_end < frame_start");
  }
  return eff;
}

inline std::vector<double> parse_double_list(const IniConfig& cfg,
                                             const std::string& section,
                                             const std::string& key) {
  const auto raw = cfg.get_list(section, key);
  std::vector<double> out;
  out.reserve(raw.size());
  for (const auto& s : raw) {
    try {
      std::size_t pos = 0;
      const double v = std::stod(s, &pos);
      if (pos != s.size()) throw std::invalid_argument("trailing");
      out.push_back(v);
    } catch (...) {
      throw std::runtime_error("failed to parse double list item for " + section + "." + key + ": '" + s + "'");
    }
  }
  if (out.empty()) {
    throw std::runtime_error("empty list for required key " + section + "." + key);
  }
  return out;
}

inline std::vector<int> parse_int_list(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& key) {
  const auto raw = cfg.get_list(section, key);
  std::vector<int> out;
  out.reserve(raw.size());
  for (const auto& s : raw) {
    try {
      std::size_t pos = 0;
      const long v = std::stol(s, &pos);
      if (pos != s.size()) throw std::invalid_argument("trailing");
      out.push_back(static_cast<int>(v));
    } catch (...) {
      throw std::runtime_error("failed to parse int list item for " + section + "." + key + ": '" + s + "'");
    }
  }
  if (out.empty()) {
    throw std::runtime_error("empty list for required key " + section + "." + key);
  }
  return out;
}

inline std::string join_doubles(const std::vector<double>& vals) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < vals.size(); ++i) {
    if (i) oss << ",";
    oss << dstr(vals[i]);
  }
  return oss.str();
}

inline std::string join_ints(const std::vector<int>& vals) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < vals.size(); ++i) {
    if (i) oss << ",";
    oss << vals[i];
  }
  return oss.str();
}

inline SelectionView get_static_combined_view(SelectionProvider& sp,
                                              const Frame& frame0,
                                              const std::string& group_ref,
                                              const std::string& topo_group_ref,
                                              const std::string& combine_expr,
                                              const std::string& measure_name) {
  if (sp.is_dynamic_spec(group_ref, topo_group_ref)) {
    throw std::runtime_error(measure_name + " requires a static selection; the selection spec depends on a dynamic group/topo_group");
  }
  return sp.get_combined_view(frame0, 0, group_ref, topo_group_ref, combine_expr);
}

inline SelectionView get_static_group_view(SelectionProvider& sp,
                                           const Frame& frame0,
                                           const std::string& group_ref,
                                           const std::string& measure_name) {
  if (sp.is_dynamic_spec(group_ref, "all")) {
    throw std::runtime_error(measure_name + " requires a static drift_group; drift_group depends on a dynamic group");
  }
  return sp.get_combined_view(frame0, 0, group_ref, "all", "A");
}

inline std::vector<std::int64_t> dfield_to_i64(std::span<const double> src,
                                               const std::string& field_name,
                                               bool strict_integer = true) {
  std::vector<std::int64_t> out(src.size(), 0);
  for (std::size_t i = 0; i < src.size(); ++i) {
    const double x = src[i];
    const std::int64_t v = static_cast<std::int64_t>(x);
    if (strict_integer && static_cast<double>(v) != x) {
      throw std::runtime_error("field '" + field_name + "' contains a non-integer value; this measure requires integer-like identifiers");
    }
    out[i] = v;
  }
  return out;
}

template <class CapsLike>
inline void append_integer_like_field_cap(CapsLike& caps,
                                          const std::string& field) {
  if (field == "id" || field == "mol") {
    caps.requires_i64fields.push_back(field);
  } else if (field == "type") {
    caps.requires_intfields.push_back(field);
  } else {
    caps.requires_dfields.push_back(field);
  }
}

inline std::vector<std::int64_t> integer_like_field_to_i64(const Frame& frame,
                                                           const std::string& field,
                                                           bool strict_integer = true);

inline std::vector<std::int64_t> chain_id_per_atom_from_config(const IniConfig& cfg,
                                                               const std::string& section,
                                                               const Frame& frame0,
                                                               const SystemContext& sysctx) {
  if (cfg.has_key(section, "chain_id_field")) {
    const std::string field = cfg.get_string(section, "chain_id_field");
    return integer_like_field_to_i64(frame0, field, true);
  }

  if (frame0.has_mol) {
    const auto mol = frame0.require_i64field("mol");
    return std::vector<std::int64_t>(mol.begin(), mol.end());
  }

  if (sysctx.topology && sysctx.topology->has_bonds()) {
    return sysctx.topology->derived_mol_ids_from_bonds(frame0);
  }

  throw std::runtime_error(
      "cannot determine chain membership: provide chain_id_field, include 'mol' in the dump, or load topology bonds");
}

struct SelectedChains {
  std::vector<std::int64_t> chain_id;
  std::vector<std::vector<std::size_t>> atom_by_chain;   // canonical atom indices
  std::vector<std::vector<std::size_t>> selpos_by_chain; // selection-local positions
  std::vector<std::size_t> chain_of_selpos;              // size = nsel

  std::size_t n_chains() const { return chain_id.size(); }
};

inline SelectedChains build_selected_chains(const SelectionView& sel,
                                            std::span<const std::int64_t> chain_id_per_atom) {
  if (sel.idx.empty()) {
    throw std::runtime_error("build_selected_chains: selection is empty");
  }

  std::unordered_map<std::int64_t, std::vector<std::size_t>> tmp_selpos;
  std::unordered_map<std::int64_t, std::vector<std::size_t>> tmp_atoms;
  tmp_selpos.reserve(sel.idx.size() / 4 + 8);
  tmp_atoms.reserve(sel.idx.size() / 4 + 8);

  for (std::size_t p = 0; p < sel.idx.size(); ++p) {
    const std::size_t i = sel.idx[p];
    if (i >= chain_id_per_atom.size()) {
      throw std::runtime_error("build_selected_chains: selection index out of range for chain_id vector");
    }
    const std::int64_t cid = chain_id_per_atom[i];
    tmp_selpos[cid].push_back(p);
    tmp_atoms[cid].push_back(i);
  }

  std::vector<std::int64_t> cids;
  cids.reserve(tmp_selpos.size());
  for (const auto& kv : tmp_selpos) cids.push_back(kv.first);
  std::sort(cids.begin(), cids.end());

  SelectedChains out;
  out.chain_id = cids;
  out.atom_by_chain.resize(cids.size());
  out.selpos_by_chain.resize(cids.size());
  out.chain_of_selpos.assign(sel.idx.size(), static_cast<std::size_t>(-1));

  for (std::size_t c = 0; c < cids.size(); ++c) {
    out.atom_by_chain[c] = std::move(tmp_atoms[cids[c]]);
    out.selpos_by_chain[c] = std::move(tmp_selpos[cids[c]]);
    for (const std::size_t p : out.selpos_by_chain[c]) {
      out.chain_of_selpos[p] = c;
    }
  }

  return out;
}

inline std::vector<std::size_t> order_linear_chain_from_topology(const Topology& topo,
                                                                 const std::vector<std::size_t>& atoms,
                                                                 bool strict = true) {
  if (!topo.has_bonds()) {
    throw std::runtime_error("order_linear_chain_from_topology: topology bonds are not loaded");
  }
  if (atoms.empty()) return {};
  if (atoms.size() == 1) return atoms;

  std::unordered_set<std::size_t> atom_set(atoms.begin(), atoms.end());
  std::unordered_map<std::size_t, std::vector<std::size_t>> local_adj;
  local_adj.reserve(atoms.size() * 2);

  std::size_t edge_count = 0;
  for (const std::size_t i : atoms) {
    if (i >= topo.adjacency.size()) {
      throw std::runtime_error("order_linear_chain_from_topology: atom index out of range for topology adjacency");
    }
    auto& nbrs = local_adj[i];
    for (const std::size_t j : topo.adjacency[i]) {
      if (atom_set.find(j) != atom_set.end()) {
        nbrs.push_back(j);
      }
    }
    std::sort(nbrs.begin(), nbrs.end());
    edge_count += nbrs.size();
  }
  edge_count /= 2;

  std::vector<std::size_t> ends;
  ends.reserve(2);
  for (const std::size_t i : atoms) {
    const std::size_t deg = local_adj[i].size();
    if (deg == 1) {
      ends.push_back(i);
    } else if (deg == 2) {
      // ok
    } else {
      if (strict) {
        throw std::runtime_error("topology-restricted chain is not linear (degree != 1/2)");
      }
      return {};
    }
  }

  if (edge_count != atoms.size() - 1 || ends.size() != 2) {
    if (strict) {
      throw std::runtime_error("topology-restricted chain is not a single linear path");
    }
    return {};
  }

  const std::size_t start = std::min(ends[0], ends[1]);
  std::vector<std::size_t> ordered;
  ordered.reserve(atoms.size());
  std::size_t prev = static_cast<std::size_t>(-1);
  std::size_t cur = start;

  for (;;) {
    ordered.push_back(cur);
    const auto& nbrs = local_adj[cur];
    std::size_t next = static_cast<std::size_t>(-1);
    for (const std::size_t nb : nbrs) {
      if (nb != prev) {
        next = nb;
        break;
      }
    }
    if (next == static_cast<std::size_t>(-1)) break;
    prev = cur;
    cur = next;
    if (ordered.size() > atoms.size()) {
      throw std::runtime_error("topology-restricted chain walk exceeded chain size");
    }
  }

  if (ordered.size() != atoms.size()) {
    if (strict) {
      throw std::runtime_error("topology-restricted chain walk did not visit all nodes");
    }
    return {};
  }

  return ordered;
}

struct OrderedChains {
  std::vector<std::int64_t> chain_id;
  std::vector<std::vector<std::size_t>> atom_by_chain;   // canonical atom indices in chain order
  std::vector<std::vector<std::size_t>> selpos_by_chain; // selection-local positions in chain order

  std::size_t n_chains() const { return chain_id.size(); }
};

inline OrderedChains build_ordered_chains_from_chain_pos(const SelectedChains& groups,
                                                         std::span<const double> chain_pos_per_atom,
                                                         bool strict = true) {
  OrderedChains out;
  out.chain_id.reserve(groups.n_chains());
  out.atom_by_chain.reserve(groups.n_chains());
  out.selpos_by_chain.reserve(groups.n_chains());

  for (std::size_t c = 0; c < groups.n_chains(); ++c) {
    std::vector<std::pair<double, std::size_t>> tmp;
    tmp.reserve(groups.atom_by_chain[c].size());
    for (std::size_t k = 0; k < groups.atom_by_chain[c].size(); ++k) {
      const std::size_t atom = groups.atom_by_chain[c][k];
      if (atom >= chain_pos_per_atom.size()) {
        throw std::runtime_error("build_ordered_chains_from_chain_pos: chain_pos vector is too small");
      }
      tmp.push_back({chain_pos_per_atom[atom], k});
    }

    std::sort(tmp.begin(), tmp.end(), [](const auto& a, const auto& b) {
      if (a.first != b.first) return a.first < b.first;
      return a.second < b.second;
    });

    if (strict) {
      for (std::size_t k = 1; k < tmp.size(); ++k) {
        if (tmp[k].first == tmp[k - 1].first) {
          throw std::runtime_error("build_ordered_chains_from_chain_pos: duplicate chain_pos inside one chain");
        }
      }
    }

    out.chain_id.push_back(groups.chain_id[c]);
    out.selpos_by_chain.emplace_back();
    out.atom_by_chain.emplace_back();
    out.selpos_by_chain.back().reserve(tmp.size());
    out.atom_by_chain.back().reserve(tmp.size());
    for (const auto& kv : tmp) {
      const std::size_t k = kv.second;
      out.selpos_by_chain.back().push_back(groups.selpos_by_chain[c][k]);
      out.atom_by_chain.back().push_back(groups.atom_by_chain[c][k]);
    }
  }
  return out;
}

inline OrderedChains build_ordered_chains_from_topology(const SelectedChains& groups,
                                                        const SelectionView& sel,
                                                        const Topology& topo,
                                                        bool strict = true) {
  std::unordered_map<std::size_t, std::size_t> atom_to_selpos;
  atom_to_selpos.reserve(sel.idx.size() * 2);
  for (std::size_t p = 0; p < sel.idx.size(); ++p) {
    atom_to_selpos.emplace(sel.idx[p], p);
  }

  OrderedChains out;
  out.chain_id.reserve(groups.n_chains());
  out.atom_by_chain.reserve(groups.n_chains());
  out.selpos_by_chain.reserve(groups.n_chains());

  for (std::size_t c = 0; c < groups.n_chains(); ++c) {
    const auto ordered_atoms = order_linear_chain_from_topology(topo, groups.atom_by_chain[c], strict);
    if (ordered_atoms.empty()) continue;

    out.chain_id.push_back(groups.chain_id[c]);
    out.atom_by_chain.push_back(ordered_atoms);
    auto& selpos = out.selpos_by_chain.emplace_back();
    selpos.reserve(ordered_atoms.size());
    for (const std::size_t i : ordered_atoms) {
      auto it = atom_to_selpos.find(i);
      if (it == atom_to_selpos.end()) {
        throw std::runtime_error("build_ordered_chains_from_topology: ordered atom is not in the selection");
      }
      selpos.push_back(it->second);
    }
  }

  return out;
}

inline OrderedChains ordered_chains_from_config(const IniConfig& cfg,
                                                const std::string& section,
                                                const Frame& frame0,
                                                const SystemContext& sysctx,
                                                const SelectionView& sel,
                                                const SelectedChains& groups) {
  const bool strict_linear = cfg.get_bool(section, "strict_linear", std::optional<bool>(true));

  if (cfg.has_key(section, "chain_pos_field")) {
    const std::string field = cfg.get_string(section, "chain_pos_field");
    const auto pos = frame0.require_dfield(field);
    return build_ordered_chains_from_chain_pos(groups, pos, strict_linear);
  }

  if (!sysctx.topology || !sysctx.topology->has_bonds()) {
    throw std::runtime_error(
        "ordered chain analysis requires chain_pos_field or topology bonds for linear-chain ordering");
  }
  return build_ordered_chains_from_topology(groups, sel, *sysctx.topology, strict_linear);
}


enum class Axis1D { X = 0, Y = 1, Z = 2 };

inline Axis1D parse_axis1d(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  if (s == "x") return Axis1D::X;
  if (s == "y") return Axis1D::Y;
  if (s == "z") return Axis1D::Z;
  throw std::runtime_error("invalid axis spec: '" + s + "' (use x|y|z)");
}

inline const char* axis1d_name(Axis1D a) {
  switch (a) {
    case Axis1D::X: return "x";
    case Axis1D::Y: return "y";
    case Axis1D::Z: return "z";
  }
  return "x";
}

inline double axis_length(const Box& box, Axis1D a) {
  switch (a) {
    case Axis1D::X: return box.lx();
    case Axis1D::Y: return box.ly();
    case Axis1D::Z: return box.lz();
  }
  return box.lx();
}

inline double box_volume(const Box& box) {
  return box.lx() * box.ly() * box.lz();
}

inline double orth_area_for_axis(const Box& box, Axis1D a) {
  switch (a) {
    case Axis1D::X: return box.ly() * box.lz();
    case Axis1D::Y: return box.lx() * box.lz();
    case Axis1D::Z: return box.lx() * box.ly();
  }
  return box.ly() * box.lz();
}

inline double primary_axis_coord(const Box& box, double x, double y, double z, Axis1D a) {
  auto lam = box.to_lambda(x, y, z);
  switch (a) {
    case Axis1D::X:
      lam[0] -= std::floor(lam[0]);
      return lam[0] * box.lx();
    case Axis1D::Y:
      lam[1] -= std::floor(lam[1]);
      return lam[1] * box.ly();
    case Axis1D::Z:
      lam[2] -= std::floor(lam[2]);
      return lam[2] * box.lz();
  }
  lam[0] -= std::floor(lam[0]);
  return lam[0] * box.lx();
}

inline double shell_volume_3d(double rlo, double rhi) {
  return (4.0 * std::numbers::pi / 3.0) * (rhi * rhi * rhi - rlo * rlo * rlo);
}

inline std::vector<std::int64_t> integer_like_field_to_i64(const Frame& frame,
                                                           const std::string& field,
                                                           bool strict_integer) {
  if (field == "id" || field == "mol") {
    const auto v = frame.require_i64field(field);
    return std::vector<std::int64_t>(v.begin(), v.end());
  }
  if (field == "type") {
    const auto v = frame.require_intfield(field);
    std::vector<std::int64_t> out(v.size(), 0);
    for (std::size_t i = 0; i < v.size(); ++i) out[i] = static_cast<std::int64_t>(v[i]);
    return out;
  }
  return dfield_to_i64(frame.require_dfield(field), field, strict_integer);
}

inline std::vector<std::int64_t> entity_id_per_atom_from_config(const IniConfig& cfg,
                                                                const std::string& section,
                                                                const Frame& frame0,
                                                                const SystemContext& sysctx) {
  if (cfg.has_key(section, "entity_id_field")) {
    return integer_like_field_to_i64(frame0, cfg.get_string(section, "entity_id_field"), true);
  }
  if (frame0.has_mol) {
    const auto mol = frame0.require_i64field("mol");
    return std::vector<std::int64_t>(mol.begin(), mol.end());
  }
  if (sysctx.topology && sysctx.topology->has_bonds()) {
    return sysctx.topology->derived_mol_ids_from_bonds(frame0);
  }
  throw std::runtime_error(
      "cannot determine entity membership: provide entity_id_field, include 'mol' in the dump, or load topology bonds");
}

inline std::vector<double> mass_by_atom_from_config(const IniConfig& cfg,
                                                    const std::string& section,
                                                    const Frame& frame0,
                                                    const SystemContext& sysctx) {
  if (cfg.has_key(section, "mass_field")) {
    const auto m = frame0.require_dfield(cfg.get_string(section, "mass_field"));
    return std::vector<double>(m.begin(), m.end());
  }
  if (!frame0.has_type) {
    throw std::runtime_error("mass_by_atom_from_config: dump is missing 'type' and no mass_field was provided");
  }
  if (!sysctx.topology || !sysctx.topology->has_masses()) {
    throw std::runtime_error("mass_by_atom_from_config: topology masses are required when mass_field is not provided");
  }
  const auto type = frame0.require_intfield("type");
  std::vector<double> out(frame0.natoms, 0.0);
  const auto& mbt = sysctx.topology->mass_by_type;
  for (std::size_t i = 0; i < frame0.natoms; ++i) {
    const int t = type[i];
    if (t <= 0 || static_cast<std::size_t>(t) >= mbt.size()) {
      throw std::runtime_error("mass_by_atom_from_config: atom type out of range for topology Masses table");
    }
    out[i] = mbt[static_cast<std::size_t>(t)];
  }
  return out;
}

inline bool same_index_set(std::span<const std::size_t> a, std::span<const std::size_t> b) {
  if (a.size() != b.size()) return false;
  if (std::equal(a.begin(), a.end(), b.begin(), b.end())) return true;
  std::vector<std::size_t> sa(a.begin(), a.end());
  std::vector<std::size_t> sb(b.begin(), b.end());
  std::sort(sa.begin(), sa.end());
  std::sort(sb.begin(), sb.end());
  return sa == sb;
}

inline bool same_index_set(const SelectionView& a, const SelectionView& b) {
  return same_index_set(a.idx, b.idx);
}

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3{a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(double s, const Vec3& a) { return Vec3{s * a.x, s * a.y, s * a.z}; }
inline Vec3 operator*(const Vec3& a, double s) { return s * a; }
inline Vec3 operator/(const Vec3& a, double s) { return Vec3{a.x / s, a.y / s, a.z / s}; }
inline Vec3& operator+=(Vec3& a, const Vec3& b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm2(const Vec3& a) { return dot(a, a); }
inline double norm(const Vec3& a) { return std::sqrt(norm2(a)); }
inline Vec3 safe_unit(const Vec3& a, double eps = 1e-18) {
  const double n = norm(a);
  if (n <= eps) return Vec3{0.0, 0.0, 0.0};
  return a / n;
}

inline Vec3 min_image_difference(const pilots::Box& box,
                                 const Vec3& a,
                                 const Vec3& b) {
  const auto d = box.min_image_displacement(b.x, b.y, b.z, a.x, a.y, a.z);
  return Vec3{d[0], d[1], d[2]};
}

inline Vec3 atom_vec3(const Frame& frame, std::size_t i) {
  if (i >= frame.natoms) throw std::runtime_error("atom_vec3: atom index out of range");
  return Vec3{frame.xu[i], frame.yu[i], frame.zu[i]};
}

inline std::int64_t resolve_finite_frame_end(bool follow,
                                             std::int64_t frame_start,
                                             std::int64_t frame_end,
                                             const fs::path& input_path,
                                             const std::string& who) {
  if (follow) {
    throw std::runtime_error(who + ": analysis requires a finite frame_end and is not supported in follow mode");
  }
  if (frame_end >= 0) {
    if (frame_end < frame_start) throw std::runtime_error(who + ": frame_end must be >= frame_start");
    return frame_end;
  }
  const std::size_t total_frames = count_lammps_frames(input_path.string());
  if (total_frames == 0) throw std::runtime_error(who + ": input contains 0 frames");
  const std::int64_t eff = static_cast<std::int64_t>(total_frames - 1);
  if (eff < frame_start) throw std::runtime_error(who + ": effective frame_end < frame_start");
  return eff;
}



inline double wrap_unit_interval(double x) {
  x -= std::floor(x);
  if (x < 0.0) x += 1.0;
  return x;
}

inline double axis_fraction_from_lambda(const std::array<double, 3>& lam,
                                        Axis1D a) {
  switch (a) {
    case Axis1D::X: return wrap_unit_interval(lam[0]);
    case Axis1D::Y: return wrap_unit_interval(lam[1]);
    case Axis1D::Z: return wrap_unit_interval(lam[2]);
  }
  return wrap_unit_interval(lam[0]);
}

inline double axis_fraction_from_xyz(const Box& box,
                                     double x,
                                     double y,
                                     double z,
                                     Axis1D a) {
  return axis_fraction_from_lambda(box.to_lambda(x, y, z), a);
}

inline double circular_mean_axis_fraction(const Frame& frame,
                                          const SelectionView& sel,
                                          Axis1D axis) {
  if (sel.idx.empty()) {
    throw std::runtime_error("circular_mean_axis_fraction: anchor selection is empty");
  }
  const auto xu = frame.require_dfield("xu");
  const auto yu = frame.require_dfield("yu");
  const auto zu = frame.require_dfield("zu");
  double csum = 0.0;
  double ssum = 0.0;
  for (const std::size_t i : sel.idx) {
    const double f = axis_fraction_from_xyz(frame.box, xu[i], yu[i], zu[i], axis);
    const double ang = 2.0 * std::numbers::pi * f;
    csum += std::cos(ang);
    ssum += std::sin(ang);
  }
  if (std::abs(csum) < 1e-18 && std::abs(ssum) < 1e-18) {
    const std::size_t i0 = sel.idx.front();
    return axis_fraction_from_xyz(frame.box, xu[i0], yu[i0], zu[i0], axis);
  }
  double ang = std::atan2(ssum, csum) / (2.0 * std::numbers::pi);
  if (ang < 0.0) ang += 1.0;
  return ang;
}

struct SlabTransformOptions {
  Axis1D axis = Axis1D::Z;
  bool slab_align_recenter = false;
  bool halfcell_fold = false;
  double target_center_frac = 0.5;
};

struct SlabFrameTransform {
  Axis1D axis = Axis1D::Z;
  double box_length = 0.0;
  double shift_frac = 0.0;
  double center_frac = 0.5;
  bool slab_align_recenter = false;
  bool halfcell_fold = false;

  double domain_length() const {
    return halfcell_fold ? (0.5 * box_length) : box_length;
  }

  double raw_to_aligned_fraction(double raw_frac) const {
    return wrap_unit_interval(raw_frac + shift_frac);
  }

  double aligned_to_raw_fraction(double aligned_frac) const {
    return wrap_unit_interval(aligned_frac - shift_frac);
  }

  double fold_aligned_coord(double aligned_coord) const {
    const double center = center_frac * box_length;
    double d = aligned_coord - center;
    d -= std::nearbyint(d / box_length) * box_length;
    return 0.5 * box_length - std::abs(d);
  }

  double point_to_domain_coord(const Box& box,
                               double x,
                               double y,
                               double z) const {
    const double raw_frac = axis_fraction_from_xyz(box, x, y, z, axis);
    const double aligned_frac = raw_to_aligned_fraction(raw_frac);
    const double aligned_coord = aligned_frac * box_length;
    return halfcell_fold ? fold_aligned_coord(aligned_coord) : aligned_coord;
  }

  double point_to_domain_fraction(const Box& box,
                                  double x,
                                  double y,
                                  double z) const {
    const double D = domain_length();
    if (D <= 0.0) return 0.0;
    return point_to_domain_coord(box, x, y, z) / D;
  }

  double aligned_fraction_from_domain_fraction(double domain_frac,
                                               int branch = 0) const {
    const double u = std::clamp(domain_frac, 0.0, 1.0);
    if (!halfcell_fold) return u;
    const double dist_frac = 0.5 * (1.0 - u);
    const double sign = (branch % 2 == 0) ? 1.0 : -1.0;
    return wrap_unit_interval(center_frac + sign * dist_frac);
  }

  std::array<double, 3> lambda_for_domain_sample(double domain_frac,
                                                 double transverse_u,
                                                 double transverse_v,
                                                 int branch = 0) const {
    std::array<double, 3> lam{0.0, 0.0, 0.0};
    const double axis_frac = aligned_to_raw_fraction(aligned_fraction_from_domain_fraction(domain_frac, branch));
    switch (axis) {
      case Axis1D::X:
        lam = {axis_frac, wrap_unit_interval(transverse_u), wrap_unit_interval(transverse_v)};
        break;
      case Axis1D::Y:
        lam = {wrap_unit_interval(transverse_u), axis_frac, wrap_unit_interval(transverse_v)};
        break;
      case Axis1D::Z:
        lam = {wrap_unit_interval(transverse_u), wrap_unit_interval(transverse_v), axis_frac};
        break;
    }
    return lam;
  }
};

inline SlabFrameTransform make_slab_frame_transform(const Frame& frame,
                                                    const SelectionView* anchor_sel,
                                                    const SlabTransformOptions& opt) {
  SlabFrameTransform tf;
  tf.axis = opt.axis;
  tf.box_length = axis_length(frame.box, opt.axis);
  tf.center_frac = opt.target_center_frac;
  tf.slab_align_recenter = opt.slab_align_recenter;
  tf.halfcell_fold = opt.halfcell_fold;
  if (opt.slab_align_recenter && anchor_sel != nullptr) {
    const double anchor_center = circular_mean_axis_fraction(frame, *anchor_sel, opt.axis);
    tf.shift_frac = wrap_unit_interval(opt.target_center_frac - anchor_center);
  } else {
    tf.shift_frac = 0.0;
  }
  return tf;
}

inline bool slab_transform_requested(const IniConfig& cfg,
                                     const std::string& section) {
  return cfg.get_bool(section, "slab_align_recenter", std::optional<bool>(false))
      || cfg.get_bool(section, "halfcell_fold", std::optional<bool>(false))
      || cfg.has_key(section, "anchor_group")
      || cfg.has_key(section, "anchor_topo_group")
      || cfg.has_key(section, "anchor_combine");
}

template <class ProviderLike>
inline std::optional<SelectionView> build_optional_anchor_selection(const IniConfig& cfg,
                                                                    const std::string& section,
                                                                    ProviderLike& sp,
                                                                    const Frame& frame0,
                                                                    const SelectionView& fallback_sel,
                                                                    bool needed,
                                                                    const std::string& measure_name) {
  const bool has_anchor_keys = cfg.has_key(section, "anchor_group")
                            || cfg.has_key(section, "anchor_topo_group")
                            || cfg.has_key(section, "anchor_combine");
  if (!needed && !has_anchor_keys) return std::nullopt;
  if (!has_anchor_keys) return fallback_sel;
  const std::string group = cfg.get_string(section, "anchor_group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "anchor_topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "anchor_combine", std::optional<std::string>("A&T"));
  return get_static_combined_view(sp, frame0, group, topo, comb, measure_name + " anchor");
}

inline double transformed_fraction_for_axis(const SlabFrameTransform* tf,
                                            const Box& box,
                                            double x,
                                            double y,
                                            double z,
                                            Axis1D axis) {
  if (tf != nullptr && axis == tf->axis) {
    return tf->point_to_domain_fraction(box, x, y, z);
  }
  return axis_fraction_from_xyz(box, x, y, z, axis);
}

inline double transformed_axis_length(const SlabFrameTransform* tf,
                                      const Box& box,
                                      Axis1D axis) {
  if (tf != nullptr && axis == tf->axis) return tf->domain_length();
  return axis_length(box, axis);
}
} // namespace pilots::measure_ext
