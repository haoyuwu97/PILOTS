#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pilots/alg/mapping/BeadFrame.hpp"
#include "pilots/alg/mapping/MappingSpec.hpp"
#include "pilots/core/Frame.hpp"
#include "pilots/topology/Topology.hpp"
#include "pilots/util/Hash.hpp"

namespace pilots::alg::mapping {
namespace fs = std::filesystem;

struct MappingFileEntry {
  std::int64_t atom_id = 0;
  std::int64_t bead_id = 0;
  double weight = 1.0;
  int bead_type = 0;
  bool has_bead_type = false;
};

inline std::string trim_ws(const std::string& s) {
  std::size_t b = 0;
  while (b < s.size() && (s[b] == ' ' || s[b] == '\t' || s[b] == '\r' || s[b] == '\n')) ++b;
  std::size_t e = s.size();
  while (e > b && (s[e - 1] == ' ' || s[e - 1] == '\t' || s[e - 1] == '\r' || s[e - 1] == '\n')) --e;
  return s.substr(b, e - b);
}

inline std::vector<std::string> split_ws(const std::string& s) {
  std::vector<std::string> out;
  std::size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r' || s[i] == '\n')) ++i;
    if (i >= s.size()) break;
    std::size_t j = i;
    while (j < s.size() && !(s[j] == ' ' || s[j] == '\t' || s[j] == '\r' || s[j] == '\n')) ++j;
    out.emplace_back(s.substr(i, j - i));
    i = j;
  }
  return out;
}

inline std::vector<MappingFileEntry> read_mapping_file(const fs::path& path) {
  std::ifstream ifs(path);
  if (!ifs) throw std::runtime_error("mapping: failed to open mapping file: " + path.string());

  std::vector<MappingFileEntry> out;
  std::string line;
  while (std::getline(ifs, line)) {
    // strip comments
    const auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    line = trim_ws(line);
    if (line.empty()) continue;

    const auto toks = split_ws(line);
    if (toks.size() < 2) continue;

    MappingFileEntry e;
    try {
      e.atom_id = std::stoll(toks[0]);
      e.bead_id = std::stoll(toks[1]);
    } catch (...) {
      throw std::runtime_error("mapping: invalid integer in mapping file line: '" + line + "'");
    }
    if (toks.size() >= 3) {
      try {
        e.weight = std::stod(toks[2]);
      } catch (...) {
        throw std::runtime_error("mapping: invalid weight in mapping file line: '" + line + "'");
      }
    }
    if (toks.size() >= 4) {
      try {
        e.bead_type = std::stoi(toks[3]);
        e.has_bead_type = true;
      } catch (...) {
        throw std::runtime_error("mapping: invalid bead_type in mapping file line: '" + line + "'");
      }
    }

    if (!(e.weight >= 0.0)) {
      throw std::runtime_error("mapping: weight must be >= 0 in mapping file line: '" + line + "'");
    }

    out.push_back(e);
  }
  return out;
}

// BeadMapping defines a stable mapping from atoms -> beads, and provides a way to compute
// bead-level positions (unwrapped) for each trajectory frame.
//
// - Bead indices are 0..nbeads-1, in a deterministic order:
//     identity: ascending atom id (within the optional source_group)
//     by_mol:   ascending mol id (within the optional source_group)
//     file:     ascending bead_id as given by the mapping file
// - bead_id is an external stable identifier (e.g., atom id / mol id / user bead id).
struct BeadMapping {
  static constexpr std::size_t kInvalid = static_cast<std::size_t>(-1);

  MappingSpec spec;

  std::size_t natoms = 0;
  std::size_t nbeads = 0;

  // External stable bead ids.
  std::vector<std::int64_t> bead_id;
  std::vector<int> bead_type; // optional; size==nbeads or empty

  // For each bead b, list of member atom indices and per-atom weights.
  std::vector<std::vector<std::size_t>> bead_atoms;
  std::vector<std::vector<double>> bead_weights;

  // atom_to_bead[i] gives the bead index for atom i (or kInvalid if unmapped).
  std::vector<std::size_t> atom_to_bead;

  // Optional masses for COM-mass mapping.
  bool has_masses = false;
  std::vector<double> atom_mass;  // size natoms (1.0 if unknown)
  std::vector<double> bead_mass;  // size nbeads (sum of atom masses)

  bool has_atom_masses() const { return has_masses; }
  bool has_masses_for_com() const { return has_masses; }

  // Construct identity mapping (beads == atoms).
  static BeadMapping build_identity(const Frame& f,
                                   const MappingSpec& spec,
                                   std::span<const std::size_t> atoms_filter = {}) {
    BeadMapping m;
    m.spec = spec;
    m.natoms = f.natoms;

    // Determine which atoms are included.
    std::vector<std::size_t> atoms;
    if (!atoms_filter.empty()) {
      atoms.assign(atoms_filter.begin(), atoms_filter.end());
    } else {
      atoms.resize(f.natoms);
      for (std::size_t i = 0; i < f.natoms; ++i) atoms[i] = i;
    }

    // Sort by atom id for determinism.
    std::sort(atoms.begin(), atoms.end(), [&](std::size_t a, std::size_t b) {
      return f.id[a] < f.id[b];
    });

    m.nbeads = atoms.size();
    m.bead_id.resize(m.nbeads);
    m.bead_atoms.assign(m.nbeads, {});
    m.bead_weights.assign(m.nbeads, {});
    m.atom_to_bead.assign(f.natoms, kInvalid);

    for (std::size_t bi = 0; bi < atoms.size(); ++bi) {
      const std::size_t ai = atoms[bi];
      m.bead_id[bi] = f.id[ai];
      m.bead_atoms[bi] = {ai};
      m.bead_weights[bi] = {1.0};
      m.atom_to_bead[ai] = bi;
    }

    return m;
  }

  // Construct mapping by mol id: each mol -> one bead.
  static BeadMapping build_by_mol_ids(const Frame& f,
                                      std::span<const std::int64_t> mol_ids,
                                      const MappingSpec& spec,
                                      std::span<const std::size_t> atoms_filter = {}) {
    if (mol_ids.size() != f.natoms) {
      throw std::runtime_error("mapping: build_by_mol_ids: mol_ids size mismatch");
    }

    BeadMapping m;
    m.spec = spec;
    m.natoms = f.natoms;
    m.atom_to_bead.assign(f.natoms, kInvalid);

    // Collect atoms (possibly filtered).
    std::vector<std::size_t> atoms;
    if (!atoms_filter.empty()) {
      atoms.assign(atoms_filter.begin(), atoms_filter.end());
    } else {
      atoms.resize(f.natoms);
      for (std::size_t i = 0; i < f.natoms; ++i) atoms[i] = i;
    }

    // Group atoms by mol id.
    std::unordered_map<std::int64_t, std::vector<std::size_t>> mol_to_atoms;
    mol_to_atoms.reserve(atoms.size() / 4 + 8);
    for (std::size_t ai : atoms) {
      const std::int64_t mol = mol_ids[ai];
      mol_to_atoms[mol].push_back(ai);
    }

    // Deterministic bead order: ascending mol id.
    std::vector<std::int64_t> mol_keys;
    mol_keys.reserve(mol_to_atoms.size());
    for (const auto& kv : mol_to_atoms) mol_keys.push_back(kv.first);
    std::sort(mol_keys.begin(), mol_keys.end());

    m.nbeads = mol_keys.size();
    m.bead_id.resize(m.nbeads);
    m.bead_atoms.resize(m.nbeads);
    m.bead_weights.resize(m.nbeads);

    for (std::size_t bi = 0; bi < mol_keys.size(); ++bi) {
      const std::int64_t mol = mol_keys[bi];
      m.bead_id[bi] = mol;
      auto& vec = mol_to_atoms.at(mol);
      std::sort(vec.begin(), vec.end());
      m.bead_atoms[bi] = vec;
      m.bead_weights[bi].assign(vec.size(), 1.0);
      for (std::size_t ai : vec) m.atom_to_bead[ai] = bi;
    }

    return m;
  }

  static BeadMapping build_by_mol(const Frame& f,
                                  const MappingSpec& spec,
                                  std::span<const std::size_t> atoms_filter = {}) {
    if (!f.has_mol) {
      throw std::runtime_error("mapping: mode=by_mol requires 'mol' field (or topology-derived mol ids at the platform layer)");
    }
    return build_by_mol_ids(f, std::span<const std::int64_t>(f.mol.data(), f.mol.size()), spec, atoms_filter);
  }

  // Construct mapping from a mapping file.
  // The mapping file lists atom_id -> bead_id (+ optional weight and bead_type).
  static BeadMapping build_from_file(const Frame& f,
                                     const MappingSpec& spec,
                                     const fs::path& path,
                                     std::span<const std::size_t> atoms_filter = {}) {
    BeadMapping m;
    m.spec = spec;
    m.natoms = f.natoms;
    m.atom_to_bead.assign(f.natoms, kInvalid);

    std::vector<unsigned char> keep;
    if (!atoms_filter.empty()) {
      keep.assign(f.natoms, 0);
      for (std::size_t ai : atoms_filter) {
        if (ai < keep.size()) keep[ai] = 1;
      }
    }

    const auto entries = read_mapping_file(path);
    if (entries.empty()) {
      throw std::runtime_error("mapping: mapping file is empty: " + path.string());
    }

    // Build atom_id -> idx map.
    std::unordered_map<std::int64_t, std::size_t> id2idx;
    id2idx.reserve(f.natoms * 2);
    for (std::size_t i = 0; i < f.natoms; ++i) {
      id2idx.emplace(f.id[i], i);
    }

    // bead_id -> (atoms, weights)
    struct BeadTmp {
      std::vector<std::size_t> atoms;
      std::vector<double> weights;
      int bead_type = 0;
      bool has_type = false;
    };
    std::unordered_map<std::int64_t, BeadTmp> bead_map;
    bead_map.reserve(entries.size() / 2 + 8);

    // Validate uniqueness of atom_id assignments.
    std::unordered_map<std::size_t, std::int64_t> atom_assigned;
    atom_assigned.reserve(entries.size() * 2);

    for (const auto& e : entries) {
      auto it = id2idx.find(e.atom_id);
      if (it == id2idx.end()) {
        throw std::runtime_error("mapping: mapping file references unknown atom_id " + std::to_string(e.atom_id));
      }
      const std::size_t ai = it->second;
      if (!keep.empty() && !keep[ai]) {
        continue; // filtered out
      }

      auto it_a = atom_assigned.find(ai);
      if (it_a != atom_assigned.end() && it_a->second != e.bead_id) {
        throw std::runtime_error("mapping: atom_id " + std::to_string(e.atom_id) + " mapped to multiple bead_ids");
      }
      atom_assigned[ai] = e.bead_id;

      auto& bt = bead_map[e.bead_id];
      bt.atoms.push_back(ai);
      bt.weights.push_back(e.weight);
      if (e.has_bead_type) {
        if (bt.has_type && bt.bead_type != e.bead_type) {
          throw std::runtime_error("mapping: bead_id " + std::to_string(e.bead_id) + " has inconsistent bead_type in mapping file");
        }
        bt.bead_type = e.bead_type;
        bt.has_type = true;
      }
    }

    // Deterministic bead order: ascending bead_id.
    std::vector<std::int64_t> bead_ids;
    bead_ids.reserve(bead_map.size());
    for (const auto& kv : bead_map) bead_ids.push_back(kv.first);
    std::sort(bead_ids.begin(), bead_ids.end());

    m.nbeads = bead_ids.size();
    m.bead_id.resize(m.nbeads);
    m.bead_atoms.resize(m.nbeads);
    m.bead_weights.resize(m.nbeads);

    bool any_type = false;
    for (const auto& kv : bead_map) {
      if (kv.second.has_type) { any_type = true; break; }
    }
    if (any_type) m.bead_type.assign(m.nbeads, 0);

    for (std::size_t bi = 0; bi < bead_ids.size(); ++bi) {
      const std::int64_t bid = bead_ids[bi];
      m.bead_id[bi] = bid;
      auto& bt = bead_map.at(bid);
      // sort atoms for determinism, and permute weights accordingly
      std::vector<std::pair<std::size_t, double>> pairs;
      pairs.reserve(bt.atoms.size());
      for (std::size_t k = 0; k < bt.atoms.size(); ++k) {
        pairs.emplace_back(bt.atoms[k], bt.weights[k]);
      }
      std::sort(pairs.begin(), pairs.end(), [](auto a, auto b) { return a.first < b.first; });
      m.bead_atoms[bi].resize(pairs.size());
      m.bead_weights[bi].resize(pairs.size());
      for (std::size_t k = 0; k < pairs.size(); ++k) {
        m.bead_atoms[bi][k] = pairs[k].first;
        m.bead_weights[bi][k] = pairs[k].second;
        m.atom_to_bead[pairs[k].first] = bi;
      }
      if (any_type && bt.has_type) {
        m.bead_type[bi] = bt.bead_type;
      }
    }

    return m;
  }

  // Populate atom_mass / bead_mass (if masses are available from topology+type).
  void compute_masses_from(const Frame& first_frame, const Topology* topo) {
    atom_mass.assign(natoms, 1.0);
    bead_mass.assign(nbeads, 0.0);
    has_masses = false;

    if (!topo || !topo->has_masses()) {
      return;
    }
    if (!first_frame.has_type) {
      return;
    }

    // Masses are 1-indexed by type.
    has_masses = true;
    for (std::size_t i = 0; i < natoms; ++i) {
      const int t = first_frame.type[i];
      double m = 1.0;
      if (t >= 0 && static_cast<std::size_t>(t) < topo->mass_by_type.size()) {
        const double mm = topo->mass_by_type[static_cast<std::size_t>(t)];
        if (mm > 0.0) m = mm;
      }
      atom_mass[i] = m;
    }

    for (std::size_t b = 0; b < nbeads; ++b) {
      double sum = 0.0;
      for (std::size_t k = 0; k < bead_atoms[b].size(); ++k) {
        sum += atom_mass[bead_atoms[b][k]];
      }
      bead_mass[b] = sum;
    }
  }

  // Compute bead positions for a frame.
  BeadFrame compute_bead_frame(const Frame& f) const {
    if (f.natoms != natoms) {
      throw std::runtime_error("mapping: frame.natoms mismatch");
    }

    BeadFrame bf;
    bf.bead_id = bead_id;
    bf.xu.assign(nbeads, 0.0);
    bf.yu.assign(nbeads, 0.0);
    bf.zu.assign(nbeads, 0.0);
    bf.bead_mass = bead_mass; // may be empty or size nbeads
    bf.bead_type = bead_type;

    for (std::size_t b = 0; b < nbeads; ++b) {
      const auto& atoms = bead_atoms[b];
      const auto& w = bead_weights[b];
      if (atoms.empty()) continue;

      if (spec.position == MappingPosition::RepresentativeAtom) {
        // stable representative: smallest atom index
        std::size_t best = atoms[0];
        for (std::size_t k = 1; k < atoms.size(); ++k) {
          if (atoms[k] < best) best = atoms[k];
        }
        bf.xu[b] = f.xu[best];
        bf.yu[b] = f.yu[best];
        bf.zu[b] = f.zu[best];
        continue;
      }

      double wx = 0.0, wy = 0.0, wz = 0.0;
      double wsum = 0.0;

      for (std::size_t k = 0; k < atoms.size(); ++k) {
        const std::size_t ai = atoms[k];
        double ww = (k < w.size() ? w[k] : 1.0);
        if (spec.position == MappingPosition::ComMass) {
          if (!has_masses) {
            throw std::runtime_error("mapping: position=com_mass requested but masses are not available (need topology masses + dump type)");
          }
          ww *= atom_mass[ai];
        }
        if (ww == 0.0) continue;
        wx += ww * f.xu[ai];
        wy += ww * f.yu[ai];
        wz += ww * f.zu[ai];
        wsum += ww;
      }

      if (wsum > 0.0) {
        bf.xu[b] = wx / wsum;
        bf.yu[b] = wy / wsum;
        bf.zu[b] = wz / wsum;
      }
    }

    return bf;
  }
};

} // namespace pilots::alg::mapping
