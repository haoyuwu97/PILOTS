#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pilots/core/Frame.hpp"
#include "pilots/alg/graph/Components.hpp"

namespace pilots {

// Topology container for polymer-/network-type analyses.
//
// P1-D goals (this package):
// - Container can hold the common LAMMPS topology sections:
//   Masses, Bonds, Angles, Dihedrals, Impropers (all optional).
// - Reading is on-demand via a section mask (Runner aggregates needs).
// - Bonds are stored with bond_type and mapped to canonical indices (by atom 'id').
// - A default all-bonds adjacency graph is built when bonds are loaded.
// - Connected components are computed on-demand and cached.
struct Topology {
  // --- Loaded section flags (what was actually populated) ---
  bool loaded_masses = false;
  bool loaded_bonds = false;
  bool loaded_angles = false;
  bool loaded_dihedrals = false;
  bool loaded_impropers = false;

  // --- Optional per-type properties ---
  // LAMMPS Masses section: 1-indexed type -> mass (index 0 unused).
  std::vector<double> mass_by_type;

  // --- Raw (id-based) topology lists (source-of-truth, before canonicalization) ---
  struct BondId {
    int type = 0;
    std::int64_t ai = 0;
    std::int64_t aj = 0;
  };
  struct AngleId {
    int type = 0;
    std::int64_t ai = 0;
    std::int64_t aj = 0;
    std::int64_t ak = 0;
  };
  struct DihedralId {
    int type = 0;
    std::int64_t ai = 0;
    std::int64_t aj = 0;
    std::int64_t ak = 0;
    std::int64_t al = 0;
  };
  using ImproperId = DihedralId;

  std::vector<BondId> bond_ids;
  std::vector<AngleId> angle_ids;
  std::vector<DihedralId> dihedral_ids;
  std::vector<ImproperId> improper_ids;

  // --- Canonicalized (index-based) topology lists (derived from first frame) ---
  struct Bond {
    int type = 0;
    std::size_t i = 0;
    std::size_t j = 0;
  };
  struct Angle {
    int type = 0;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
  };
  struct Dihedral {
    int type = 0;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
    std::size_t l = 0;
  };
  using Improper = Dihedral;

  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Dihedral> dihedrals;
  std::vector<Improper> impropers;

  // Default all-bonds adjacency graph (built iff loaded_bonds).
  std::vector<std::vector<std::size_t>> adjacency; // adjacency[i] -> neighbor atom indices

  std::size_t natoms = 0;
  bool finalized = false;

  // --- Cached connected components for the default all-bonds graph ---
  // component_id[i] gives the component index for atom i.
  // component_sizes[c] gives the size of component c.
  mutable bool components_valid = false;
  mutable std::vector<std::size_t> component_id;
  mutable std::vector<std::size_t> component_sizes;

  void clear() {
    loaded_masses = loaded_bonds = loaded_angles = loaded_dihedrals = loaded_impropers = false;
    mass_by_type.clear();
    bond_ids.clear();
    angle_ids.clear();
    dihedral_ids.clear();
    improper_ids.clear();
    bonds.clear();
    angles.clear();
    dihedrals.clear();
    impropers.clear();
    adjacency.clear();
    natoms = 0;
    finalized = false;
    components_valid = false;
    component_id.clear();
    component_sizes.clear();
  }

  bool has_masses() const { return loaded_masses; }
  bool has_bonds() const { return loaded_bonds; }
  bool has_angles() const { return loaded_angles; }
  bool has_dihedrals() const { return loaded_dihedrals; }
  bool has_impropers() const { return loaded_impropers; }

  std::vector<std::string> loaded_sections() const {
    std::vector<std::string> out;
    if (loaded_masses) out.emplace_back("masses");
    if (loaded_bonds) out.emplace_back("bonds");
    if (loaded_angles) out.emplace_back("angles");
    if (loaded_dihedrals) out.emplace_back("dihedrals");
    if (loaded_impropers) out.emplace_back("impropers");
    return out;
  }

  // Derive a per-atom "mol id" vector from bond connectivity (connected components of TopologyGraph).
  //
  // Semantics:
  // - Requires bonds loaded and finalized mapping to current Frame.
  // - Component ids are mapped to 1..n_components (stable for a given graph).
  std::vector<std::int64_t> derived_mol_ids_from_bonds(const Frame& frame) const {
    if (!finalized) {
      throw std::runtime_error("Topology: cannot derive mol ids before finalize_from_frame");
    }
    if (frame.natoms != natoms) {
      throw std::runtime_error("Topology: frame.natoms mismatch in derived_mol_ids_from_bonds");
    }
    if (!loaded_bonds) {
      throw std::runtime_error("Topology: bonds not loaded; cannot derive mol ids from bonds");
    }
    ensure_components();
    std::vector<std::int64_t> mol(natoms, 0);
    for (std::size_t i = 0; i < natoms; ++i) {
      mol[i] = static_cast<std::int64_t>(component_id[i] + 1);
    }
    return mol;
  }

  // Build canonical index topology lists + adjacency from the first frame's id list.
  // Atom identity is determined by the LAMMPS 'id' field and is mapped to the reader's
  // canonical ordering (first frame file order).
  void finalize_from_frame(const Frame& frame) {
    natoms = frame.natoms;
    if (natoms == 0) throw std::runtime_error("Topology: frame.natoms is 0");
    if (frame.id.size() != natoms) throw std::runtime_error("Topology: frame.id missing/invalid");

    std::unordered_map<std::int64_t, std::size_t> id2idx;
    id2idx.reserve(natoms * 2);
    for (std::size_t i = 0; i < natoms; ++i) {
      auto [it, ok] = id2idx.emplace(frame.id[i], i);
      if (!ok) {
        throw std::runtime_error("Topology: duplicate atom id in frame: " + std::to_string(frame.id[i]));
      }
    }

    // Bonds + adjacency
    bonds.clear();
    adjacency.clear();
    if (loaded_bonds) {
      bonds.reserve(bond_ids.size());
      adjacency.assign(natoms, {});

      for (const auto& b : bond_ids) {
        auto it_i = id2idx.find(b.ai);
        auto it_j = id2idx.find(b.aj);
        if (it_i == id2idx.end() || it_j == id2idx.end()) {
          throw std::runtime_error("Topology: bond references unknown atom id(s): (" + std::to_string(b.ai) + ", " + std::to_string(b.aj) + ")");
        }
        const std::size_t i = it_i->second;
        const std::size_t j = it_j->second;
        if (i == j) continue;
        bonds.push_back(Bond{b.type, i, j});
        adjacency[i].push_back(j);
        adjacency[j].push_back(i);
      }

      // Ensure determinism + correct degree/end-point semantics even if the input contains
      // duplicate bond records. We do not attempt to uniquify the bonds list itself (which
      // might carry distinct bond types in unusual cases), but the adjacency should be a
      // simple graph for component/degree-based analyses.
      for (auto& nbrs : adjacency) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
      }
    }

    // Angles
    angles.clear();
    if (loaded_angles) {
      angles.reserve(angle_ids.size());
      for (const auto& a : angle_ids) {
        auto it_i = id2idx.find(a.ai);
        auto it_j = id2idx.find(a.aj);
        auto it_k = id2idx.find(a.ak);
        if (it_i == id2idx.end() || it_j == id2idx.end() || it_k == id2idx.end()) {
          throw std::runtime_error("Topology: angle references unknown atom id(s)");
        }
        angles.push_back(Angle{a.type, it_i->second, it_j->second, it_k->second});
      }
    }

    // Dihedrals
    dihedrals.clear();
    if (loaded_dihedrals) {
      dihedrals.reserve(dihedral_ids.size());
      for (const auto& d : dihedral_ids) {
        auto it_i = id2idx.find(d.ai);
        auto it_j = id2idx.find(d.aj);
        auto it_k = id2idx.find(d.ak);
        auto it_l = id2idx.find(d.al);
        if (it_i == id2idx.end() || it_j == id2idx.end() || it_k == id2idx.end() || it_l == id2idx.end()) {
          throw std::runtime_error("Topology: dihedral references unknown atom id(s)");
        }
        dihedrals.push_back(Dihedral{d.type, it_i->second, it_j->second, it_k->second, it_l->second});
      }
    }

    // Impropers
    impropers.clear();
    if (loaded_impropers) {
      impropers.reserve(improper_ids.size());
      for (const auto& d : improper_ids) {
        auto it_i = id2idx.find(d.ai);
        auto it_j = id2idx.find(d.aj);
        auto it_k = id2idx.find(d.ak);
        auto it_l = id2idx.find(d.al);
        if (it_i == id2idx.end() || it_j == id2idx.end() || it_k == id2idx.end() || it_l == id2idx.end()) {
          throw std::runtime_error("Topology: improper references unknown atom id(s)");
        }
        impropers.push_back(Improper{d.type, it_i->second, it_j->second, it_k->second, it_l->second});
      }
    }

    components_valid = false;
    component_id.clear();
    component_sizes.clear();

    finalized = true;
  }

  void ensure_components() const {
    if (components_valid) return;
    if (!finalized) throw std::runtime_error("Topology: ensure_components() requires finalized topology");

    // Reuse the platform algorithm layer for connectivity.
    // This ensures TopoGroups and future cluster/percolation measures share identical
    // connected-component semantics.
    pilots::alg::graph::GraphView gv(adjacency);
    const auto comps = pilots::alg::graph::compute_components(gv);

    component_id = comps.component_id;
    component_sizes = comps.component_sizes;
    components_valid = true;
  }
};

} // namespace pilots
