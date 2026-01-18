#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pilots/core/Frame.hpp"

namespace pilots {

// MoleculeIndex: mol_id -> list of atom indices (canonical order).
//
// Design notes:
// - Built from the trajectory (dump) mol field, if present.
// - Kept in SystemContext for mol-first chain/molecule aggregation in measures.
struct MoleculeIndex {
  std::vector<std::int64_t> mol_ids;                   // unique mol ids (stable order)
  std::vector<std::vector<std::size_t>> mol_atoms;     // parallel to mol_ids
  std::unordered_map<std::int64_t, std::size_t> mol_id_to_slot;

  std::size_t natoms = 0;
  bool built = false;

  void clear() {
    mol_ids.clear();
    mol_atoms.clear();
    mol_id_to_slot.clear();
    natoms = 0;
    built = false;
  }

  void build_from_frame(const Frame& frame) {
    clear();
    if (!frame.has_mol) {
      throw std::runtime_error("MoleculeIndex: frame has no 'mol' field");
    }
    if (frame.mol.size() != frame.natoms) {
      throw std::runtime_error("MoleculeIndex: invalid 'mol' array length");
    }
    natoms = frame.natoms;

    // Build in first-seen order for reproducibility.
    mol_id_to_slot.reserve(natoms * 2);
    mol_ids.reserve(natoms);
    mol_atoms.reserve(natoms);

    for (std::size_t i = 0; i < natoms; ++i) {
      const std::int64_t mid = frame.mol[i];
      auto it = mol_id_to_slot.find(mid);
      if (it == mol_id_to_slot.end()) {
        const std::size_t slot = mol_ids.size();
        mol_id_to_slot.emplace(mid, slot);
        mol_ids.push_back(mid);
        mol_atoms.emplace_back();
        mol_atoms.back().push_back(i);
      } else {
        mol_atoms[it->second].push_back(i);
      }
    }
    built = true;
  }

  // Build from a per-atom mol id vector (canonical atom index order).
  // This is useful when the mol field is not present in the dump but derived elsewhere (e.g., from bond connectivity).
  void build_from_mol_vector(const std::vector<std::int64_t>& mol, std::size_t natoms_in) {
    clear();
    if (mol.size() != natoms_in) {
      throw std::runtime_error("MoleculeIndex: invalid mol vector length");
    }
    natoms = natoms_in;

    mol_id_to_slot.reserve(natoms * 2);
    mol_ids.reserve(natoms);
    mol_atoms.reserve(natoms);

    for (std::size_t i = 0; i < natoms; ++i) {
      const std::int64_t mid = mol[i];
      auto it = mol_id_to_slot.find(mid);
      if (it == mol_id_to_slot.end()) {
        const std::size_t slot = mol_ids.size();
        mol_id_to_slot.emplace(mid, slot);
        mol_ids.push_back(mid);
        mol_atoms.emplace_back();
        mol_atoms.back().push_back(i);
      } else {
        mol_atoms[it->second].push_back(i);
      }
    }
    built = true;
  }

  std::size_t nmol() const { return mol_ids.size(); }
};

} // namespace pilots
