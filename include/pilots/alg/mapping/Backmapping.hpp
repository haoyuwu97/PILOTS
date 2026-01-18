#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "pilots/alg/mapping/BeadMapping.hpp"

namespace pilots::alg::mapping {

// Backmapping utilities: scatter bead-level values back to per-atom arrays.
//
// This is intentionally "tool-level" in v1: it does not write trajectory files; it only
// constructs a per-atom array that measures can write as an output column.
struct BackMapper {
  // Broadcast an intensive bead value to all atoms belonging to that bead.
  //
  // Unmapped atoms are set to fill_value.
  static std::vector<double> scatter_intensive(const BeadMapping& m,
                                              std::span<const double> bead_value,
                                              double fill_value = 0.0) {
    std::vector<double> out(m.natoms, fill_value);
    for (std::size_t i = 0; i < m.natoms; ++i) {
      const std::size_t b = m.atom_to_bead[i];
      if (b == BeadMapping::kInvalid) continue;
      if (b >= bead_value.size()) continue;
      out[i] = bead_value[b];
    }
    return out;
  }

  // Scatter an extensive bead value to atoms proportionally to the mapping weights.
  //
  // Example: distributing a bead's total mass/charge to atoms.
  static std::vector<double> scatter_extensive_by_weight(const BeadMapping& m,
                                                         std::span<const double> bead_value,
                                                         double fill_value = 0.0) {
    std::vector<double> out(m.natoms, fill_value);
    for (std::size_t b = 0; b < m.nbeads; ++b) {
      if (b >= bead_value.size()) break;
      const auto& atoms = m.bead_atoms[b];
      const auto& w = m.bead_weights[b];
      double wsum = 0.0;
      for (double ww : w) wsum += ww;
      if (wsum <= 0.0) {
        // fall back to uniform
        const double v = atoms.empty() ? 0.0 : (bead_value[b] / static_cast<double>(atoms.size()));
        for (std::size_t ai : atoms) out[ai] = v;
        continue;
      }
      for (std::size_t k = 0; k < atoms.size(); ++k) {
        const double ww = (k < w.size() ? w[k] : 1.0);
        out[atoms[k]] = bead_value[b] * (ww / wsum);
      }
    }
    return out;
  }
};

} // namespace pilots::alg::mapping
