#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pilots {

// An analysis graph built from a Topology's bond list, with user-configurable filtering.
//
// This is intentionally a light-weight container:
// - adjacency is built once (mol-first workflows typically need fast neighbor queries)
// - connected components are computed lazily and cached
// - no "backbone semantics" are imposed; the user controls what bonds/nodes are included via config
struct BondGraph {
  std::string mode = "all_bonds";
  std::size_t natoms = 0;

  // Auditable definition
  std::vector<int> bond_types_filter; // empty = all
  std::string atom_group_filter;      // empty = no atom restriction

  // Graph storage
  std::vector<std::vector<std::size_t>> adjacency; // size natoms

  // Cached connected components
  mutable bool comps_valid = false;
  mutable std::vector<std::size_t> comp_id;   // per atom, component index
  mutable std::vector<std::size_t> comp_size; // per component

  std::size_t degree(std::size_t i) const { return (i < adjacency.size()) ? adjacency[i].size() : 0; }

  void ensure_components() const {
    if (comps_valid) return;
    const std::size_t N = natoms;
    comp_id.assign(N, static_cast<std::size_t>(-1));
    comp_size.clear();
    comp_size.reserve(N);

    std::deque<std::size_t> q;
    std::size_t cid = 0;
    for (std::size_t s = 0; s < N; ++s) {
      if (comp_id[s] != static_cast<std::size_t>(-1)) continue;
      comp_id[s] = cid;
      q.clear();
      q.push_back(s);
      std::size_t sz = 0;
      while (!q.empty()) {
        const std::size_t u = q.front();
        q.pop_front();
        ++sz;
        if (u >= adjacency.size()) continue;
        for (const std::size_t v : adjacency[u]) {
          if (v >= N) continue;
          if (comp_id[v] == static_cast<std::size_t>(-1)) {
            comp_id[v] = cid;
            q.push_back(v);
          }
        }
      }
      comp_size.push_back(sz);
      ++cid;
    }

    comps_valid = true;
  }
};

} // namespace pilots
