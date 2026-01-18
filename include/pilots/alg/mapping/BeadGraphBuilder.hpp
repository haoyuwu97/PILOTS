#pragma once

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "pilots/alg/graph/EdgeList.hpp"
#include "pilots/alg/mapping/BeadMapping.hpp"
#include "pilots/topology/Topology.hpp"

namespace pilots::alg::mapping {

// Build a bead-level undirected simple graph by compressing atom-level bonds.
//
// - Each bond (i,j) maps to (bi,bj) using atom_to_bead.
// - Bonds with unmapped endpoints or self-edges (bi==bj) are dropped.
// - Multiple bonds between the same bead pair are uniquified.
inline pilots::alg::graph::EdgeList build_bead_graph_from_bonds(const pilots::Topology& topo,
                                                                const BeadMapping& map) {
  pilots::alg::graph::EdgeList el;
  el.n_nodes = map.nbeads;
  if (!topo.has_bonds() || map.nbeads == 0) return el;

  el.edges.reserve(topo.bonds.size());

  for (const auto& b : topo.bonds) {
    if (b.i >= map.atom_to_bead.size() || b.j >= map.atom_to_bead.size()) continue;
    const std::size_t bi = map.atom_to_bead[b.i];
    const std::size_t bj = map.atom_to_bead[b.j];
    if (bi == BeadMapping::kInvalid || bj == BeadMapping::kInvalid) continue;
    if (bi == bj) continue;
    const std::size_t u = (bi < bj) ? bi : bj;
    const std::size_t v = (bi < bj) ? bj : bi;
    el.edges.emplace_back(u, v);
  }

  std::sort(el.edges.begin(), el.edges.end());
  el.edges.erase(std::unique(el.edges.begin(), el.edges.end()), el.edges.end());
  return el;
}

} // namespace pilots::alg::mapping
