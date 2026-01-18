#pragma once

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include "pilots/alg/geom/GeometryView.hpp"
#include "pilots/alg/graph/EdgeList.hpp"

namespace pilots::alg::neighbor {

struct ContactGraph {
  pilots::alg::graph::EdgeList edge_list;
  // map node id -> underlying particle index (atom idx in Frame)
  std::vector<std::size_t> node_to_particle;
};

// Naive O(N^2) contact graph builder (v1: correctness first).
//
// Builds a graph over the provided particle indices. Nodes are 0..(m-1)
// corresponding to entries of `particles`.
inline ContactGraph build_contact_graph_naive(
    const pilots::alg::geom::GeometryView& geo,
    std::span<const std::size_t> particles,
    double cutoff,
    bool include_self = false) {

  if (cutoff < 0.0) throw std::runtime_error("build_contact_graph_naive: cutoff must be nonnegative");
  const double rc2 = cutoff * cutoff;

  ContactGraph cg;
  cg.node_to_particle.assign(particles.begin(), particles.end());
  cg.edge_list.n_nodes = cg.node_to_particle.size();
  cg.edge_list.edges.clear();

  const std::size_t m = cg.node_to_particle.size();
  // Upper bound of edges is m*(m-1)/2; reserve a small fraction to avoid huge alloc.
  if (m > 0) {
    cg.edge_list.edges.reserve(std::min<std::size_t>(m * 4, (m * (m - 1)) / 2));
  }

  for (std::size_t a = 0; a < m; ++a) {
    const std::size_t i = cg.node_to_particle[a];
    if (i >= geo.size()) throw std::runtime_error("build_contact_graph_naive: particle index out of range");

    const std::size_t b0 = include_self ? a : (a + 1);
    for (std::size_t b = b0; b < m; ++b) {
      const std::size_t j = cg.node_to_particle[b];
      if (j >= geo.size()) throw std::runtime_error("build_contact_graph_naive: particle index out of range");

      if (!include_self && i == j) continue;
      const double d2 = geo.min_image_distance_sq(i, j);
      if (d2 <= rc2) {
        cg.edge_list.edges.emplace_back(a, b);
      }
    }
  }

  return cg;
}

} // namespace pilots::alg::neighbor
