#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "pilots/alg/graph/Components.hpp"

namespace pilots::alg::graph {

struct GraphSummary {
  std::size_t node_count = 0;
  std::size_t edge_count = 0;
  std::size_t n_components = 0;

  std::size_t n_ends = 0;          // degree==1
  std::size_t n_branch_points = 0; // degree>=3

  // cycle_rank = E - V + C for undirected graphs
  std::int64_t cycle_rank = 0;

  // degree -> count
  std::unordered_map<std::size_t, std::size_t> degree_hist;
};

inline std::unordered_map<std::size_t, std::size_t> degree_distribution(const GraphView& g) {
  std::unordered_map<std::size_t, std::size_t> h;
  const auto& deg = g.degree();
  h.reserve(deg.size());
  for (std::size_t d : deg) {
    ++h[d];
  }
  return h;
}

inline std::vector<NodeId> ends(const GraphView& g) {
  const auto& deg = g.degree();
  std::vector<NodeId> out;
  out.reserve(deg.size());
  for (std::size_t i = 0; i < deg.size(); ++i) {
    if (deg[i] == 1) out.push_back(i);
  }
  return out;
}

inline std::vector<NodeId> branch_points(const GraphView& g, std::size_t min_degree = 3) {
  const auto& deg = g.degree();
  std::vector<NodeId> out;
  out.reserve(deg.size());
  for (std::size_t i = 0; i < deg.size(); ++i) {
    if (deg[i] >= min_degree) out.push_back(i);
  }
  return out;
}

inline GraphSummary graph_summary(const GraphView& g, const ComponentsResult& comps) {
  GraphSummary s;
  s.node_count = g.node_count();
  s.edge_count = g.edge_count();
  s.n_components = comps.component_sizes.size();

  const auto& deg = g.degree();
  for (std::size_t d : deg) {
    if (d == 1) ++s.n_ends;
    if (d >= 3) ++s.n_branch_points;
  }

  // cycle_rank = E - V + C
  s.cycle_rank = static_cast<std::int64_t>(s.edge_count)
               - static_cast<std::int64_t>(s.node_count)
               + static_cast<std::int64_t>(s.n_components);

  s.degree_hist = degree_distribution(g);
  return s;
}

// Convenience: compute components internally.
inline GraphSummary graph_summary(const GraphView& g) {
  return graph_summary(g, compute_components(g));
}

} // namespace pilots::alg::graph
