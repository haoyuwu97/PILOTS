#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace pilots::alg::graph {

using NodeId = std::size_t;
using Edge = std::pair<NodeId, NodeId>; // undirected edge (u,v)

struct EdgeList {
  std::size_t n_nodes = 0;
  std::vector<Edge> edges;

  std::size_t node_count() const { return n_nodes; }
  std::size_t edge_count() const { return edges.size(); }
};

} // namespace pilots::alg::graph
