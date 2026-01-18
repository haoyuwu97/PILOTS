#pragma once

#include <cstddef>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "pilots/alg/graph/EdgeList.hpp"

namespace pilots::alg::graph {

// Adjacency list representation expected to be an undirected simple graph.
// (In practice, duplicates won't break most algorithms but will skew degrees.)
using Adjacency = std::vector<std::vector<NodeId>>;

// A lightweight, non-owning graph view.
//
// Supported backends:
//  - EdgeList: (n_nodes, edges)
//  - Adjacency: vector<vector<NodeId>>
//
// This object intentionally does NOT couple to Runner/Measure. It's a pure data view.
class GraphView {
public:
  GraphView() = default;

  explicit GraphView(const EdgeList& el)
      : el_(&el), adj_(nullptr) {
    if (el.n_nodes == 0 && !el.edges.empty()) {
      throw std::runtime_error("GraphView: EdgeList has edges but n_nodes=0");
    }
  }

  explicit GraphView(const Adjacency& adj)
      : el_(nullptr), adj_(&adj) {}

  bool has_edge_list() const { return el_ != nullptr; }
  bool has_adjacency() const { return adj_ != nullptr; }

  std::size_t node_count() const {
    if (el_) return el_->n_nodes;
    if (adj_) return adj_->size();
    return 0;
  }

  // Undirected edge count.
  // For Adjacency we assume symmetric adjacency lists and return sum(deg)/2.
  std::size_t edge_count() const {
    if (el_) return el_->edges.size();
    if (adj_) {
      std::size_t sum = 0;
      for (const auto& nbrs : *adj_) sum += nbrs.size();
      return sum / 2;
    }
    return 0;
  }

  const Adjacency& adjacency() const {
    if (!adj_) throw std::runtime_error("GraphView: adjacency() requested but backend is EdgeList");
    return *adj_;
  }

  // Iterate each undirected edge once.
  // For Adjacency, this assumes symmetric adjacency lists.
  template <class F>
  void for_each_edge(F&& f) const {
    if (el_) {
      for (const auto& e : el_->edges) {
        if (e.first == e.second) continue;
        f(e.first, e.second);
      }
      return;
    }
    if (adj_) {
      const std::size_t n = adj_->size();
      for (std::size_t u = 0; u < n; ++u) {
        for (NodeId v : (*adj_)[u]) {
          if (u < v) {
            f(u, v);
          }
        }
      }
      return;
    }
  }

  // Degree vector (cached). For Adjacency this is just nbr-list sizes.
  const std::vector<std::size_t>& degree() const {
    if (degree_cache_.has_value()) return *degree_cache_;
    degree_cache_.emplace();
    auto& deg = *degree_cache_;
    const std::size_t n = node_count();
    deg.assign(n, 0);

    if (adj_) {
      for (std::size_t i = 0; i < n; ++i) {
        deg[i] = (*adj_)[i].size();
      }
      return deg;
    }

    if (el_) {
      for (const auto& e : el_->edges) {
        if (e.first >= n || e.second >= n) {
          throw std::runtime_error("GraphView: edge endpoint out of range");
        }
        if (e.first == e.second) continue;
        ++deg[e.first];
        ++deg[e.second];
      }
      return deg;
    }

    return deg;
  }

  void clear_caches() const {
    degree_cache_.reset();
  }

private:
  const EdgeList* el_ = nullptr;
  const Adjacency* adj_ = nullptr;
  mutable std::optional<std::vector<std::size_t>> degree_cache_;
};

} // namespace pilots::alg::graph
