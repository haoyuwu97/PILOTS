#pragma once

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pilots/alg/graph/GraphView.hpp"

namespace pilots::alg::graph {

// Disjoint-set / union-find for undirected connectivity.
class UnionFind {
public:
  explicit UnionFind(std::size_t n = 0) { reset(n); }

  void reset(std::size_t n) {
    parent_.resize(n);
    rank_.assign(n, 0);
    std::iota(parent_.begin(), parent_.end(), 0);
  }

  std::size_t find(std::size_t x) {
    if (x >= parent_.size()) throw std::runtime_error("UnionFind: find() out of range");
    while (parent_[x] != x) {
      parent_[x] = parent_[parent_[x]];
      x = parent_[x];
    }
    return x;
  }

  void unite(std::size_t a, std::size_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return;
    if (rank_[a] < rank_[b]) std::swap(a, b);
    parent_[b] = a;
    if (rank_[a] == rank_[b]) ++rank_[a];
  }

private:
  std::vector<std::size_t> parent_;
  std::vector<unsigned char> rank_;
};

struct ComponentsResult {
  std::vector<std::size_t> component_id;    // size = n_nodes
  std::vector<std::size_t> component_sizes; // indexed by component_id
  std::size_t largest_component_size = 0;

  std::size_t n_components() const { return component_sizes.size(); }

  // Histogram: size -> count
  std::unordered_map<std::size_t, std::size_t> size_histogram() const {
    std::unordered_map<std::size_t, std::size_t> h;
    h.reserve(component_sizes.size());
    for (auto s : component_sizes) {
      ++h[s];
    }
    return h;
  }
};

inline ComponentsResult compute_components(const GraphView& g) {
  const std::size_t n = g.node_count();
  UnionFind uf(n);

  g.for_each_edge([&](NodeId u, NodeId v) {
    if (u >= n || v >= n) throw std::runtime_error("compute_components: edge endpoint out of range");
    uf.unite(u, v);
  });

  // Compress and remap root -> component id [0..C-1]
  std::unordered_map<std::size_t, std::size_t> root_to_cid;
  root_to_cid.reserve(n);

  ComponentsResult r;
  r.component_id.assign(n, static_cast<std::size_t>(-1));

  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t root = uf.find(i);
    auto it = root_to_cid.find(root);
    if (it == root_to_cid.end()) {
      const std::size_t cid = root_to_cid.size();
      root_to_cid.emplace(root, cid);
      r.component_sizes.push_back(0);
      r.component_id[i] = cid;
      ++r.component_sizes[cid];
    } else {
      const std::size_t cid = it->second;
      r.component_id[i] = cid;
      ++r.component_sizes[cid];
    }
  }

  r.largest_component_size = 0;
  for (auto s : r.component_sizes) r.largest_component_size = std::max(r.largest_component_size, s);

  return r;
}

struct ClusterMoments {
  // Raw moments Mk = sum over clusters of (size^k)
  // m0 = number of clusters.
  double m0 = 0.0;
  double m1 = 0.0;
  double m2 = 0.0;
  double m3 = 0.0;

  double weight_average_size() const {
    return (m1 > 0.0) ? (m2 / m1) : 0.0;
  }
};

inline ClusterMoments cluster_moments_0_3(const ComponentsResult& comps) {
  ClusterMoments cm;
  cm.m0 = static_cast<double>(comps.component_sizes.size());
  for (std::size_t s : comps.component_sizes) {
    const double sd = static_cast<double>(s);
    cm.m1 += sd;
    cm.m2 += sd * sd;
    cm.m3 += sd * sd * sd;
  }
  return cm;
}

// Non-PBC spanning check.
//
// v1 intent: provide a simple, clearly-audited spanning criterion that does NOT
// claim full periodic percolation correctness.
struct NonPBCSpanningResult {
  // per-component
  std::vector<bool> spans_x;
  std::vector<bool> spans_y;
  std::vector<bool> spans_z;

  bool any_spans_x = false;
  bool any_spans_y = false;
  bool any_spans_z = false;

  bool any_spans() const { return any_spans_x || any_spans_y || any_spans_z; }
};

} // namespace pilots::alg::graph
