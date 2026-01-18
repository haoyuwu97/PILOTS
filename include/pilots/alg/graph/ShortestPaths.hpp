#pragma once

#include <algorithm>
#include <cstddef>
#include <deque>
#include <stdexcept>
#include <vector>

#include "pilots/alg/graph/GraphView.hpp"

namespace pilots::alg::graph {

// Build an adjacency list from an EdgeList-backed graph.
inline Adjacency build_adjacency_from_edges(const GraphView& g) {
  const std::size_t n = g.node_count();
  Adjacency adj;
  adj.assign(n, {});
  g.for_each_edge([&](NodeId u, NodeId v) {
    if (u >= n || v >= n) throw std::runtime_error("build_adjacency_from_edges: out of range");
    adj[u].push_back(v);
    adj[v].push_back(u);
  });
  // Sort/unique for deterministic behavior.
  for (auto& nbrs : adj) {
    std::sort(nbrs.begin(), nbrs.end());
    nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
  }
  return adj;
}

// Unweighted BFS shortest path distances from a single source.
// Returns dist[i] = number of edges on shortest path, or -1 if unreachable.
inline std::vector<int> single_source_dist(const GraphView& g, NodeId source) {
  const std::size_t n = g.node_count();
  if (source >= n) throw std::runtime_error("single_source_dist: source out of range");

  const Adjacency* adj_ptr = g.has_adjacency() ? &g.adjacency() : nullptr;
  Adjacency tmp;
  if (!adj_ptr) {
    tmp = build_adjacency_from_edges(g);
    adj_ptr = &tmp;
  }

  std::vector<int> dist(n, -1);
  std::deque<NodeId> q;
  dist[source] = 0;
  q.push_back(source);

  while (!q.empty()) {
    NodeId u = q.front();
    q.pop_front();
    const int du = dist[u];
    for (NodeId v : (*adj_ptr)[u]) {
      if (dist[v] == -1) {
        dist[v] = du + 1;
        q.push_back(v);
      }
    }
  }
  return dist;
}

// Unweighted BFS distances from multiple sources.
inline std::vector<int> multi_source_dist(const GraphView& g, const std::vector<NodeId>& sources) {
  const std::size_t n = g.node_count();

  const Adjacency* adj_ptr = g.has_adjacency() ? &g.adjacency() : nullptr;
  Adjacency tmp;
  if (!adj_ptr) {
    tmp = build_adjacency_from_edges(g);
    adj_ptr = &tmp;
  }

  std::vector<int> dist(n, -1);
  std::deque<NodeId> q;
  for (NodeId s : sources) {
    if (s >= n) throw std::runtime_error("multi_source_dist: source out of range");
    if (dist[s] == 0) continue;
    dist[s] = 0;
    q.push_back(s);
  }

  while (!q.empty()) {
    NodeId u = q.front();
    q.pop_front();
    const int du = dist[u];
    for (NodeId v : (*adj_ptr)[u]) {
      if (dist[v] == -1) {
        dist[v] = du + 1;
        q.push_back(v);
      }
    }
  }

  return dist;
}

// Binned accumulation by graph distance.
struct GraphDistanceBinned {
  std::vector<std::size_t> count;
  std::vector<double> sum;
  std::vector<double> sumsq;

  std::size_t max_distance() const {
    return count.empty() ? 0 : (count.size() - 1);
  }
};

// Bin a scalar field by graph distance.
// For each node i with dist[i] >= 0, accumulates value[i] into bin dist[i].
inline GraphDistanceBinned graph_distance_histogram(const std::vector<int>& dist,
                                                    const std::vector<double>& value) {
  if (dist.size() != value.size()) {
    throw std::runtime_error("graph_distance_histogram: dist/value size mismatch");
  }

  int dmax = -1;
  for (int d : dist) if (d > dmax) dmax = d;
  const std::size_t nbins = (dmax >= 0) ? static_cast<std::size_t>(dmax + 1) : 0;

  GraphDistanceBinned out;
  out.count.assign(nbins, 0);
  out.sum.assign(nbins, 0.0);
  out.sumsq.assign(nbins, 0.0);

  for (std::size_t i = 0; i < dist.size(); ++i) {
    const int d = dist[i];
    if (d < 0) continue;
    const std::size_t bin = static_cast<std::size_t>(d);
    ++out.count[bin];
    out.sum[bin] += value[i];
    out.sumsq[bin] += value[i] * value[i];
  }

  return out;
}

} // namespace pilots::alg::graph
