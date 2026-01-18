#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pilots/alg/graph/Components.hpp"
#include "pilots/alg/graph/GraphView.hpp"

namespace pilots::alg::index {

// ChainIndex: maps chain_id -> ordered node list (typically bead indices).
//
// This is a pure algorithm-level construct. It does NOT prescribe a physical meaning
// of "chain"; it only provides an auditable, reusable ordering.
struct ChainIndex {
  // External stable chain ids (e.g., mol id, user chain id). Same length as chains.
  std::vector<std::int64_t> chain_id;

  // chains[c] gives an ordered list of node indices for chain c.
  std::vector<std::vector<std::size_t>> chains;

  std::size_t n_chains() const { return chains.size(); }

  std::size_t total_nodes() const {
    std::size_t s = 0;
    for (const auto& ch : chains) s += ch.size();
    return s;
  }
};

inline ChainIndex build_chain_index_from_chain_pos(
    std::span<const std::int64_t> chain_id_per_node,
    std::span<const std::int64_t> chain_pos_per_node,
    bool strict = true) {

  if (chain_id_per_node.size() != chain_pos_per_node.size()) {
    throw std::runtime_error("ChainIndex: chain_id and chain_pos size mismatch");
  }

  // Group nodes by chain_id.
  std::unordered_map<std::int64_t, std::vector<std::pair<std::int64_t, std::size_t>>> tmp;
  tmp.reserve(chain_id_per_node.size() / 4 + 8);

  for (std::size_t i = 0; i < chain_id_per_node.size(); ++i) {
    const std::int64_t cid = chain_id_per_node[i];
    const std::int64_t pos = chain_pos_per_node[i];
    tmp[cid].push_back({pos, i});
  }

  // Deterministic chain order: ascending chain_id.
  std::vector<std::int64_t> chain_ids;
  chain_ids.reserve(tmp.size());
  for (const auto& kv : tmp) chain_ids.push_back(kv.first);
  std::sort(chain_ids.begin(), chain_ids.end());

  ChainIndex out;
  out.chain_id = chain_ids;
  out.chains.resize(chain_ids.size());

  for (std::size_t c = 0; c < chain_ids.size(); ++c) {
    auto& v = tmp[chain_ids[c]];
    std::sort(v.begin(), v.end(), [](auto a, auto b) {
      if (a.first != b.first) return a.first < b.first;
      return a.second < b.second;
    });

    // Validate uniqueness of positions if strict.
    if (strict) {
      for (std::size_t k = 1; k < v.size(); ++k) {
        if (v[k].first == v[k - 1].first) {
          throw std::runtime_error("ChainIndex: duplicate chain_pos=" + std::to_string(v[k].first) +
                                   " in chain_id=" + std::to_string(chain_ids[c]));
        }
      }
    }

    out.chains[c].reserve(v.size());
    for (auto [pos, idx] : v) {
      (void)pos;
      out.chains[c].push_back(idx);
    }
  }

  return out;
}

// Convenience overload: chain_pos may come from a double-valued extra field.
inline ChainIndex build_chain_index_from_chain_pos(
    std::span<const std::int64_t> chain_id_per_node,
    std::span<const double> chain_pos_per_node,
    bool strict = true) {

  std::vector<std::int64_t> pos(chain_pos_per_node.size(), 0);
  for (std::size_t i = 0; i < chain_pos_per_node.size(); ++i) {
    const double x = chain_pos_per_node[i];
    const std::int64_t p = static_cast<std::int64_t>(x);
    if (strict && (static_cast<double>(p) != x)) {
      throw std::runtime_error("ChainIndex: chain_pos contains non-integer value");
    }
    pos[i] = p;
  }
  return build_chain_index_from_chain_pos(chain_id_per_node,
                                          std::span<const std::int64_t>(pos.data(), pos.size()),
                                          strict);
}

// Build ChainIndex by extracting ordered paths from linear connected components.
//
// Rule (v1):
// - If any component is non-linear (ring/branched), strict=true throws.
// - If strict=false, non-linear components are skipped.
inline ChainIndex build_chain_index_from_linear_components(
    const pilots::alg::graph::GraphView& g,
    bool strict = true) {

  const auto comps = pilots::alg::graph::compute_components(g);
  const auto& deg = g.degree();
  const std::size_t n = g.node_count();

  // Count edges per component once (undirected edges).
  std::vector<std::size_t> edges_in_comp(comps.component_sizes.size(), 0);
  g.for_each_edge([&](pilots::alg::graph::NodeId u, pilots::alg::graph::NodeId v) {
    const std::size_t cu = comps.component_id[u];
    const std::size_t cv = comps.component_id[v];
    if (cu == cv) ++edges_in_comp[cu];
  });

  // Collect nodes per component.
  std::vector<std::vector<std::size_t>> nodes_in_comp(comps.component_sizes.size());
  for (std::size_t i = 0; i < n; ++i) {
    nodes_in_comp[comps.component_id[i]].push_back(i);
  }

  ChainIndex out;
  out.chain_id.reserve(nodes_in_comp.size());
  out.chains.reserve(nodes_in_comp.size());

  for (std::size_t cid = 0; cid < nodes_in_comp.size(); ++cid) {
    const auto& nodes = nodes_in_comp[cid];
    if (nodes.empty()) continue;

    // Determine if component is linear.
    std::size_t n_end = 0;
    std::size_t n_branch = 0;
    for (std::size_t v : nodes) {
      const std::size_t d = deg[v];
      if (d == 1) ++n_end;
      if (d >= 3) ++n_branch;
    }

    // Compute cycle_rank for this component: E - V + 1.
    const std::size_t E = edges_in_comp[cid];
    const std::int64_t cycle_rank = static_cast<std::int64_t>(E) - static_cast<std::int64_t>(nodes.size()) + 1;

    const bool is_linear = (cycle_rank == 0) && (n_branch == 0) && ((nodes.size() == 1) || (n_end == 2));
    if (!is_linear) {
      if (strict) {
        throw std::runtime_error("ChainIndex: non-linear component encountered (cid=" + std::to_string(cid) + ")");
      }
      continue;
    }

    // Find an end node (or smallest node if size==1).
    std::size_t start = nodes[0];
    if (nodes.size() > 1) {
      bool found = false;
      for (std::size_t v : nodes) {
        if (deg[v] == 1) {
          start = v;
          found = true;
          break;
        }
      }
      if (!found) {
        // size==2 case still has ends; should not happen.
        start = nodes[0];
      }
    }

    // Walk the path.
    std::vector<std::size_t> chain;
    chain.reserve(nodes.size());
    std::size_t prev = static_cast<std::size_t>(-1);
    std::size_t cur = start;

    for (;;) {
      chain.push_back(cur);
      // choose next neighbor not equal to prev
      std::optional<std::size_t> next;
      if (g.has_adjacency()) {
        const auto& adj = g.adjacency();
        for (std::size_t nb : adj[cur]) {
          if (nb == prev) continue;
          next = nb;
          break;
        }
      } else {
        // EdgeList backend: fallback O(E) neighbor search (rare for large; most topology graphs are adjacency).
        std::size_t found = 0;
        std::size_t nb = 0;
        g.for_each_edge([&](pilots::alg::graph::NodeId u, pilots::alg::graph::NodeId v) {
          if (u == cur && v != prev) { nb = v; ++found; }
          else if (v == cur && u != prev) { nb = u; ++found; }
        });
        if (found > 0) next = nb;
      }

      if (!next.has_value()) {
        break;
      }
      prev = cur;
      cur = *next;
      if (chain.size() > nodes.size()) {
        throw std::runtime_error("ChainIndex: path walk exceeded component size (graph not linear?)");
      }
    }

    if (chain.size() != nodes.size()) {
      if (strict) {
        throw std::runtime_error("ChainIndex: linear component walk did not visit all nodes");
      }
      continue;
    }

    out.chain_id.push_back(static_cast<std::int64_t>(cid + 1));
    out.chains.push_back(std::move(chain));
  }

  return out;
}

} // namespace pilots::alg::index
