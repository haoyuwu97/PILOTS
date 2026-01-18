#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pilots/alg/graph/Components.hpp"
#include "pilots/alg/graph/Descriptors.hpp"
#include "pilots/alg/graph/GraphView.hpp"

namespace pilots::alg::polymer {

struct ComponentSummary {
  std::size_t id = 0;
  std::size_t size = 0;
  std::size_t n_ends = 0;          // degree==1 within component
  std::size_t n_branch_points = 0; // degree>=3 within component
  std::int64_t cycle_rank = 0;     // E - V + 1 (connected component)
  std::string tag;                // linear|ring|branched|network|other
};

struct PolymerClassification {
  std::string graph_source; // topology_graph|bond_graph|bead_graph
  std::vector<ComponentSummary> components;
  std::vector<std::string> tags; // global tags: linear/ring/branched/network/mixed
  std::string hint;

  std::size_t n_components() const { return components.size(); }
  std::size_t largest_component_size() const {
    std::size_t m = 0;
    for (const auto& c : components) m = std::max(m, c.size);
    return m;
  }
};

inline std::string to_lower(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  return s;
}

inline std::string classify_component_tag(const ComponentSummary& c) {
  // Note: this is a coarse, graph-theoretic tag. Interpretation depends on the graph abstraction.
  if (c.size <= 1) return "other";
  if (c.cycle_rank == 0) {
    if (c.n_branch_points == 0 && c.n_ends == 2) return "linear";
    if (c.n_branch_points > 0) return "branched";
    return "other";
  }
  // cycle_rank > 0
  if (c.n_branch_points == 0 && c.n_ends == 0 && c.cycle_rank == 1) return "ring";
  return "network";
}

inline PolymerClassification classify_graph(const pilots::alg::graph::GraphView& g,
                                           std::string graph_source,
                                           std::optional<std::string> model_scale_name = std::nullopt) {
  PolymerClassification out;
  out.graph_source = std::move(graph_source);

  const auto comps = pilots::alg::graph::compute_components(g);
  const auto& deg = g.degree();

  // Count edges per component.
  std::vector<std::size_t> edges_in_comp(comps.component_sizes.size(), 0);
  g.for_each_edge([&](pilots::alg::graph::NodeId u, pilots::alg::graph::NodeId v) {
    const std::size_t cu = comps.component_id[u];
    const std::size_t cv = comps.component_id[v];
    if (cu == cv) {
      ++edges_in_comp[cu];
    }
  });

  // Count ends/branch points per component in one pass.
  std::vector<std::size_t> ends_in_comp(comps.component_sizes.size(), 0);
  std::vector<std::size_t> branch_in_comp(comps.component_sizes.size(), 0);
  for (std::size_t i = 0; i < comps.component_id.size(); ++i) {
    const std::size_t cid = comps.component_id[i];
    const std::size_t d = deg[i];
    if (d == 1) ++ends_in_comp[cid];
    if (d >= 3) ++branch_in_comp[cid];
  }

  out.components.reserve(comps.component_sizes.size());
  for (std::size_t cid = 0; cid < comps.component_sizes.size(); ++cid) {
    ComponentSummary c;
    c.id = cid;
    c.size = comps.component_sizes[cid];

    c.n_ends = ends_in_comp[cid];
    c.n_branch_points = branch_in_comp[cid];

    const std::size_t V = c.size;
    const std::size_t E = edges_in_comp[cid];
    // For a connected component, cycle_rank = E - V + 1.
    c.cycle_rank = static_cast<std::int64_t>(E) - static_cast<std::int64_t>(V) + 1;

    c.tag = classify_component_tag(c);
    out.components.push_back(std::move(c));
  }

  // Global tags: if mixed component types exist, mark mixed.
  std::unordered_map<std::string, std::size_t> cnt;
  for (const auto& c : out.components) {
    ++cnt[c.tag];
  }

  // Order of preference: if only one meaningful tag, use it; else mixed.
  std::vector<std::string> keys;
  keys.reserve(cnt.size());
  for (const auto& kv : cnt) {
    if (kv.first == "other") continue;
    keys.push_back(kv.first);
  }
  std::sort(keys.begin(), keys.end());

  if (keys.empty()) {
    out.tags = {"other"};
  } else if (keys.size() == 1) {
    out.tags = {keys[0]};
  } else {
    out.tags = {"mixed"};
    // keep component-type list for debugging
    for (const auto& k : keys) out.tags.push_back(k);
  }

  // Applicability hint (heuristic, best-effort).
  // We avoid encoding heavy semantics; the intent is to warn against over-interpreting
  // atom-level covalent graphs in AA/UA.
  if (model_scale_name.has_value()) {
    const std::string ms = to_lower(*model_scale_name);
    if ((ms == "aa" || ms == "ua") && (out.graph_source == "topology_graph" || out.graph_source == "bond_graph")) {
      out.hint =
          "Hint: model.scale=" + *model_scale_name + " and graph_source=" + out.graph_source +
          ". Degree/branching on an atom-level covalent graph often reflects chemistry/valence, not polymer architecture. "
          "Consider using mapping (bead_graph) or a filtered bond_graph to represent the architectural backbone.";
    }
  }

  return out;
}

} // namespace pilots::alg::polymer
