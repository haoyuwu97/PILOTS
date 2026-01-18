#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include "pilots/alg/graph/EdgeList.hpp"
#include "pilots/alg/graph/GraphView.hpp"
#include "pilots/alg/mapping/BeadGraphBuilder.hpp"
#include "pilots/core/SystemContext.hpp"
#include "pilots/topology/BondGraph.hpp"
#include "pilots/topology/Topology.hpp"

namespace pilots::sdk {

enum class GraphSource {
  TopologyGraph,
  BondGraph,
  BeadGraph,
};

inline std::string graph_source_name(GraphSource s) {
  switch (s) {
    case GraphSource::TopologyGraph: return "topology_graph";
    case GraphSource::BondGraph: return "bond_graph";
    case GraphSource::BeadGraph: return "bead_graph";
  }
  return "topology_graph";
}

// GraphHandle: a stable holder for a graph view (topology_graph/bond_graph/bead_graph).
//
// - For adjacency-based graphs, the view points into Topology/BondGraph storage.
// - For bead_graph, the EdgeList is owned by the handle to keep the view valid.
class GraphHandle {
public:
  GraphHandle() = default;

  static GraphHandle from_system(std::string measure_name,
                                 const pilots::SystemContext& ctx,
                                 GraphSource src) {
    GraphHandle h;
    h.measure_ = std::move(measure_name);
    h.source_ = src;

    if (src == GraphSource::TopologyGraph) {
      if (!ctx.topology) throw std::runtime_error(h.qualified_("requires topology_graph but topology is null"));
      if (!ctx.topology->has_bonds()) {
        throw std::runtime_error(h.qualified_("requires topology_graph but topology bonds are not loaded"));
      }
      h.gv_ = pilots::alg::graph::GraphView(ctx.topology->adjacency);
      return h;
    }

    if (src == GraphSource::BondGraph) {
      if (!ctx.bond_graph) throw std::runtime_error(h.qualified_("requires bond_graph but ctx.bond_graph is null"));
      h.gv_ = pilots::alg::graph::GraphView(ctx.bond_graph->adjacency);
      return h;
    }

    if (src == GraphSource::BeadGraph) {
      if (!ctx.mapping) throw std::runtime_error(h.qualified_("requires bead_graph but mapping is disabled"));
      if (!ctx.topology || !ctx.topology->has_bonds()) {
        throw std::runtime_error(h.qualified_("requires bead_graph but topology bonds are not loaded"));
      }
      h.owned_edges_ = pilots::alg::mapping::build_bead_graph_from_bonds(*ctx.topology, *ctx.mapping);
      h.gv_ = pilots::alg::graph::GraphView(*h.owned_edges_);
      return h;
    }

    throw std::runtime_error(h.qualified_("unknown graph source"));
  }

  // Convenience: choose a reasonable default without imposing backbone semantics.
  // Prefer bead_graph when mapping is enabled, otherwise bond_graph if available,
  // otherwise topology_graph.
  static GraphHandle auto_from_system(std::string measure_name,
                                      const pilots::SystemContext& ctx) {
    if (ctx.mapping && ctx.topology && ctx.topology->has_bonds()) {
      return from_system(std::move(measure_name), ctx, GraphSource::BeadGraph);
    }
    if (ctx.bond_graph) {
      return from_system(std::move(measure_name), ctx, GraphSource::BondGraph);
    }
    return from_system(std::move(measure_name), ctx, GraphSource::TopologyGraph);
  }

  GraphSource source() const { return source_; }
  std::string source_name() const { return graph_source_name(source_); }

  const pilots::alg::graph::GraphView& view() const { return gv_; }

private:
  std::string measure_;
  GraphSource source_ = GraphSource::TopologyGraph;
  std::optional<pilots::alg::graph::EdgeList> owned_edges_;
  pilots::alg::graph::GraphView gv_;

  std::string qualified_(const std::string& msg) const {
    if (measure_.empty()) return std::string("GraphHandle: ") + msg;
    return std::string("GraphHandle[") + measure_ + "]: " + msg;
  }
};

} // namespace pilots::sdk
