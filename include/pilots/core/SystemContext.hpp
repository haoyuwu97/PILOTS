#pragma once

#include <cstddef>

namespace pilots {

struct Topology;
struct MoleculeIndex;
struct BondGraph;
namespace alg::mapping { struct BeadMapping; }

// SystemContext carries trajectory-invariant information into measures.
//
// P0: provide topology injection so polymer/topology-dependent measures can
// be written without custom I/O.
struct SystemContext {
  const Topology* topology = nullptr; // may be null if no topology configured
  const MoleculeIndex* molecules = nullptr; // may be null if 'mol' not present or not built
  const BondGraph* bond_graph = nullptr; // may be null if topology has no bonds or graph not constructed
  const alg::mapping::BeadMapping* mapping = nullptr; // may be null if mapping is disabled

  bool has_topology() const { return topology != nullptr; }
  bool has_molecules() const { return molecules != nullptr; }
  bool has_bond_graph() const { return bond_graph != nullptr; }
  bool has_mapping() const { return mapping != nullptr; }
};

} // namespace pilots
