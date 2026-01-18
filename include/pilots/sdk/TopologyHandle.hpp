#pragma once

#include <stdexcept>
#include <string>

#include "pilots/core/SystemContext.hpp"
#include "pilots/topology/Topology.hpp"

namespace pilots::sdk {

class TopologyHandle {
public:
  explicit TopologyHandle(std::string measure_name, const pilots::SystemContext& ctx)
      : measure_(std::move(measure_name)), topo_(ctx.topology) {}

  bool has() const { return topo_ != nullptr; }

  const pilots::Topology& require() const {
    if (!topo_) throw std::runtime_error(qualified_("requires topology, but none was loaded"));
    return *topo_;
  }

  void require_section(const std::string& sec_name) const {
    const auto& t = require();
    if (sec_name == "masses") {
      if (!t.has_masses()) throw std::runtime_error(qualified_("requires topology masses, but masses were not loaded"));
      return;
    }
    if (sec_name == "bonds") {
      if (!t.has_bonds()) throw std::runtime_error(qualified_("requires topology bonds, but bonds were not loaded"));
      return;
    }
    if (sec_name == "angles") {
      if (!t.has_angles()) throw std::runtime_error(qualified_("requires topology angles, but angles were not loaded"));
      return;
    }
    if (sec_name == "dihedrals") {
      if (!t.has_dihedrals()) throw std::runtime_error(qualified_("requires topology dihedrals, but dihedrals were not loaded"));
      return;
    }
    if (sec_name == "impropers") {
      if (!t.has_impropers()) throw std::runtime_error(qualified_("requires topology impropers, but impropers were not loaded"));
      return;
    }
    throw std::runtime_error(qualified_("unknown topology section '" + sec_name + "'"));
  }

private:
  std::string measure_;
  const pilots::Topology* topo_ = nullptr;

  std::string qualified_(const std::string& msg) const {
    if (measure_.empty()) return std::string("TopologyHandle: ") + msg;
    return std::string("TopologyHandle[") + measure_ + "]: " + msg;
  }
};

} // namespace pilots::sdk
