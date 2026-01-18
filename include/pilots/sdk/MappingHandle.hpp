#pragma once

#include <cstdint>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

#include "pilots/alg/mapping/BeadMapping.hpp"
#include "pilots/core/SystemContext.hpp"

namespace pilots::sdk {

class MappingHandle {
public:
  explicit MappingHandle(std::string measure_name, const pilots::SystemContext& ctx)
      : measure_(std::move(measure_name)), map_(ctx.mapping) {}

  bool enabled() const { return map_ != nullptr && map_->spec.enabled(); }
  bool has() const { return map_ != nullptr; }

  const pilots::alg::mapping::BeadMapping& require() const {
    if (!map_) throw std::runtime_error(qualified_("requires mapping, but mapping is disabled"));
    return *map_;
  }

  // Convenience: compute a bead-level frame for the current atom frame.
  pilots::alg::mapping::BeadFrame compute_bead_frame(const pilots::Frame& f) const {
    return require().compute_bead_frame(f);
  }

  std::uint64_t spec_hash_fnv1a64() const {
    return require().spec.spec_hash_fnv1a64();
  }

  std::string spec_hash_hex() const {
    std::ostringstream oss;
    oss << std::hex << std::setw(16) << std::setfill('0') << spec_hash_fnv1a64();
    return oss.str();
  }

  const pilots::alg::mapping::MappingSpec& spec() const {
    return require().spec;
  }

private:
  std::string measure_;
  const pilots::alg::mapping::BeadMapping* map_ = nullptr;

  std::string qualified_(const std::string& msg) const {
    if (measure_.empty()) return std::string("MappingHandle: ") + msg;
    return std::string("MappingHandle[") + measure_ + "]: " + msg;
  }
};

} // namespace pilots::sdk
