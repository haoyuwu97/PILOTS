#pragma once

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include "pilots/select/SelectionProvider.hpp"

namespace pilots::sdk {

// A thin wrapper around SelectionProvider::get_combined_view...
//
// It centralizes:
//   - canonical combined selection key string (A=...;T=...;combine=...)
//   - a stable place to hang audit access
class SelectionHandle {
public:
  SelectionHandle(std::string measure_name,
                 pilots::SelectionProvider* provider,
                 std::string group_ref,
                 std::string topo_group_ref,
                 std::string combine_expr)
      : measure_(std::move(measure_name)),
        provider_(provider),
        group_ref_(std::move(group_ref)),
        topo_ref_(std::move(topo_group_ref)),
        combine_(std::move(combine_expr)) {
    if (!provider_) throw std::runtime_error(qualified_("selection_provider is null"));
  }

  const std::string& group_ref() const { return group_ref_; }
  const std::string& topo_group_ref() const { return topo_ref_; }
  const std::string& combine_expr() const { return combine_; }

  // Canonical key that matches SelectionProvider's internal cache key format.
  std::string key() const {
    const std::string g = group_ref_.empty() ? std::string("all") : group_ref_;
    const std::string t = topo_ref_.empty() ? std::string("all") : topo_ref_;
    const std::string c = combine_.empty() ? std::string("A&T") : combine_;
    std::string k;
    k.reserve(g.size() + t.size() + c.size() + 64);
    k += "A="; k += g;
    k += ";T="; k += t;
    k += ";combine="; k += c;
    return k;
  }

  bool is_dynamic_spec() const {
    const std::string g = group_ref_.empty() ? std::string("all") : group_ref_;
    const std::string t = topo_ref_.empty() ? std::string("all") : topo_ref_;
    return provider_->is_dynamic_spec(g, t);
  }

  pilots::SelectionView view(const pilots::Frame& f, std::size_t frame_index) {
    return provider_->get_combined_view(f, frame_index, group_ref_, topo_ref_, combine_);
  }

  // Best-effort audit lookup. For dynamic specs this returns size stats collected so far.
  std::optional<pilots::SelectionProvider::SelectionAuditInfo> audit_info() const {
    const std::string k = key();
    for (const auto& a : provider_->audit()) {
      if (a.key == k) return a;
    }
    return std::nullopt;
  }

private:
  std::string measure_;
  pilots::SelectionProvider* provider_ = nullptr;
  std::string group_ref_;
  std::string topo_ref_;
  std::string combine_;

  std::string qualified_(const std::string& msg) const {
    if (measure_.empty()) return std::string("SelectionHandle: ") + msg;
    return std::string("SelectionHandle[") + measure_ + "]: " + msg;
  }
};

} // namespace pilots::sdk
