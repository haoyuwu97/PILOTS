#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>

#include "pilots/core/Frame.hpp"

namespace pilots::sdk {

// ExtraFieldHandle<T>: resolve an extra per-atom field once (on_start), then
// access it by field_id on every frame (hot path without string lookup).
//
// Current PILOTS extra storage is double-only, so T is constrained to double.
// The template is kept for future dtype expansion.
template <typename T>
class ExtraFieldHandle;

template <>
class ExtraFieldHandle<double> {
public:
  explicit ExtraFieldHandle(std::string field_name, std::string measure_name = {})
      : name_(std::move(field_name)), measure_(std::move(measure_name)) {}

  const std::string& name() const { return name_; }

  // Resolve field_id. Call in on_start(first_frame).
  void resolve(const pilots::Frame& first_frame) {
    fid_ = first_frame.require_extra_id(name_);
    resolved_ = true;
  }

  bool resolved() const { return resolved_; }
  std::size_t field_id() const {
    if (!resolved_) throw std::runtime_error(qualified_("field_id requested before resolve()"));
    return fid_;
  }

  // Get a span for a frame. Throws if the field is missing in this frame's ATOMS header.
  std::span<const double> get(const pilots::Frame& f) const {
    if (!resolved_) {
      throw std::runtime_error(qualified_("get() called before resolve(); call resolve(first_frame) in on_start"));
    }
    return f.require_dfield_by_id(fid_);
  }

private:
  std::string name_;
  std::string measure_;
  std::size_t fid_ = static_cast<std::size_t>(-1);
  bool resolved_ = false;

  std::string qualified_(const std::string& msg) const {
    if (measure_.empty()) return std::string("ExtraFieldHandle: ") + msg + " (field='" + name_ + "')";
    return std::string("ExtraFieldHandle[") + measure_ + "]: " + msg + " (field='" + name_ + "')";
  }
};

} // namespace pilots::sdk
