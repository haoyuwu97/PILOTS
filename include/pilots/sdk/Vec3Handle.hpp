#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>

#include "pilots/core/Frame.hpp"

namespace pilots::sdk {

struct Vec3View {
  std::span<const double> x;
  std::span<const double> y;
  std::span<const double> z;
};

// Vec3Handle: canonical access to unwrapped coordinates.
//
// PILOTS normalizes LAMMPS dumps into Frame::xu/yu/zu, regardless of whether
// the dump provides x/y/z with image flags or direct xu/yu/zu.
//
// This handle mainly exists to:
//   - centralize error messages
//   - provide a consistent view type (Vec3View)
//   - keep future extensions (e.g., bead coordinates) consistent
class Vec3Handle {
public:
  explicit Vec3Handle(std::string measure_name = {})
      : measure_(std::move(measure_name)) {}

  Vec3View get(const pilots::Frame& f) const {
    if (f.xu.size() != f.natoms || f.yu.size() != f.natoms || f.zu.size() != f.natoms) {
      throw std::runtime_error(qualified_("positions xu/yu/zu are missing or wrong length"));
    }
    return Vec3View{std::span<const double>(f.xu.data(), f.xu.size()),
                    std::span<const double>(f.yu.data(), f.yu.size()),
                    std::span<const double>(f.zu.data(), f.zu.size())};
  }

private:
  std::string measure_;

  std::string qualified_(const std::string& msg) const {
    if (measure_.empty()) return std::string("Vec3Handle: ") + msg;
    return std::string("Vec3Handle[") + measure_ + "]: " + msg;
  }
};

} // namespace pilots::sdk
