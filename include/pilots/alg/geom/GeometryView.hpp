#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <stdexcept>

#include "pilots/core/Box.hpp"

namespace pilots::alg::geom {

// A light, non-owning view over particle positions + simulation box.
//
// Positions are expected to be unwrapped coordinates (xu/yu/zu). For distance
// computations under PBC we use Box::min_image_displacement, which supports
// orthorhombic and triclinic boxes.
struct GeometryView {
  std::span<const double> x;
  std::span<const double> y;
  std::span<const double> z;
  const pilots::Box* box = nullptr;

  GeometryView() = default;

  GeometryView(std::span<const double> xu,
               std::span<const double> yu,
               std::span<const double> zu,
               const pilots::Box& b)
      : x(xu), y(yu), z(zu), box(&b) {
    if (x.size() != y.size() || x.size() != z.size()) {
      throw std::runtime_error("GeometryView: position spans must have same length");
    }
  }

  std::size_t size() const { return x.size(); }

  std::array<double,3> position(std::size_t i) const {
    if (i >= size()) throw std::runtime_error("GeometryView: index out of range");
    return {x[i], y[i], z[i]};
  }

  // Minimum-image displacement rj - ri.
  std::array<double,3> min_image_displacement(std::size_t i, std::size_t j) const {
    if (!box) throw std::runtime_error("GeometryView: box is null");
    if (i >= size() || j >= size()) throw std::runtime_error("GeometryView: index out of range");
    return box->min_image_displacement(x[i], y[i], z[i], x[j], y[j], z[j]);
  }

  double min_image_distance_sq(std::size_t i, std::size_t j) const {
    const auto d = min_image_displacement(i, j);
    return d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
  }

  // Wrap a point into the primary box.
  std::array<double,3> wrap_position(const std::array<double,3>& r) const {
    if (!box) throw std::runtime_error("GeometryView: box is null");
    return box->wrap(r[0], r[1], r[2]);
  }
};

} // namespace pilots::alg::geom
