#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace pilots::alg::mapping {

struct BeadFrameView {
  std::size_t nbeads = 0;
  std::span<const std::int64_t> bead_id;
  std::span<const double> xu;
  std::span<const double> yu;
  std::span<const double> zu;

  // Optional (may be empty).
  std::span<const double> bead_mass;
  std::span<const int> bead_type;
};

// Owning container for a bead-level "frame".
struct BeadFrame {
  std::vector<std::int64_t> bead_id;
  std::vector<double> xu;
  std::vector<double> yu;
  std::vector<double> zu;
  std::vector<double> bead_mass;
  std::vector<int> bead_type;

  std::size_t nbeads() const { return bead_id.size(); }

  BeadFrameView view() const {
    return BeadFrameView{
        bead_id.size(),
        std::span<const std::int64_t>(bead_id.data(), bead_id.size()),
        std::span<const double>(xu.data(), xu.size()),
        std::span<const double>(yu.data(), yu.size()),
        std::span<const double>(zu.data(), zu.size()),
        std::span<const double>(bead_mass.data(), bead_mass.size()),
        std::span<const int>(bead_type.data(), bead_type.size())};
  }
};

} // namespace pilots::alg::mapping
