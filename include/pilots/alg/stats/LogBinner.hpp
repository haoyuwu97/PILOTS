#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace pilots::alg::stats {

// Logarithmic binner for positive values.
//
// Defines bin edges spaced uniformly in log space between [min_value, max_value].
// Values <=0 are invalid.
class LogBinner {
public:
  LogBinner() = default;

  LogBinner(double min_value, double max_value, std::size_t n_bins)
      : min_(min_value), max_(max_value), n_bins_(n_bins) {
    if (!(max_ > min_)) throw std::runtime_error("LogBinner: max must be > min");
    if (min_ <= 0.0) throw std::runtime_error("LogBinner: min must be > 0");
    if (n_bins_ == 0) throw std::runtime_error("LogBinner: n_bins must be > 0");
    log_min_ = std::log(min_);
    log_max_ = std::log(max_);
    inv_dlog_ = static_cast<double>(n_bins_) / (log_max_ - log_min_);
  }

  std::size_t n_bins() const { return n_bins_; }

  // Returns bin index in [0, n_bins-1], or n_bins if out of range.
  std::size_t bin_index(double x) const {
    if (x < min_ || x >= max_) return n_bins_;
    const double lx = std::log(x);
    const std::size_t i = static_cast<std::size_t>((lx - log_min_) * inv_dlog_);
    return std::min(i, n_bins_ - 1);
  }

  // Bin edges (size n_bins+1).
  std::vector<double> edges() const {
    std::vector<double> e(n_bins_ + 1);
    for (std::size_t i = 0; i <= n_bins_; ++i) {
      const double t = static_cast<double>(i) / static_cast<double>(n_bins_);
      e[i] = std::exp(log_min_ + t * (log_max_ - log_min_));
    }
    return e;
  }

private:
  double min_ = 1.0;
  double max_ = 10.0;
  std::size_t n_bins_ = 0;
  double log_min_ = 0.0;
  double log_max_ = 0.0;
  double inv_dlog_ = 1.0;
};

} // namespace pilots::alg::stats
