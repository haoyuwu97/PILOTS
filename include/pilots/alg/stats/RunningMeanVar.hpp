#pragma once

#include <cmath>
#include <cstddef>

namespace pilots::alg::stats {

// Welford online mean/variance.
// Tracks population variance by default; sample variance available too.
class RunningMeanVar {
public:
  void reset() {
    n_ = 0;
    mean_ = 0.0;
    m2_ = 0.0;
  }

  void add(double x) {
    ++n_;
    const double delta = x - mean_;
    mean_ += delta / static_cast<double>(n_);
    const double delta2 = x - mean_;
    m2_ += delta * delta2;
  }

  std::size_t count() const { return n_; }
  double mean() const { return mean_; }

  // Population variance.
  double var_pop() const { return (n_ > 0) ? (m2_ / static_cast<double>(n_)) : 0.0; }

  // Unbiased sample variance.
  double var_sample() const { return (n_ > 1) ? (m2_ / static_cast<double>(n_ - 1)) : 0.0; }

  double std_pop() const { return std::sqrt(var_pop()); }
  double std_sample() const { return std::sqrt(var_sample()); }

private:
  std::size_t n_ = 0;
  double mean_ = 0.0;
  double m2_ = 0.0;
};

} // namespace pilots::alg::stats
