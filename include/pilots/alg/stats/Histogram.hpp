#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pilots::alg::stats {

class Histogram1D {
public:
  Histogram1D() = default;

  Histogram1D(double min_value, double max_value, std::size_t n_bins)
      : min_(min_value), max_(max_value), bins_(n_bins, 0.0) {
    if (!(max_ > min_)) throw std::runtime_error("Histogram1D: max must be > min");
    if (n_bins == 0) throw std::runtime_error("Histogram1D: n_bins must be > 0");
    inv_dx_ = static_cast<double>(n_bins) / (max_ - min_);
  }

  void reset() {
    std::fill(bins_.begin(), bins_.end(), 0.0);
    underflow_ = 0.0;
    overflow_ = 0.0;
  }

  std::size_t n_bins() const { return bins_.size(); }
  double min() const { return min_; }
  double max() const { return max_; }

  void add(double x, double w = 1.0) {
    if (x < min_) { underflow_ += w; return; }
    if (x >= max_) { overflow_ += w; return; }
    const std::size_t i = static_cast<std::size_t>((x - min_) * inv_dx_);
    bins_[std::min(i, bins_.size() - 1)] += w;
  }

  const std::vector<double>& bins() const { return bins_; }
  double underflow() const { return underflow_; }
  double overflow() const { return overflow_; }

  // Bin center coordinate.
  double center(std::size_t i) const {
    if (i >= bins_.size()) throw std::runtime_error("Histogram1D: center index out of range");
    const double dx = (max_ - min_) / static_cast<double>(bins_.size());
    return min_ + (static_cast<double>(i) + 0.5) * dx;
  }

private:
  double min_ = 0.0;
  double max_ = 1.0;
  double inv_dx_ = 1.0;
  std::vector<double> bins_;
  double underflow_ = 0.0;
  double overflow_ = 0.0;
};

class Histogram2D {
public:
  Histogram2D() = default;

  Histogram2D(double xmin, double xmax, std::size_t nx,
              double ymin, double ymax, std::size_t ny)
      : xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax),
        nx_(nx), ny_(ny), bins_(nx * ny, 0.0) {
    if (!(xmax_ > xmin_)) throw std::runtime_error("Histogram2D: xmax must be > xmin");
    if (!(ymax_ > ymin_)) throw std::runtime_error("Histogram2D: ymax must be > ymin");
    if (nx_ == 0 || ny_ == 0) throw std::runtime_error("Histogram2D: bins must be > 0");
    inv_dx_ = static_cast<double>(nx_) / (xmax_ - xmin_);
    inv_dy_ = static_cast<double>(ny_) / (ymax_ - ymin_);
  }

  void reset() {
    std::fill(bins_.begin(), bins_.end(), 0.0);
    underflow_ = 0.0;
    overflow_ = 0.0;
  }

  std::pair<std::size_t,std::size_t> shape() const { return {nx_, ny_}; }

  void add(double x, double y, double w = 1.0) {
    if (x < xmin_ || y < ymin_) { underflow_ += w; return; }
    if (x >= xmax_ || y >= ymax_) { overflow_ += w; return; }
    const std::size_t ix = static_cast<std::size_t>((x - xmin_) * inv_dx_);
    const std::size_t iy = static_cast<std::size_t>((y - ymin_) * inv_dy_);
    const std::size_t iix = std::min(ix, nx_ - 1);
    const std::size_t iiy = std::min(iy, ny_ - 1);
    bins_[iix + nx_ * iiy] += w;
  }

  const std::vector<double>& bins() const { return bins_; }
  double underflow() const { return underflow_; }
  double overflow() const { return overflow_; }

private:
  double xmin_ = 0.0, xmax_ = 1.0;
  double ymin_ = 0.0, ymax_ = 1.0;
  std::size_t nx_ = 0, ny_ = 0;
  double inv_dx_ = 1.0, inv_dy_ = 1.0;
  std::vector<double> bins_;
  double underflow_ = 0.0;
  double overflow_ = 0.0;
};

} // namespace pilots::alg::stats
