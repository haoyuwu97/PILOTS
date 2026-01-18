#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "pilots/alg/stats/RunningMeanVar.hpp"

namespace pilots::alg::stats {

// Simple block averaging accumulator.
//
// Feed a sequence of scalar samples x_t. The class groups samples into blocks of
// fixed size B. Each completed block yields a block mean. Statistics of block means
// can be used to estimate SEM under weak correlation assumptions.
class BlockAverager {
public:
  explicit BlockAverager(std::size_t block_size = 1) : block_size_(block_size) {
    if (block_size_ == 0) throw std::runtime_error("BlockAverager: block_size must be > 0");
  }

  void reset() {
    n_in_block_ = 0;
    sum_in_block_ = 0.0;
    blocks_.clear();
    stats_.reset();
  }

  std::size_t block_size() const { return block_size_; }
  std::size_t n_blocks() const { return blocks_.size(); }

  void add(double x) {
    sum_in_block_ += x;
    ++n_in_block_;
    if (n_in_block_ == block_size_) {
      const double mean = sum_in_block_ / static_cast<double>(block_size_);
      blocks_.push_back(mean);
      stats_.add(mean);
      n_in_block_ = 0;
      sum_in_block_ = 0.0;
    }
  }

  // Returns block means (completed blocks only).
  const std::vector<double>& blocks() const { return blocks_; }

  // Mean of block means.
  double mean() const { return stats_.mean(); }

  // Unbiased variance of block means.
  double var() const { return stats_.var_sample(); }

  // Standard error of the mean estimated from block means.
  double sem() const {
    const std::size_t n = stats_.count();
    if (n == 0) return 0.0;
    return std::sqrt(stats_.var_sample() / static_cast<double>(n));
  }

private:
  std::size_t block_size_ = 1;
  std::size_t n_in_block_ = 0;
  double sum_in_block_ = 0.0;
  std::vector<double> blocks_;
  RunningMeanVar stats_;
};

} // namespace pilots::alg::stats
