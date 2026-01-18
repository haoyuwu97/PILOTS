#pragma once

#include <chrono>

namespace pilots {

class WallTimer {
public:
  using clock = std::chrono::steady_clock;

  WallTimer() : t0_(clock::now()) {}

  void reset() { t0_ = clock::now(); }

  double elapsed_seconds() const {
    const auto dt = clock::now() - t0_;
    return std::chrono::duration_cast<std::chrono::duration<double>>(dt).count();
  }

private:
  clock::time_point t0_;
};

// Accumulates elapsed time into a referenced double.
class ScopedTimer {
public:
  explicit ScopedTimer(double* accumulator) : acc_(accumulator) {
    if (acc_) t_.reset();
  }
  ~ScopedTimer() {
    if (acc_) *acc_ += t_.elapsed_seconds();
  }

  ScopedTimer(const ScopedTimer&) = delete;
  ScopedTimer& operator=(const ScopedTimer&) = delete;

private:
  double* acc_ = nullptr;
  WallTimer t_;
};

} // namespace pilots
