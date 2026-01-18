#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>

#include "pilots/correlate/CorrelationSeries.hpp"
#include "pilots/correlate/ICorrelatorT6.hpp"

namespace pilots {

// Placeholder FFT-based correlator (tensor6). Not implemented in this step.
// We keep the class to freeze the API shape and allow future drop-in implementation.

template <typename PairOp>
class FFTCorrelatorStubT6 final : public ICorrelatorT6 {
public:
  struct Options {
    double dt = 0.0;
  };

  FFTCorrelatorStubT6(std::size_t, Options opt, PairOp)
  : opt_(opt) {
    if (!(opt_.dt > 0.0)) {
      throw std::runtime_error("FFTCorrelatorStubT6: dt must be > 0");
    }
  }

  void start() override {
    throw std::runtime_error("FFT correlator is a stub in this step (planned). Use exact or multitau.");
  }

  void push(std::int64_t,
            std::span<const double>,
            std::span<const double>,
            std::span<const double>,
            std::span<const double>,
            std::span<const double>,
            std::span<const double>) override {
    throw std::runtime_error("FFT correlator is a stub in this step (planned). Use exact or multitau.");
  }

  CorrelationSeriesT6 finalize() override {
    throw std::runtime_error("FFT correlator is a stub in this step (planned). Use exact or multitau.");
  }

private:
  Options opt_;
};

} // namespace pilots
