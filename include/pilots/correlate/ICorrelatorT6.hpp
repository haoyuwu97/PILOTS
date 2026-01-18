#pragma once

#include <cstdint>
#include <istream>
#include <ostream>
#include <span>
#include <stdexcept>

#include "pilots/correlate/CorrelationSeries.hpp"

namespace pilots {

// Type-erased correlator interface for 6-channel tensor signals.
// Concrete correlators (exact / multi-tau / FFT) implement this interface.
//
// P0 extension:
// - snapshot(): non-destructive readout for online/follow mode
// - save_state()/load_state(): checkpoint/resume
class ICorrelatorT6 {
public:
  virtual ~ICorrelatorT6() = default;

  virtual void start() = 0;

  // Push one sample in chronological order.
  // All spans must have the same length (nsel).
  virtual void push(std::int64_t timestep,
                    std::span<const double> xx,
                    std::span<const double> yy,
                    std::span<const double> zz,
                    std::span<const double> xy,
                    std::span<const double> xz,
                    std::span<const double> yz) = 0;

  // Finalize and return the tensor-6 correlation series.
  // NOTE: In the current implementations, finalize() is non-destructive and may be called multiple times.
  virtual CorrelationSeriesT6 finalize() = 0;

  // Online snapshot (default: calls finalize()).
  virtual CorrelationSeriesT6 snapshot() const {
    // finalize() is not const; current correlators are logically const here.
    return const_cast<ICorrelatorT6*>(this)->finalize();
  }

  // Checkpoint/resume hooks.
  // Default implementations throw: concrete correlators should override if they claim to support checkpointing.
  virtual void save_state(std::ostream&) const {
    throw std::runtime_error("ICorrelatorT6: save_state not implemented for this correlator");
  }

  virtual void load_state(std::istream&) {
    throw std::runtime_error("ICorrelatorT6: load_state not implemented for this correlator");
  }
};

} // namespace pilots
