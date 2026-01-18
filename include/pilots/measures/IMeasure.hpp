#pragma once

#include <cstddef>
#include <istream>
#include <ostream>
#include <string>

#include "pilots/core/Frame.hpp"
#include "pilots/core/SystemContext.hpp"
#include "pilots/output/ResultsIndex.hpp"

namespace pilots {

// Minimal Measure interface (Step2):
// - A single trajectory-read loop calls on_frame(...) for all enabled measures.
// - Each measure writes its own output(s) and manages its own internal buffering.
//
// P0 extension:
// - flush_partial(): allow safe intermediate output in follow/online mode
// - save_state()/load_state(): enable checkpoint/resume for long HPC jobs
// - on_start(first_frame, ctx): inject topology/system context without custom I/O
class IMeasure {
public:
  virtual ~IMeasure() = default;

  // Measure type (e.g. "msd", "isf").
  virtual std::string type() const = 0;

  // Instance name from config (e.g. section [measure.foo] -> "foo").
  // Must be unique within a run.
  virtual std::string instance_name() const = 0;

  // Backward-compatible alias: prefer instance_name().
  virtual std::string name() const { return instance_name(); }

  // Called once, after the first frame has been read and group selections are available.
  virtual void on_start(const Frame& first_frame) = 0;

  // Context-aware startup (default: ignore ctx and call legacy on_start).
  virtual void on_start(const Frame& first_frame, const SystemContext& ctx) {
    (void)ctx;
    on_start(first_frame);
  }

  // Called for every frame in the input trajectory.
  virtual void on_frame(const Frame& frame, std::size_t frame_index) = 0;

  // Called once after the input is exhausted (or when the run terminates).
  virtual void finalize() = 0;

  // Optional: write intermediate output that is safe to call multiple times.
  virtual void flush_partial() {}

  // Optional: checkpoint interface (default = stateless).
  virtual void save_state(std::ostream&) const {}
  virtual void load_state(std::istream&) {}

  // Machine-readable descriptor for results indexing (e.g., results.json).
  virtual output::MeasureDescriptor describe() const = 0;
};

} // namespace pilots
