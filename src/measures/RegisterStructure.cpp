#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if PILOTS_HAS_OPENMP
#include <omp.h>
#endif

#include "MeasureCommon.hpp"
#include "pilots/measures/IMeasure.hpp"
#include "pilots/measures/MeasureRegistry.hpp"
#include "pilots/select/SelectionView.hpp"
#include "pilots/util/AtomicFile.hpp"

namespace fs = std::filesystem;

namespace pilots {
namespace {

using measure_ext::box_volume;
using measure_ext::count_diag_dims;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::get_static_group_view;
using measure_ext::join_doubles;
using measure_ext::parse_double_list;
using measure_ext::resolve_path;
using measure_ext::same_index_set;
using measure_ext::shell_volume_3d;
using measure_ext::x_unit_for_axis;

constexpr double kPi = 3.141592653589793238462643383279502884;

inline double spherical_bessel_j0(double x) {
  if (std::abs(x) < 1e-14) return 1.0;
  return std::sin(x) / x;
}

struct RangeOptions {
  std::int64_t frame_start = 0;
  std::int64_t frame_end = -1;
  bool dry_run = false;
};

struct LagOptions {
  std::int64_t max_lag_frames = -1;
  std::int64_t lag_stride = 1;
  std::int64_t origin_stride = 1;
};

struct StoredFrame {
  std::int64_t timestep = 0;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
};

inline std::vector<double> q_values_from_config(const IniConfig& cfg,
                                                const std::string& section) {
  if (cfg.has_key(section, "q_values")) {
    return parse_double_list(cfg, section, "q_values");
  }
  if (cfg.has_key(section, "q_min") && cfg.has_key(section, "q_max") && cfg.has_key(section, "q_step")) {
    const double qmin = cfg.get_double(section, "q_min");
    const double qmax = cfg.get_double(section, "q_max");
    const double dq = cfg.get_double(section, "q_step");
    if (!(dq > 0.0)) throw std::runtime_error(section + ": q_step must be > 0");
    if (!(qmax >= qmin)) throw std::runtime_error(section + ": q_max must be >= q_min");
    std::vector<double> q;
    for (double x = qmin; x <= qmax + 0.5 * dq; x += dq) q.push_back(x);
    if (q.empty()) throw std::runtime_error(section + ": q grid is empty");
    return q;
  }
  throw std::runtime_error(section + ": provide q_values=... or q_min/q_max/q_step");
}

inline LagOptions lag_options_from_config(const IniConfig& cfg,
                                          const std::string& section) {
  LagOptions opt;
  opt.max_lag_frames = cfg.get_int64(section, "max_lag_frames", std::optional<std::int64_t>(-1));
  opt.lag_stride = cfg.get_int64(section, "lag_stride", std::optional<std::int64_t>(1));
  opt.origin_stride = cfg.get_int64(section, "origin_stride", std::optional<std::int64_t>(1));
  if (opt.max_lag_frames == 0) throw std::runtime_error(section + ": max_lag_frames must be -1 or >= 1");
  if (opt.max_lag_frames < -1) throw std::runtime_error(section + ": max_lag_frames must be -1 or >= 1");
  if (opt.lag_stride <= 0) throw std::runtime_error(section + ": lag_stride must be >= 1");
  if (opt.origin_stride <= 0) throw std::runtime_error(section + ": origin_stride must be >= 1");
  return opt;
}

inline std::vector<std::size_t> build_lag_list(std::size_t nframes,
                                               const LagOptions& opt,
                                               bool include_zero = true) {
  std::vector<std::size_t> out;
  if (nframes == 0) return out;
  const std::size_t max_lag = (opt.max_lag_frames < 0)
      ? (nframes - 1)
      : std::min<std::size_t>(static_cast<std::size_t>(opt.max_lag_frames), nframes - 1);
  std::size_t start = include_zero ? 0 : 1;
  for (std::size_t lag = start; lag <= max_lag; lag += static_cast<std::size_t>(opt.lag_stride)) {
    out.push_back(lag);
  }
  if (include_zero && out.empty()) out.push_back(0);
  return out;
}

inline double sem_from_sums(double sum, double sumsq, std::size_t n) {
  if (n < 2) return 0.0;
  const double mean = sum / static_cast<double>(n);
  const double sample_var = std::max(0.0, (sumsq - static_cast<double>(n) * mean * mean) / static_cast<double>(n - 1));
  return std::sqrt(sample_var / static_cast<double>(n));
}

inline double mean_dtimestep_for_lag(const std::vector<StoredFrame>& frames,
                                     std::size_t lag,
                                     std::size_t origin_stride) {
  if (frames.empty()) return 0.0;
  if (lag >= frames.size()) return 0.0;
  double sum = 0.0;
  std::size_t n = 0;
  for (std::size_t o = 0; o + lag < frames.size(); o += origin_stride) {
    sum += static_cast<double>(frames[o + lag].timestep - frames[o].timestep);
    ++n;
  }
  return (n > 0) ? (sum / static_cast<double>(n)) : 0.0;
}

inline bool frame_in_range(std::size_t frame_index, const RangeOptions& opt) {
  const std::int64_t fi = static_cast<std::int64_t>(frame_index);
  if (fi < opt.frame_start) return false;
  if (opt.frame_end >= 0 && fi > opt.frame_end) return false;
  return true;
}

inline void fill_centered_positions(const Frame& frame,
                                    SelectionView sel,
                                    bool remove_drift,
                                    SelectionView drift_sel,
                                    std::vector<double>& outx,
                                    std::vector<double>& outy,
                                    std::vector<double>& outz) {
  const auto xu = frame.require_dfield("xu");
  const auto yu = frame.require_dfield("yu");
  const auto zu = frame.require_dfield("zu");

  double comx = 0.0, comy = 0.0, comz = 0.0;
  if (remove_drift) {
    double sx = 0.0, sy = 0.0, sz = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sx,sy,sz)
#endif
    for (std::size_t k = 0; k < drift_sel.idx.size(); ++k) {
      const std::size_t i = drift_sel.idx[k];
      sx += xu[i];
      sy += yu[i];
      sz += zu[i];
    }
    const double inv = 1.0 / static_cast<double>(drift_sel.idx.size());
    comx = sx * inv;
    comy = sy * inv;
    comz = sz * inv;
  }

  outx.resize(sel.idx.size());
  outy.resize(sel.idx.size());
  outz.resize(sel.idx.size());
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
  for (std::size_t p = 0; p < sel.idx.size(); ++p) {
    const std::size_t i = sel.idx[p];
    outx[p] = xu[i] - comx;
    outy[p] = yu[i] - comy;
    outz[p] = zu[i] - comz;
  }
}

inline void append_common_selection_caps(const IniConfig& cfg,
                                         const std::string& section,
                                         MeasureCapabilities& caps,
                                         bool add_drift) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  if (cfg.has_key(section, "group_b")) {
    caps.group_refs.push_back(cfg.get_string(section, "group_b"));
  }
  if (add_drift && cfg.get_bool(section, "remove_drift", std::optional<bool>(true))) {
    caps.group_refs.push_back(cfg.get_string(section, "drift_group", std::optional<std::string>("all")));
  }
}

class RDFMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    double dr = 0.0;
    double r_max = 0.0;
  };

  RDFMeasure(std::string instance_name,
             std::string output_path,
             SelectionView sel_a,
             SelectionView sel_b,
             Options opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(opt) {
    sel_a_name_owned_ = std::string(sel_a.name);
    sel_b_name_owned_ = std::string(sel_b.name);
    sel_a_ = SelectionView{sel_a_name_owned_, sel_a.idx};
    sel_b_ = SelectionView{sel_b_name_owned_, sel_b.idx};
    same_sel_ = same_index_set(sel_a_, sel_b_);

    if (opt_.range.frame_start < 0) throw std::runtime_error("rdf: frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error("rdf: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.dr > 0.0)) throw std::runtime_error("rdf: dr must be > 0");
    if (!(opt_.r_max > 0.0)) throw std::runtime_error("rdf: r_max must be > 0");
    if (sel_a_.idx.empty() || sel_b_.idx.empty()) throw std::runtime_error("rdf: selection is empty");
    if (same_sel_ && sel_a_.idx.size() < 2) throw std::runtime_error("rdf: same-selection RDF requires at least 2 atoms");

    const std::size_t nb = static_cast<std::size_t>(std::ceil(opt_.r_max / opt_.dr));
    g_sum_.assign(nb, 0.0);
  }

  std::string type() const override { return "rdf"; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_a_.name);
    md.n_selected = sel_a_.idx.size();

    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.x_axis = "r";
    od.x_unit = "trajectory_length";
    od.columns = {"r_lo", "r_hi", "r", "g_r", "coordination", "n_frames"};
    md.outputs.push_back(std::move(od));

    md.params["selection_b"] = std::string(sel_b_.name);
    md.params["same_selection"] = same_sel_ ? "1" : "0";
    md.params["dr"] = dstr(opt_.dr);
    md.params["r_max"] = dstr(opt_.r_max);
    md.params["frame_start"] = std::to_string(opt_.range.frame_start);
    md.params["frame_end"] = std::to_string(opt_.range.frame_end);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;

    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    const double V = box_volume(frame.box);
    const std::size_t nA = sel_a_.idx.size();
    const std::size_t nB = sel_b_.idx.size();

    const double rhoB = same_sel_
        ? ((nA > 1) ? (static_cast<double>(nA - 1) / V) : 0.0)
        : (static_cast<double>(nB) / V);
    if (!(rhoB > 0.0)) return;

    std::vector<double> hist(g_sum_.size(), 0.0);
    if (same_sel_) {
      for (std::size_t ia = 0; ia < nA; ++ia) {
        const std::size_t i = sel_a_.idx[ia];
        for (std::size_t ja = ia + 1; ja < nA; ++ja) {
          const std::size_t j = sel_a_.idx[ja];
          const auto d = frame.box.min_image_displacement(xu[i], yu[i], zu[i], xu[j], yu[j], zu[j]);
          const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
          if (r >= opt_.r_max) continue;
          const std::size_t b = static_cast<std::size_t>(r / opt_.dr);
          if (b < hist.size()) hist[b] += 2.0;
        }
      }
    } else {
      for (std::size_t ia = 0; ia < nA; ++ia) {
        const std::size_t i = sel_a_.idx[ia];
        for (std::size_t jb = 0; jb < nB; ++jb) {
          const std::size_t j = sel_b_.idx[jb];
          if (i == j) continue;
          const auto d = frame.box.min_image_displacement(xu[i], yu[i], zu[i], xu[j], yu[j], zu[j]);
          const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
          if (r >= opt_.r_max) continue;
          const std::size_t b = static_cast<std::size_t>(r / opt_.dr);
          if (b < hist.size()) hist[b] += 1.0;
        }
      }
    }

    const double denom_pref = static_cast<double>(nA) * rhoB;
    for (std::size_t b = 0; b < hist.size(); ++b) {
      const double r_lo = static_cast<double>(b) * opt_.dr;
      const double r_hi = std::min(opt_.r_max, r_lo + opt_.dr);
      const double shell = shell_volume_3d(r_lo, r_hi);
      if (shell > 0.0) {
        g_sum_[b] += hist[b] / (denom_pref * shell);
      }
    }
    rho_b_sum_ += rhoB;
    ++n_frames_;
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      const double mean_rho_b = (n_frames_ > 0) ? (rho_b_sum_ / static_cast<double>(n_frames_)) : 0.0;
      double coord = 0.0;
      for (std::size_t b = 0; b < g_sum_.size(); ++b) {
        const double r_lo = static_cast<double>(b) * opt_.dr;
        const double r_hi = std::min(opt_.r_max, r_lo + opt_.dr);
        const double r = 0.5 * (r_lo + r_hi);
        const double g = (n_frames_ > 0) ? (g_sum_[b] / static_cast<double>(n_frames_)) : 0.0;
        coord += mean_rho_b * g * shell_volume_3d(r_lo, r_hi);
        ofs << std::setprecision(17) << r_lo << ' '
            << std::setprecision(17) << r_hi << ' '
            << std::setprecision(17) << r << ' '
            << std::setprecision(17) << g << ' '
            << std::setprecision(17) << coord << ' '
            << n_frames_ << '\n';
      }
    });
  }

  void finalize() override {
    if (opt_.range.dry_run) return;
    flush_partial();
  }

private:
  std::string instance_name_;
  std::string output_path_;
  std::string sel_a_name_owned_;
  std::string sel_b_name_owned_;
  SelectionView sel_a_;
  SelectionView sel_b_;
  Options opt_;
  bool same_sel_ = false;
  bool started_ = false;
  std::vector<double> g_sum_;
  double rho_b_sum_ = 0.0;
  std::size_t n_frames_ = 0;

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: radial distribution function g(r)\n";
    ofs << "# selection_a: " << sel_a_.name << " (n=" << sel_a_.idx.size() << ")\n";
    ofs << "# selection_b: " << sel_b_.name << " (n=" << sel_b_.idx.size() << ")\n";
    ofs << "# same_selection: " << (same_sel_ ? 1 : 0) << "\n";
    ofs << "# dr: " << std::setprecision(17) << opt_.dr << "\n";
    ofs << "# r_max: " << std::setprecision(17) << opt_.r_max << "\n";
    ofs << "# frame_range: [" << opt_.range.frame_start << ", " << opt_.range.frame_end << "]\n";
    ofs << "# columns: r_lo  r_hi  r  g_r  coordination  n_frames\n";
  }
};

class StaticStructureFactorMeasure final : public IMeasure {
public:
  StaticStructureFactorMeasure(std::string instance_name,
                               std::string output_path,
                               SelectionView sel,
                               RangeOptions range,
                               std::vector<double> q_values)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        range_(range),
        q_values_(std::move(q_values)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    if (range_.frame_start < 0) throw std::runtime_error("static_structure_factor: frame_start must be >= 0");
    if (range_.frame_end >= 0 && range_.frame_end < range_.frame_start) {
      throw std::runtime_error("static_structure_factor: frame_end must be -1 or >= frame_start");
    }
    if (sel_.idx.empty()) throw std::runtime_error("static_structure_factor: selection is empty");
    if (q_values_.empty()) throw std::runtime_error("static_structure_factor: q_values is empty");
    sq_sum_.assign(q_values_.size(), 0.0);
    sq_sumsq_.assign(q_values_.size(), 0.0);
  }

  std::string type() const override { return "static_structure_factor"; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();

    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.x_axis = "q";
    od.x_unit = "1/trajectory_length";
    od.columns = {"q", "sq", "sem", "n_frames"};
    md.outputs.push_back(std::move(od));
    md.params["q_values"] = join_doubles(q_values_);
    md.params["frame_start"] = std::to_string(range_.frame_start);
    md.params["frame_end"] = std::to_string(range_.frame_end);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!range_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, range_)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    const std::size_t N = sel_.idx.size();

    std::vector<double> sq(q_values_.size(), 1.0);
    for (std::size_t ia = 0; ia < N; ++ia) {
      const std::size_t i = sel_.idx[ia];
      for (std::size_t ja = ia + 1; ja < N; ++ja) {
        const std::size_t j = sel_.idx[ja];
        const auto d = frame.box.min_image_displacement(xu[i], yu[i], zu[i], xu[j], yu[j], zu[j]);
        const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        for (std::size_t iq = 0; iq < q_values_.size(); ++iq) {
          sq[iq] += (2.0 / static_cast<double>(N)) * spherical_bessel_j0(q_values_[iq] * r);
        }
      }
    }

    for (std::size_t iq = 0; iq < q_values_.size(); ++iq) {
      sq_sum_[iq] += sq[iq];
      sq_sumsq_[iq] += sq[iq] * sq[iq];
    }
    ++n_frames_;
  }

  void flush_partial() override {
    if (range_.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: static structure factor S(q) from direct Debye sum\n";
      ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
      ofs << "# frame_range: [" << range_.frame_start << ", " << range_.frame_end << "]\n";
      ofs << "# q_values: " << join_doubles(q_values_) << "\n";
      ofs << "# columns: q  sq  sem  n_frames\n";
      for (std::size_t iq = 0; iq < q_values_.size(); ++iq) {
        const double mean = (n_frames_ > 0) ? (sq_sum_[iq] / static_cast<double>(n_frames_)) : 0.0;
        const double sem = sem_from_sums(sq_sum_[iq], sq_sumsq_[iq], n_frames_);
        ofs << std::setprecision(17) << q_values_[iq] << ' '
            << std::setprecision(17) << mean << ' '
            << std::setprecision(17) << sem << ' '
            << n_frames_ << '\n';
      }
    });
  }

  void finalize() override {
    if (range_.dry_run) return;
    flush_partial();
  }

private:
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  SelectionView sel_;
  RangeOptions range_;
  std::vector<double> q_values_;
  std::vector<double> sq_sum_;
  std::vector<double> sq_sumsq_;
  std::size_t n_frames_ = 0;
  bool started_ = false;
};

class DynamicStructureFactorMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    LagOptions lag;
    bool remove_drift = true;
  };

  DynamicStructureFactorMeasure(std::string measure_type,
                                std::string instance_name,
                                std::string output_path,
                                SelectionView sel,
                                SelectionView drift_sel,
                                double dt,
                                Options opt,
                                std::vector<double> q_values)
      : measure_type_(std::move(measure_type)),
        instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        dt_(dt),
        opt_(opt),
        q_values_(std::move(q_values)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};

    if (!(dt_ > 0.0)) throw std::runtime_error(measure_type_ + ": dt must be > 0");
    if (opt_.range.frame_start < 0) throw std::runtime_error(measure_type_ + ": frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error(measure_type_ + ": frame_end must be -1 or >= frame_start");
    }
    if (sel_.idx.empty()) throw std::runtime_error(measure_type_ + ": selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) throw std::runtime_error(measure_type_ + ": drift selection is empty");
    if (q_values_.empty()) throw std::runtime_error(measure_type_ + ": q_values is empty");
  }

  std::string type() const override { return measure_type_; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();

    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.x_axis = "lag";
    od.x_unit = "frames";
    od.columns = {"lag", "time", "q", "f_qt", "fs_qt", "fd_qt", "sem_f_qt", "count_origins", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    md.params["dt"] = dstr(dt_);
    md.params["q_values"] = join_doubles(q_values_);
    md.params["frame_start"] = std::to_string(opt_.range.frame_start);
    md.params["frame_end"] = std::to_string(opt_.range.frame_end);
    md.params["lag_stride"] = std::to_string(opt_.lag.lag_stride);
    md.params["max_lag_frames"] = std::to_string(opt_.lag.max_lag_frames);
    md.params["origin_stride"] = std::to_string(opt_.lag.origin_stride);
    md.params["remove_drift"] = opt_.remove_drift ? "1" : "0";
    md.params["drift_group"] = std::string(drift_sel_.name);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    StoredFrame f;
    f.timestep = frame.timestep;
    fill_centered_positions(frame, sel_, opt_.remove_drift, drift_sel_, f.x, f.y, f.z);
    frames_.push_back(std::move(f));
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      const auto lags = build_lag_list(frames_.size(), opt_.lag, true);
      const std::size_t N = sel_.idx.size();
      for (const double q : q_values_) {
        for (const std::size_t lag : lags) {
          double sum = 0.0, sumsq = 0.0;
          double sum_self = 0.0;
          std::size_t n_origins = 0;
          for (std::size_t o = 0; o + lag < frames_.size(); o += static_cast<std::size_t>(opt_.lag.origin_stride)) {
            const auto& f0 = frames_[o];
            const auto& f1 = frames_[o + lag];
            double val = 0.0;
            double val_self = 0.0;
            for (std::size_t i = 0; i < N; ++i) {
              const double dxi = f1.x[i] - f0.x[i];
              const double dyi = f1.y[i] - f0.y[i];
              const double dzi = f1.z[i] - f0.z[i];
              val_self += spherical_bessel_j0(q * std::sqrt(dxi * dxi + dyi * dyi + dzi * dzi));
              for (std::size_t j = 0; j < N; ++j) {
                const double dx = f1.x[j] - f0.x[i];
                const double dy = f1.y[j] - f0.y[i];
                const double dz = f1.z[j] - f0.z[i];
                val += spherical_bessel_j0(q * std::sqrt(dx * dx + dy * dy + dz * dz));
              }
            }
            val /= static_cast<double>(N);
            val_self /= static_cast<double>(N);
            sum += val;
            sumsq += val * val;
            sum_self += val_self;
            ++n_origins;
          }
          const double mean = (n_origins > 0) ? (sum / static_cast<double>(n_origins)) : 0.0;
          const double mean_self = (n_origins > 0) ? (sum_self / static_cast<double>(n_origins)) : 0.0;
          const double mean_dist = mean - mean_self;
          const double sem = sem_from_sums(sum, sumsq, n_origins);
          const double mean_dts = mean_dtimestep_for_lag(frames_, lag, static_cast<std::size_t>(opt_.lag.origin_stride));
          ofs << lag << ' '
              << std::setprecision(17) << (mean_dts * dt_) << ' '
              << std::setprecision(17) << q << ' '
              << std::setprecision(17) << mean << ' '
              << std::setprecision(17) << mean_self << ' '
              << std::setprecision(17) << mean_dist << ' '
              << std::setprecision(17) << sem << ' '
              << n_origins << ' '
              << std::setprecision(17) << mean_dts << '\n';
        }
      }
    });
  }

  void finalize() override {
    if (opt_.range.dry_run) return;
    flush_partial();
  }

private:
  std::string measure_type_;
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  double dt_ = 0.0;
  Options opt_;
  std::vector<double> q_values_;
  std::vector<StoredFrame> frames_;
  bool started_ = false;

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: coherent intermediate scattering function / time-domain dynamic structure factor\n";
    ofs << "# type: " << measure_type_ << "\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# remove_drift: " << (opt_.remove_drift ? 1 : 0) << " (drift_group=" << drift_sel_.name << ")\n";
    ofs << "# frame_range: [" << opt_.range.frame_start << ", " << opt_.range.frame_end << "]\n";
    ofs << "# q_values: " << join_doubles(q_values_) << "\n";
    ofs << "# lag_stride: " << opt_.lag.lag_stride << ", max_lag_frames: " << opt_.lag.max_lag_frames
        << ", origin_stride: " << opt_.lag.origin_stride << "\n";
    ofs << "# dt: " << std::setprecision(17) << dt_ << " (seconds per timestep)\n";
    ofs << "# columns: lag  time  q  f_qt  fs_qt  fd_qt  sem_f_qt  count_origins  mean_dtimestep\n";
  }
};

class VanHoveMeasure final : public IMeasure {
public:
  enum class Mode {
    Self,
    Distinct,
  };

  struct Options {
    RangeOptions range;
    LagOptions lag;
    bool remove_drift = true;
    double dr = 0.0;
    double r_max = 0.0;
    Mode mode = Mode::Self;
  };

  VanHoveMeasure(std::string instance_name,
                 std::string output_path,
                 SelectionView sel,
                 SelectionView drift_sel,
                 double dt,
                 Options opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        dt_(dt),
        opt_(opt) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};

    if (!(dt_ > 0.0)) throw std::runtime_error(type() + ": dt must be > 0");
    if (opt_.range.frame_start < 0) throw std::runtime_error(type() + ": frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error(type() + ": frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.dr > 0.0)) throw std::runtime_error(type() + ": dr must be > 0");
    if (!(opt_.r_max > 0.0)) throw std::runtime_error(type() + ": r_max must be > 0");
    if (sel_.idx.empty()) throw std::runtime_error(type() + ": selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) throw std::runtime_error(type() + ": drift selection is empty");
  }

  std::string type() const override {
    return (opt_.mode == Mode::Self) ? "van_hove_self" : "van_hove_distinct";
  }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();

    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.x_axis = "lag";
    od.x_unit = "frames";
    od.columns = {"lag", "time", "r_lo", "r_hi", "r", (opt_.mode == Mode::Self) ? "gs_rt" : "gd_rt", "count_origins", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    md.params["dt"] = dstr(dt_);
    md.params["dr"] = dstr(opt_.dr);
    md.params["r_max"] = dstr(opt_.r_max);
    md.params["frame_start"] = std::to_string(opt_.range.frame_start);
    md.params["frame_end"] = std::to_string(opt_.range.frame_end);
    md.params["lag_stride"] = std::to_string(opt_.lag.lag_stride);
    md.params["max_lag_frames"] = std::to_string(opt_.lag.max_lag_frames);
    md.params["origin_stride"] = std::to_string(opt_.lag.origin_stride);
    md.params["remove_drift"] = opt_.remove_drift ? "1" : "0";
    md.params["drift_group"] = std::string(drift_sel_.name);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    StoredFrame f;
    f.timestep = frame.timestep;
    fill_centered_positions(frame, sel_, opt_.remove_drift, drift_sel_, f.x, f.y, f.z);
    frames_.push_back(std::move(f));
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      const auto lags = build_lag_list(frames_.size(), opt_.lag, true);
      const std::size_t N = sel_.idx.size();
      const std::size_t nb = static_cast<std::size_t>(std::ceil(opt_.r_max / opt_.dr));
      for (const std::size_t lag : lags) {
        std::vector<double> hist(nb, 0.0);
        std::size_t n_origins = 0;
        for (std::size_t o = 0; o + lag < frames_.size(); o += static_cast<std::size_t>(opt_.lag.origin_stride)) {
          const auto& f0 = frames_[o];
          const auto& f1 = frames_[o + lag];
          if (opt_.mode == Mode::Self) {
            for (std::size_t i = 0; i < N; ++i) {
              const double dx = f1.x[i] - f0.x[i];
              const double dy = f1.y[i] - f0.y[i];
              const double dz = f1.z[i] - f0.z[i];
              const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
              if (r >= opt_.r_max) continue;
              const std::size_t b = static_cast<std::size_t>(r / opt_.dr);
              if (b < hist.size()) hist[b] += 1.0;
            }
          } else {
            for (std::size_t i = 0; i < N; ++i) {
              for (std::size_t j = 0; j < N; ++j) {
                if (i == j) continue;
                const double dx = f1.x[j] - f0.x[i];
                const double dy = f1.y[j] - f0.y[i];
                const double dz = f1.z[j] - f0.z[i];
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r >= opt_.r_max) continue;
                const std::size_t b = static_cast<std::size_t>(r / opt_.dr);
                if (b < hist.size()) hist[b] += 1.0;
              }
            }
          }
          ++n_origins;
        }

        const double mean_dts = mean_dtimestep_for_lag(frames_, lag, static_cast<std::size_t>(opt_.lag.origin_stride));
        for (std::size_t b = 0; b < hist.size(); ++b) {
          const double r_lo = static_cast<double>(b) * opt_.dr;
          const double r_hi = std::min(opt_.r_max, r_lo + opt_.dr);
          const double r = 0.5 * (r_lo + r_hi);
          const double shell = shell_volume_3d(r_lo, r_hi);
          const double denom = (opt_.mode == Mode::Self)
              ? (static_cast<double>(n_origins) * static_cast<double>(N) * shell)
              : (static_cast<double>(n_origins) * static_cast<double>(N) * shell);
          const double g = (denom > 0.0) ? (hist[b] / denom) : 0.0;
          ofs << lag << ' '
              << std::setprecision(17) << (mean_dts * dt_) << ' '
              << std::setprecision(17) << r_lo << ' '
              << std::setprecision(17) << r_hi << ' '
              << std::setprecision(17) << r << ' '
              << std::setprecision(17) << g << ' '
              << n_origins << ' '
              << std::setprecision(17) << mean_dts << '\n';
        }
      }
    });
  }

  void finalize() override {
    if (opt_.range.dry_run) return;
    flush_partial();
  }

private:
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  double dt_ = 0.0;
  Options opt_;
  std::vector<StoredFrame> frames_;
  bool started_ = false;

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: " << type() << " radial distribution in space and time\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# remove_drift: " << (opt_.remove_drift ? 1 : 0) << " (drift_group=" << drift_sel_.name << ")\n";
    ofs << "# frame_range: [" << opt_.range.frame_start << ", " << opt_.range.frame_end << "]\n";
    ofs << "# lag_stride: " << opt_.lag.lag_stride << ", max_lag_frames: " << opt_.lag.max_lag_frames
        << ", origin_stride: " << opt_.lag.origin_stride << "\n";
    ofs << "# dr: " << std::setprecision(17) << opt_.dr << ", r_max: " << std::setprecision(17) << opt_.r_max << "\n";
    ofs << "# dt: " << std::setprecision(17) << dt_ << " (seconds per timestep)\n";
    ofs << "# columns: lag  time  r_lo  r_hi  r  " << ((opt_.mode == Mode::Self) ? "gs_rt" : "gd_rt")
        << "  count_origins  mean_dtimestep\n";
  }
};

MeasureCapabilities rdf_caps(const IniConfig& cfg,
                             const std::string& section,
                             const std::string& instance,
                             const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_common_selection_caps(cfg, section, caps, false);
  return caps;
}

std::unique_ptr<IMeasure> rdf_create(const IniConfig& cfg,
                                     const std::string& section,
                                     const std::string& instance,
                                     const MeasureBuildEnv& env,
                                     const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("rdf factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;

  const std::string group_a = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo_a = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb_a = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));

  const std::string group_b = cfg.get_string(section, "group_b", std::optional<std::string>(group_a));
  const std::string topo_b = cfg.get_string(section, "topo_group_b", std::optional<std::string>(topo_a));
  const std::string comb_b = cfg.get_string(section, "combine_b", std::optional<std::string>(comb_a));

  SelectionView sel_a = get_static_combined_view(*env.selection_provider, frame0, group_a, topo_a, comb_a, "rdf");
  SelectionView sel_b = get_static_combined_view(*env.selection_provider, frame0, group_b, topo_b, comb_b, "rdf");

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("rdf.dat"))).lexically_normal();

  RDFMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.dr = cfg.get_double(section, "dr");
  opt.r_max = cfg.get_double(section, "r_max");
  return std::make_unique<RDFMeasure>(instance, out.string(), sel_a, sel_b, opt);
}

MeasureCapabilities sq_static_caps(const IniConfig& cfg,
                                   const std::string& section,
                                   const std::string& instance,
                                   const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_common_selection_caps(cfg, section, caps, false);
  return caps;
}

std::unique_ptr<IMeasure> sq_static_create(const IniConfig& cfg,
                                           const std::string& section,
                                           const std::string& instance,
                                           const MeasureBuildEnv& env,
                                           const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("static_structure_factor factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;

  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "static_structure_factor");

  const auto q_values = q_values_from_config(cfg, section);
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("static_structure_factor.dat"))).lexically_normal();

  RangeOptions range;
  range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  range.dry_run = env.dry_run;

  return std::make_unique<StaticStructureFactorMeasure>(instance, out.string(), sel, range, q_values);
}

MeasureCapabilities dynsf_caps(const IniConfig& cfg,
                               const std::string& section,
                               const std::string& instance,
                               const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_common_selection_caps(cfg, section, caps, true);
  caps.requires_identity_consistent = true;
  return caps;
}

std::unique_ptr<IMeasure> dynsf_create(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& instance,
                                       const MeasureBuildEnv& env,
                                       const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("dynamic_structure_factor factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;

  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));

  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "dynamic_structure_factor");
  SelectionView drift_sel = remove_drift
      ? get_static_group_view(*env.selection_provider, frame0, drift_group, "dynamic_structure_factor")
      : env.selection_provider->get_combined_view(frame0, 0, "all", "all", "A");

  const auto q_values = q_values_from_config(cfg, section);
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const std::string default_name = (cfg.get_string(section, "type") == "coherent_isf") ? "coherent_isf.dat" : "dynamic_structure_factor.dat";
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>(default_name))).lexically_normal();

  DynamicStructureFactorMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.lag = lag_options_from_config(cfg, section);
  opt.remove_drift = remove_drift;

  return std::make_unique<DynamicStructureFactorMeasure>(
      cfg.get_string(section, "type"), instance, out.string(), sel, drift_sel, env.dt, opt, q_values);
}

MeasureCapabilities vanhove_caps(const IniConfig& cfg,
                                 const std::string& section,
                                 const std::string& instance,
                                 const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_common_selection_caps(cfg, section, caps, true);
  caps.requires_identity_consistent = true;
  return caps;
}

std::unique_ptr<IMeasure> vanhove_create(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env,
                                         const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("van_hove factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;

  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));

  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "van_hove");
  SelectionView drift_sel = remove_drift
      ? get_static_group_view(*env.selection_provider, frame0, drift_group, "van_hove")
      : env.selection_provider->get_combined_view(frame0, 0, "all", "all", "A");

  const std::string type_s = cfg.get_string(section, "type");
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const std::string default_name = (type_s == "van_hove_self") ? "van_hove_self.dat" : "van_hove_distinct.dat";
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>(default_name))).lexically_normal();

  VanHoveMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.lag = lag_options_from_config(cfg, section);
  opt.remove_drift = remove_drift;
  opt.dr = cfg.get_double(section, "dr");
  opt.r_max = cfg.get_double(section, "r_max");
  opt.mode = (type_s == "van_hove_self") ? VanHoveMeasure::Mode::Self : VanHoveMeasure::Mode::Distinct;

  return std::make_unique<VanHoveMeasure>(instance, out.string(), sel, drift_sel, env.dt, opt);
}

static MeasureRegistrar g_register_rdf("rdf", &rdf_caps, &rdf_create);
static MeasureRegistrar g_register_sq_static("static_structure_factor", &sq_static_caps, &sq_static_create);
static MeasureRegistrar g_register_sq_static_alias("sq_static", &sq_static_caps, &sq_static_create);
static MeasureRegistrar g_register_dynsf("dynamic_structure_factor", &dynsf_caps, &dynsf_create);
static MeasureRegistrar g_register_dynsf_alias("coherent_isf", &dynsf_caps, &dynsf_create);
static MeasureRegistrar g_register_vhs("van_hove_self", &vanhove_caps, &vanhove_create);
static MeasureRegistrar g_register_vhd("van_hove_distinct", &vanhove_caps, &vanhove_create);

} // namespace
} // namespace pilots
