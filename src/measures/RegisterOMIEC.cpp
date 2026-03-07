#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "MeasureCommon.hpp"
#include "pilots/measures/IMeasure.hpp"
#include "pilots/measures/MeasureRegistry.hpp"
#include "pilots/select/SelectionView.hpp"
#include "pilots/util/AtomicFile.hpp"

namespace fs = std::filesystem;

namespace pilots {
namespace {

using measure_ext::Axis1D;
using measure_ext::axis1d_name;
using measure_ext::axis_length;
using measure_ext::box_volume;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::mass_by_atom_from_config;
using measure_ext::parse_axis1d;
using measure_ext::primary_axis_coord;
using measure_ext::resolve_path;
using measure_ext::same_index_set;

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

inline bool frame_in_range(std::size_t frame_index, const RangeOptions& opt) {
  const std::int64_t fi = static_cast<std::int64_t>(frame_index);
  if (fi < opt.frame_start) return false;
  if (opt.frame_end >= 0 && fi > opt.frame_end) return false;
  return true;
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

inline double mean_dtimestep_for_lag(const std::vector<std::int64_t>& timesteps,
                                     std::size_t lag,
                                     std::size_t origin_stride) {
  if (timesteps.empty() || lag >= timesteps.size()) return 0.0;
  double sum = 0.0;
  std::size_t n = 0;
  for (std::size_t o = 0; o + lag < timesteps.size(); o += origin_stride) {
    sum += static_cast<double>(timesteps[o + lag] - timesteps[o]);
    ++n;
  }
  return (n > 0) ? (sum / static_cast<double>(n)) : 0.0;
}

inline double sem_from_sums(double sum, double sumsq, std::size_t n) {
  if (n < 2) return 0.0;
  const double mean = sum / static_cast<double>(n);
  const double sample_var = std::max(0.0, (sumsq - static_cast<double>(n) * mean * mean) / static_cast<double>(n - 1));
  return std::sqrt(sample_var / static_cast<double>(n));
}

enum class ProfileMode {
  Number,
  Charge,
  Mass,
};

inline ProfileMode parse_profile_mode(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  if (s == "number" || s == "count") return ProfileMode::Number;
  if (s == "charge" || s == "q") return ProfileMode::Charge;
  if (s == "mass") return ProfileMode::Mass;
  throw std::runtime_error("invalid profile mode: '" + s + "' (use number|charge|mass)");
}

inline const char* profile_mode_name(ProfileMode m) {
  switch (m) {
    case ProfileMode::Number: return "number";
    case ProfileMode::Charge: return "charge";
    case ProfileMode::Mass: return "mass";
  }
  return "number";
}

class Profile1DMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    Axis1D axis = Axis1D::Z;
    std::size_t n_bins = 0;
    ProfileMode mode = ProfileMode::Number;
    std::string charge_field = "q";
  };

  Profile1DMeasure(std::string instance_name,
                   std::string output_path,
                   SelectionView sel,
                   Options opt,
                   std::vector<double> mass_by_atom)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        mass_by_atom_(std::move(mass_by_atom)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    if (opt_.range.frame_start < 0) throw std::runtime_error("profile1d: frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error("profile1d: frame_end must be -1 or >= frame_start");
    }
    if (opt_.n_bins == 0) throw std::runtime_error("profile1d: n_bins must be >= 1");
    if (sel_.idx.empty()) throw std::runtime_error("profile1d: selection is empty");
    weight_sum_.assign(opt_.n_bins, 0.0);
    density_sum_.assign(opt_.n_bins, 0.0);
  }

  std::string type() const override { return "profile1d"; }
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
    od.x_axis = "bin";
    od.x_unit = "fractional_coordinate";
    od.columns = {"bin", "coord_lo_frac", "coord_hi_frac", "coord_center_frac", "coord_center_mean", "mean_weight", "mean_density", "n_frames"};
    md.outputs.push_back(std::move(od));

    md.params["axis"] = axis1d_name(opt_.axis);
    md.params["n_bins"] = std::to_string(opt_.n_bins);
    md.params["mode"] = profile_mode_name(opt_.mode);
    md.params["charge_field"] = opt_.charge_field;
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
    const auto q = (opt_.mode == ProfileMode::Charge)
        ? frame.require_dfield(opt_.charge_field)
        : std::span<const double>();
    const double vbin = box_volume(frame.box) / static_cast<double>(opt_.n_bins);

    std::vector<double> local_weight(opt_.n_bins, 0.0);
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t i = sel_.idx[p];
      double w = 1.0;
      if (opt_.mode == ProfileMode::Charge) {
        w = q[i];
      } else if (opt_.mode == ProfileMode::Mass) {
        w = mass_by_atom_.at(i);
      }
      const double s = primary_axis_coord(frame.box, xu[i], yu[i], zu[i], opt_.axis) / axis_length(frame.box, opt_.axis);
      std::size_t b = static_cast<std::size_t>(std::floor(s * static_cast<double>(opt_.n_bins)));
      if (b >= opt_.n_bins) b = opt_.n_bins - 1;
      local_weight[b] += w;
    }

    const double coord_scale = axis_length(frame.box, opt_.axis);
    mean_axis_length_ += coord_scale;
    for (std::size_t b = 0; b < opt_.n_bins; ++b) {
      weight_sum_[b] += local_weight[b];
      density_sum_[b] += local_weight[b] / vbin;
    }
    ++n_frames_;
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: 1D profile for OMIEC/interfacial analysis\n";
      ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
      ofs << "# axis: " << axis1d_name(opt_.axis) << ", n_bins: " << opt_.n_bins << "\n";
      ofs << "# mode: " << profile_mode_name(opt_.mode) << "\n";
      if (opt_.mode == ProfileMode::Charge) ofs << "# charge_field: " << opt_.charge_field << "\n";
      ofs << "# frame_range: [" << opt_.range.frame_start << ", " << opt_.range.frame_end << "]\n";
      ofs << "# columns: bin  coord_lo_frac  coord_hi_frac  coord_center_frac  coord_center_mean  mean_weight  mean_density  n_frames\n";
      const double mean_L = (n_frames_ > 0) ? (mean_axis_length_ / static_cast<double>(n_frames_)) : 0.0;
      for (std::size_t b = 0; b < opt_.n_bins; ++b) {
        const double lo_frac = static_cast<double>(b) / static_cast<double>(opt_.n_bins);
        const double hi_frac = static_cast<double>(b + 1) / static_cast<double>(opt_.n_bins);
        const double center_frac = 0.5 * (lo_frac + hi_frac);
        const double mean_weight = (n_frames_ > 0) ? (weight_sum_[b] / static_cast<double>(n_frames_)) : 0.0;
        const double mean_density = (n_frames_ > 0) ? (density_sum_[b] / static_cast<double>(n_frames_)) : 0.0;
        ofs << b << ' '
            << std::setprecision(17) << lo_frac << ' '
            << std::setprecision(17) << hi_frac << ' '
            << std::setprecision(17) << center_frac << ' '
            << std::setprecision(17) << center_frac * mean_L << ' '
            << std::setprecision(17) << mean_weight << ' '
            << std::setprecision(17) << mean_density << ' '
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
  std::string sel_name_owned_;
  SelectionView sel_;
  Options opt_;
  std::vector<double> mass_by_atom_;
  std::vector<double> weight_sum_;
  std::vector<double> density_sum_;
  double mean_axis_length_ = 0.0;
  std::size_t n_frames_ = 0;
  bool started_ = false;
};

class CoordinationMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    double r_cut = 0.0;
  };

  struct Row {
    std::size_t frame = 0;
    std::int64_t timestep = 0;
    double time = 0.0;
    double mean_cn = 0.0;
    double std_cn = 0.0;
    double bound_fraction = 0.0;
    std::size_t total_contacts = 0;
  };

  CoordinationMeasure(std::string instance_name,
                      std::string output_path,
                      SelectionView sel_a,
                      SelectionView sel_b,
                      double dt,
                      Options opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        dt_(dt),
        opt_(opt) {
    sel_a_name_owned_ = std::string(sel_a.name);
    sel_b_name_owned_ = std::string(sel_b.name);
    sel_a_ = SelectionView{sel_a_name_owned_, sel_a.idx};
    sel_b_ = SelectionView{sel_b_name_owned_, sel_b.idx};
    same_sel_ = same_index_set(sel_a_, sel_b_);

    if (!(dt_ > 0.0)) throw std::runtime_error("coordination: dt must be > 0");
    if (opt_.range.frame_start < 0) throw std::runtime_error("coordination: frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error("coordination: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.r_cut > 0.0)) throw std::runtime_error("coordination: r_cut must be > 0");
    if (sel_a_.idx.empty() || sel_b_.idx.empty()) throw std::runtime_error("coordination: selection is empty");
  }

  std::string type() const override { return "coordination"; }
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
    od.x_axis = "frame";
    od.x_unit = "frames";
    od.columns = {"frame", "timestep", "time", "mean_cn", "std_cn", "bound_fraction", "total_contacts"};
    md.outputs.push_back(std::move(od));
    md.params["selection_b"] = std::string(sel_b_.name);
    md.params["same_selection"] = same_sel_ ? "1" : "0";
    md.params["r_cut"] = dstr(opt_.r_cut);
    md.params["dt"] = dstr(dt_);
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
    const double r2_cut = opt_.r_cut * opt_.r_cut;

    std::vector<std::size_t> cn(sel_a_.idx.size(), 0);
    for (std::size_t ia = 0; ia < sel_a_.idx.size(); ++ia) {
      const std::size_t i = sel_a_.idx[ia];
      for (std::size_t jb = 0; jb < sel_b_.idx.size(); ++jb) {
        const std::size_t j = sel_b_.idx[jb];
        if (i == j) continue;
        const auto d = frame.box.min_image_displacement(xu[i], yu[i], zu[i], xu[j], yu[j], zu[j]);
        const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if (r2 <= r2_cut) ++cn[ia];
      }
    }

    double sum = 0.0, sum2 = 0.0;
    std::size_t bound = 0;
    std::size_t total = 0;
    for (const std::size_t x : cn) {
      const double xd = static_cast<double>(x);
      sum += xd;
      sum2 += xd * xd;
      total += x;
      if (x > 0) ++bound;
    }
    const double invn = 1.0 / static_cast<double>(cn.size());
    const double mean = sum * invn;
    const double var = std::max(0.0, sum2 * invn - mean * mean);
    rows_.push_back(Row{frame_index, frame.timestep, static_cast<double>(frame.timestep) * dt_, mean, std::sqrt(var), static_cast<double>(bound) * invn, total});
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: coordination statistics within a fixed cutoff\n";
      ofs << "# selection_a: " << sel_a_.name << " (n=" << sel_a_.idx.size() << ")\n";
      ofs << "# selection_b: " << sel_b_.name << " (n=" << sel_b_.idx.size() << ")\n";
      ofs << "# same_selection: " << (same_sel_ ? 1 : 0) << "\n";
      ofs << "# r_cut: " << std::setprecision(17) << opt_.r_cut << "\n";
      ofs << "# dt: " << std::setprecision(17) << dt_ << " (seconds per timestep)\n";
      ofs << "# frame_range: [" << opt_.range.frame_start << ", " << opt_.range.frame_end << "]\n";
      ofs << "# columns: frame  timestep  time  mean_cn  std_cn  bound_fraction  total_contacts\n";
      for (const auto& r : rows_) {
        ofs << r.frame << ' ' << r.timestep << ' '
            << std::setprecision(17) << r.time << ' '
            << std::setprecision(17) << r.mean_cn << ' '
            << std::setprecision(17) << r.std_cn << ' '
            << std::setprecision(17) << r.bound_fraction << ' '
            << r.total_contacts << '\n';
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
  double dt_ = 0.0;
  Options opt_;
  bool same_sel_ = false;
  std::vector<Row> rows_;
  bool started_ = false;
};

class ResidenceACFMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    LagOptions lag;
    double r_cut = 0.0;
  };

  ResidenceACFMeasure(std::string instance_name,
                      std::string output_path,
                      SelectionView sel_a,
                      SelectionView sel_b,
                      double dt,
                      Options opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        dt_(dt),
        opt_(opt) {
    sel_a_name_owned_ = std::string(sel_a.name);
    sel_b_name_owned_ = std::string(sel_b.name);
    sel_a_ = SelectionView{sel_a_name_owned_, sel_a.idx};
    sel_b_ = SelectionView{sel_b_name_owned_, sel_b.idx};
    same_sel_ = same_index_set(sel_a_, sel_b_);

    if (!(dt_ > 0.0)) throw std::runtime_error("residence_acf: dt must be > 0");
    if (opt_.range.frame_start < 0) throw std::runtime_error("residence_acf: frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error("residence_acf: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.r_cut > 0.0)) throw std::runtime_error("residence_acf: r_cut must be > 0");
    if (sel_a_.idx.empty() || sel_b_.idx.empty()) throw std::runtime_error("residence_acf: selection is empty");
  }

  std::string type() const override { return "residence_acf"; }
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
    od.x_axis = "lag";
    od.x_unit = "frames";
    od.columns = {"lag", "time", "c_res", "sem", "count_origins", "mean_dtimestep", "mean_ref_contacts"};
    md.outputs.push_back(std::move(od));
    md.params["selection_b"] = std::string(sel_b_.name);
    md.params["same_selection"] = same_sel_ ? "1" : "0";
    md.params["r_cut"] = dstr(opt_.r_cut);
    md.params["dt"] = dstr(dt_);
    md.params["lag_stride"] = std::to_string(opt_.lag.lag_stride);
    md.params["max_lag_frames"] = std::to_string(opt_.lag.max_lag_frames);
    md.params["origin_stride"] = std::to_string(opt_.lag.origin_stride);
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
    const double r2_cut = opt_.r_cut * opt_.r_cut;

    std::vector<unsigned char> h(sel_a_.idx.size() * sel_b_.idx.size(), 0);
    for (std::size_t ia = 0; ia < sel_a_.idx.size(); ++ia) {
      const std::size_t i = sel_a_.idx[ia];
      for (std::size_t jb = 0; jb < sel_b_.idx.size(); ++jb) {
        const std::size_t j = sel_b_.idx[jb];
        if (i == j) continue;
        const auto d = frame.box.min_image_displacement(xu[i], yu[i], zu[i], xu[j], yu[j], zu[j]);
        const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if (r2 <= r2_cut) h[ia * sel_b_.idx.size() + jb] = 1;
      }
    }
    contacts_.push_back(std::move(h));
    timesteps_.push_back(frame.timestep);
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: intermittent residence/contact autocorrelation function\n";
      ofs << "# selection_a: " << sel_a_.name << " (n=" << sel_a_.idx.size() << ")\n";
      ofs << "# selection_b: " << sel_b_.name << " (n=" << sel_b_.idx.size() << ")\n";
      ofs << "# same_selection: " << (same_sel_ ? 1 : 0) << "\n";
      ofs << "# r_cut: " << std::setprecision(17) << opt_.r_cut << "\n";
      ofs << "# dt: " << std::setprecision(17) << dt_ << " (seconds per timestep)\n";
      ofs << "# lag_stride: " << opt_.lag.lag_stride << ", max_lag_frames: " << opt_.lag.max_lag_frames
          << ", origin_stride: " << opt_.lag.origin_stride << "\n";
      ofs << "# frame_range: [" << opt_.range.frame_start << ", " << opt_.range.frame_end << "]\n";
      ofs << "# columns: lag  time  c_res  sem  count_origins  mean_dtimestep  mean_ref_contacts\n";
      const auto lags = build_lag_list(contacts_.size(), opt_.lag, true);
      for (const std::size_t lag : lags) {
        double sum = 0.0, sumsq = 0.0;
        std::size_t n_valid = 0;
        double mean_ref_contacts = 0.0;
        for (std::size_t o = 0; o + lag < contacts_.size(); o += static_cast<std::size_t>(opt_.lag.origin_stride)) {
          const auto& h0 = contacts_[o];
          const auto& ht = contacts_[o + lag];
          std::size_t ref = 0;
          std::size_t overlap = 0;
          for (std::size_t k = 0; k < h0.size(); ++k) {
            if (h0[k] != 0) {
              ++ref;
              if (ht[k] != 0) ++overlap;
            }
          }
          if (ref == 0) continue;
          const double c = static_cast<double>(overlap) / static_cast<double>(ref);
          sum += c;
          sumsq += c * c;
          mean_ref_contacts += static_cast<double>(ref);
          ++n_valid;
        }
        const double mean = (n_valid > 0) ? (sum / static_cast<double>(n_valid)) : 0.0;
        const double sem = sem_from_sums(sum, sumsq, n_valid);
        const double mean_ref = (n_valid > 0) ? (mean_ref_contacts / static_cast<double>(n_valid)) : 0.0;
        const double mean_dts = mean_dtimestep_for_lag(timesteps_, lag, static_cast<std::size_t>(opt_.lag.origin_stride));
        ofs << lag << ' '
            << std::setprecision(17) << (mean_dts * dt_) << ' '
            << std::setprecision(17) << mean << ' '
            << std::setprecision(17) << sem << ' '
            << n_valid << ' '
            << std::setprecision(17) << mean_dts << ' '
            << std::setprecision(17) << mean_ref << '\n';
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
  double dt_ = 0.0;
  Options opt_;
  bool same_sel_ = false;
  std::vector<std::vector<unsigned char>> contacts_;
  std::vector<std::int64_t> timesteps_;
  bool started_ = false;
};

class BoxStatsMeasure final : public IMeasure {
public:
  struct Row {
    std::size_t frame = 0;
    std::int64_t timestep = 0;
    double time = 0.0;
    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
    double volume = 0.0;
    double rel_lx = 0.0;
    double rel_ly = 0.0;
    double rel_lz = 0.0;
    double rel_volume = 0.0;
  };

  BoxStatsMeasure(std::string instance_name,
                  std::string output_path,
                  double dt,
                  RangeOptions range)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        dt_(dt),
        range_(range) {
    if (!(dt_ > 0.0)) throw std::runtime_error("box_stats: dt must be > 0");
    if (range_.frame_start < 0) throw std::runtime_error("box_stats: frame_start must be >= 0");
    if (range_.frame_end >= 0 && range_.frame_end < range_.frame_start) {
      throw std::runtime_error("box_stats: frame_end must be -1 or >= frame_start");
    }
  }

  std::string type() const override { return "box_stats"; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = "all";
    md.n_selected = 0;

    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.x_axis = "frame";
    od.x_unit = "frames";
    od.columns = {"frame", "timestep", "time", "lx", "ly", "lz", "volume", "rel_lx", "rel_ly", "rel_lz", "rel_volume"};
    md.outputs.push_back(std::move(od));
    md.params["dt"] = dstr(dt_);
    md.params["frame_start"] = std::to_string(range_.frame_start);
    md.params["frame_end"] = std::to_string(range_.frame_end);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    ref_lx_ = first_frame.box.lx();
    ref_ly_ = first_frame.box.ly();
    ref_lz_ = first_frame.box.lz();
    ref_volume_ = box_volume(first_frame.box);
    started_ = true;
    if (!range_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, range_)) return;
    const double lx = frame.box.lx();
    const double ly = frame.box.ly();
    const double lz = frame.box.lz();
    const double vol = box_volume(frame.box);
    rows_.push_back(Row{frame_index, frame.timestep, static_cast<double>(frame.timestep) * dt_, lx, ly, lz, vol, lx / ref_lx_, ly / ref_ly_, lz / ref_lz_, vol / ref_volume_});
  }

  void flush_partial() override {
    if (range_.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: box-length and swelling statistics\n";
      ofs << "# dt: " << std::setprecision(17) << dt_ << " (seconds per timestep)\n";
      ofs << "# frame_range: [" << range_.frame_start << ", " << range_.frame_end << "]\n";
      ofs << "# reference_lengths: " << std::setprecision(17) << ref_lx_ << ' ' << ref_ly_ << ' ' << ref_lz_ << "\n";
      ofs << "# reference_volume: " << std::setprecision(17) << ref_volume_ << "\n";
      ofs << "# columns: frame  timestep  time  lx  ly  lz  volume  rel_lx  rel_ly  rel_lz  rel_volume\n";
      for (const auto& r : rows_) {
        ofs << r.frame << ' ' << r.timestep << ' '
            << std::setprecision(17) << r.time << ' '
            << std::setprecision(17) << r.lx << ' '
            << std::setprecision(17) << r.ly << ' '
            << std::setprecision(17) << r.lz << ' '
            << std::setprecision(17) << r.volume << ' '
            << std::setprecision(17) << r.rel_lx << ' '
            << std::setprecision(17) << r.rel_ly << ' '
            << std::setprecision(17) << r.rel_lz << ' '
            << std::setprecision(17) << r.rel_volume << '\n';
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
  double dt_ = 0.0;
  RangeOptions range_;
  double ref_lx_ = 1.0;
  double ref_ly_ = 1.0;
  double ref_lz_ = 1.0;
  double ref_volume_ = 1.0;
  std::vector<Row> rows_;
  bool started_ = false;
};

void append_static_selection_caps(const IniConfig& cfg,
                                  const std::string& section,
                                  MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  if (cfg.has_key(section, "group_b")) caps.group_refs.push_back(cfg.get_string(section, "group_b"));
}

MeasureCapabilities profile_caps(const IniConfig& cfg,
                                 const std::string& section,
                                 const std::string& instance,
                                 const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  append_static_selection_caps(cfg, section, caps);
  const ProfileMode mode = parse_profile_mode(cfg.get_string(section, "mode", std::optional<std::string>("number")));
  if (mode == ProfileMode::Charge) {
    caps.requires_dfields.push_back(cfg.get_string(section, "charge_field", std::optional<std::string>("q")));
  } else if (mode == ProfileMode::Mass) {
    if (cfg.has_key(section, "mass_field")) {
      caps.requires_dfields.push_back(cfg.get_string(section, "mass_field"));
    } else {
      caps.requires_intfields.push_back("type");
      caps.requires_topology_sections.push_back("masses");
    }
  }
  return caps;
}

std::unique_ptr<IMeasure> profile_create(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env,
                                         const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("profile1d factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "profile1d");

  const ProfileMode mode = parse_profile_mode(cfg.get_string(section, "mode", std::optional<std::string>("number")));
  std::vector<double> mass_by_atom(frame0.natoms, 1.0);
  if (mode == ProfileMode::Mass) mass_by_atom = mass_by_atom_from_config(cfg, section, frame0, sysctx);

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("profile1d.dat"))).lexically_normal();

  Profile1DMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.axis = parse_axis1d(cfg.get_string(section, "axis", std::optional<std::string>("z")));
  opt.n_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_bins", std::optional<std::int64_t>(100)));
  opt.mode = mode;
  opt.charge_field = cfg.get_string(section, "charge_field", std::optional<std::string>("q"));

  return std::make_unique<Profile1DMeasure>(instance, out.string(), sel, opt, std::move(mass_by_atom));
}

MeasureCapabilities coordination_caps(const IniConfig& cfg,
                                      const std::string& section,
                                      const std::string& instance,
                                      const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_static_selection_caps(cfg, section, caps);
  return caps;
}

std::unique_ptr<IMeasure> coordination_create(const IniConfig& cfg,
                                              const std::string& section,
                                              const std::string& instance,
                                              const MeasureBuildEnv& env,
                                              const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("coordination factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group_a = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo_a = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb_a = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  const std::string group_b = cfg.get_string(section, "group_b", std::optional<std::string>(group_a));
  const std::string topo_b = cfg.get_string(section, "topo_group_b", std::optional<std::string>(topo_a));
  const std::string comb_b = cfg.get_string(section, "combine_b", std::optional<std::string>(comb_a));
  SelectionView sel_a = get_static_combined_view(*env.selection_provider, frame0, group_a, topo_a, comb_a, "coordination");
  SelectionView sel_b = get_static_combined_view(*env.selection_provider, frame0, group_b, topo_b, comb_b, "coordination");

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("coordination.dat"))).lexically_normal();

  CoordinationMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.r_cut = cfg.get_double(section, "r_cut");
  return std::make_unique<CoordinationMeasure>(instance, out.string(), sel_a, sel_b, env.dt, opt);
}

MeasureCapabilities residence_caps(const IniConfig& cfg,
                                   const std::string& section,
                                   const std::string& instance,
                                   const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_static_selection_caps(cfg, section, caps);
  caps.requires_identity_consistent = true;
  return caps;
}

std::unique_ptr<IMeasure> residence_create(const IniConfig& cfg,
                                           const std::string& section,
                                           const std::string& instance,
                                           const MeasureBuildEnv& env,
                                           const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("residence_acf factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group_a = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo_a = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb_a = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  const std::string group_b = cfg.get_string(section, "group_b", std::optional<std::string>(group_a));
  const std::string topo_b = cfg.get_string(section, "topo_group_b", std::optional<std::string>(topo_a));
  const std::string comb_b = cfg.get_string(section, "combine_b", std::optional<std::string>(comb_a));
  SelectionView sel_a = get_static_combined_view(*env.selection_provider, frame0, group_a, topo_a, comb_a, "residence_acf");
  SelectionView sel_b = get_static_combined_view(*env.selection_provider, frame0, group_b, topo_b, comb_b, "residence_acf");

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("residence_acf.dat"))).lexically_normal();

  ResidenceACFMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.lag = lag_options_from_config(cfg, section);
  opt.r_cut = cfg.get_double(section, "r_cut");
  return std::make_unique<ResidenceACFMeasure>(instance, out.string(), sel_a, sel_b, env.dt, opt);
}

MeasureCapabilities box_stats_caps(const IniConfig& cfg,
                                   const std::string& section,
                                   const std::string& instance,
                                   const MeasureBuildEnv& env) {
  (void)cfg;
  (void)section;
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  caps.scale = ScaleCompatibility{true, true, true};
  return caps;
}

std::unique_ptr<IMeasure> box_stats_create(const IniConfig& cfg,
                                           const std::string& section,
                                           const std::string& instance,
                                           const MeasureBuildEnv& env,
                                           const SystemContext& sysctx) {
  (void)sysctx;
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("box_stats.dat"))).lexically_normal();

  RangeOptions range;
  range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  range.dry_run = env.dry_run;
  return std::make_unique<BoxStatsMeasure>(instance, out.string(), env.dt, range);
}

static MeasureRegistrar g_register_profile("profile1d", &profile_caps, &profile_create);
static MeasureRegistrar g_register_coordination("coordination", &coordination_caps, &coordination_create);
static MeasureRegistrar g_register_residence("residence_acf", &residence_caps, &residence_create);
static MeasureRegistrar g_register_box_stats("box_stats", &box_stats_caps, &box_stats_create);
static MeasureRegistrar g_register_swelling_stats_alias("swelling_stats", &box_stats_caps, &box_stats_create);

} // namespace
} // namespace pilots
