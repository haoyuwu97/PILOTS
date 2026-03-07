#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <optional>
#include <span>
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

using measure_ext::count_diag_dims;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::get_static_group_view;
using measure_ext::parse_diag_mask;
using measure_ext::parse_double_list;
using measure_ext::parse_int_list;
using measure_ext::resolve_finite_frame_end;
using measure_ext::resolve_path;

struct RangeOptions {
  std::int64_t frame_start = 0;
  std::int64_t frame_end = -1;
  bool dry_run = false;
};

inline bool frame_in_range(std::size_t frame_index, const RangeOptions& opt) {
  const std::int64_t fi = static_cast<std::int64_t>(frame_index);
  if (fi < opt.frame_start) return false;
  if (opt.frame_end >= 0 && fi > opt.frame_end) return false;
  return true;
}

enum class TwoTimeMode {
  MSD,
  FS,
  Q,
  CHI4,
};

inline TwoTimeMode parse_two_time_mode(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  if (s == "msd" || s == "msd_tw") return TwoTimeMode::MSD;
  if (s == "fs" || s == "self_isf" || s == "fs_tw") return TwoTimeMode::FS;
  if (s == "q" || s == "overlap" || s == "q_tw") return TwoTimeMode::Q;
  if (s == "chi4" || s == "chi4_tw" || s == "chi4_overlap") return TwoTimeMode::CHI4;
  throw std::runtime_error("invalid two-time observable: '" + s + "' (use msd|fs|q|chi4)");
}

inline const char* two_time_mode_name(TwoTimeMode m) {
  switch (m) {
    case TwoTimeMode::MSD: return "msd";
    case TwoTimeMode::FS: return "fs";
    case TwoTimeMode::Q: return "q";
    case TwoTimeMode::CHI4: return "chi4";
  }
  return "msd";
}

inline double spherical_bessel_j0(double x) {
  if (std::abs(x) < 1e-14) return 1.0;
  return std::sin(x) / x;
}

inline std::vector<std::size_t> parse_nonnegative_index_list(const IniConfig& cfg,
                                                             const std::string& section,
                                                             const std::string& key) {
  const auto vals = parse_int_list(cfg, section, key);
  std::vector<std::size_t> out;
  out.reserve(vals.size());
  for (const int v : vals) {
    if (v < 0) {
      throw std::runtime_error(section + ": " + key + " must contain non-negative integers only");
    }
    out.push_back(static_cast<std::size_t>(v));
  }
  return out;
}

inline std::vector<std::size_t> build_waiting_anchor_list(const IniConfig& cfg,
                                                          const std::string& section,
                                                          std::int64_t frame_start,
                                                          std::int64_t frame_end) {
  std::vector<std::size_t> out;
  if (cfg.has_key(section, "waiting_frames")) {
    out = parse_nonnegative_index_list(cfg, section, "waiting_frames");
  } else {
    const std::int64_t wstart = cfg.get_int64(section, "waiting_start_frame", std::optional<std::int64_t>(frame_start));
    const std::int64_t wend = cfg.get_int64(section, "waiting_end_frame", std::optional<std::int64_t>(frame_end));
    const std::int64_t wstride = cfg.get_int64(section, "waiting_stride_frames", std::optional<std::int64_t>(1));
    if (wstart < frame_start) throw std::runtime_error(section + ": waiting_start_frame must be >= frame_start");
    if (wend > frame_end) throw std::runtime_error(section + ": waiting_end_frame must be <= frame_end");
    if (wend < wstart) throw std::runtime_error(section + ": waiting_end_frame must be >= waiting_start_frame");
    if (wstride <= 0) throw std::runtime_error(section + ": waiting_stride_frames must be >= 1");
    for (std::int64_t w = wstart; w <= wend; w += wstride) {
      out.push_back(static_cast<std::size_t>(w));
    }
  }
  if (out.empty()) throw std::runtime_error(section + ": waiting anchor list is empty");
  for (const std::size_t w : out) {
    if (static_cast<std::int64_t>(w) < frame_start || static_cast<std::int64_t>(w) > frame_end) {
      throw std::runtime_error(section + ": a waiting frame lies outside [frame_start, frame_end]");
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

inline std::vector<std::size_t> build_lag_list(const IniConfig& cfg,
                                               const std::string& section,
                                               std::int64_t frame_start,
                                               std::int64_t frame_end) {
  std::vector<std::size_t> out;
  if (cfg.has_key(section, "lag_frames")) {
    out = parse_nonnegative_index_list(cfg, section, "lag_frames");
  } else {
    const std::int64_t max_lag = cfg.get_int64(section, "max_lag_frames", std::optional<std::int64_t>(frame_end - frame_start));
    const std::int64_t lag_stride = cfg.get_int64(section, "lag_stride_frames", std::optional<std::int64_t>(1));
    if (max_lag < 0) throw std::runtime_error(section + ": max_lag_frames must be >= 0");
    if (lag_stride <= 0) throw std::runtime_error(section + ": lag_stride_frames must be >= 1");
    const std::size_t max_lag_eff = static_cast<std::size_t>(std::min<std::int64_t>(max_lag, frame_end - frame_start));
    for (std::size_t lag = 0; lag <= max_lag_eff; lag += static_cast<std::size_t>(lag_stride)) {
      out.push_back(lag);
    }
  }
  if (out.empty()) throw std::runtime_error(section + ": lag list is empty");
  const std::size_t max_allowed = static_cast<std::size_t>(frame_end - frame_start);
  for (const std::size_t lag : out) {
    if (lag > max_allowed) {
      throw std::runtime_error(section + ": lag_frames contains a lag larger than the available analysis window");
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

inline void append_static_caps(const IniConfig& cfg,
                               const std::string& section,
                               MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  if (remove_drift) {
    caps.group_refs.push_back(cfg.get_string(section, "drift_group", std::optional<std::string>("all")));
  }
}

struct StoredFrame {
  std::size_t frame_index_abs = 0;
  std::int64_t timestep = 0;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
};

struct TwoTimeOptions {
  RangeOptions range;
  std::vector<std::size_t> waiting_frames_abs;
  std::size_t waiting_window_frames = 1;
  std::size_t origin_stride_frames = 1;
  std::vector<std::size_t> lag_frames;
  int diag_mask = 7;
  bool remove_drift = true;
  double dt = 0.0;
  std::vector<double> q_values;
  std::vector<double> a_values;
  std::string observable_name;
};

class TwoTimeMeasure final : public IMeasure {
public:
  TwoTimeMeasure(std::string type_name,
                 std::string instance_name,
                 std::string output_path,
                 SelectionView sel,
                 SelectionView drift_sel,
                 TwoTimeMode mode,
                 TwoTimeOptions opt)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        mode_(mode),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error(type_name_ + ": selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) {
      throw std::runtime_error(type_name_ + ": remove_drift requested but drift selection is empty");
    }
    if (opt_.range.frame_start < 0) throw std::runtime_error(type_name_ + ": frame_start must be >= 0");
    if (opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error(type_name_ + ": frame_end must be >= frame_start");
    }
    if (opt_.waiting_window_frames == 0) {
      throw std::runtime_error(type_name_ + ": waiting_window_frames must be >= 1");
    }
    if (opt_.origin_stride_frames == 0) {
      throw std::runtime_error(type_name_ + ": origin_stride_frames must be >= 1");
    }
    if (mode_ == TwoTimeMode::CHI4 && opt_.waiting_window_frames < 2) {
      throw std::runtime_error(type_name_ + ": chi4_tw requires waiting_window_frames >= 2 to estimate fluctuations across time origins");
    }
    if (mode_ == TwoTimeMode::FS && opt_.q_values.empty()) {
      throw std::runtime_error(type_name_ + ": q_values must not be empty");
    }
    if ((mode_ == TwoTimeMode::Q || mode_ == TwoTimeMode::CHI4) && opt_.a_values.empty()) {
      throw std::runtime_error(type_name_ + ": a_values must not be empty");
    }
  }

  std::string type() const override { return type_name_; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type_name_;
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();

    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.x_axis = "waiting_lag_grid";
    od.x_unit = "mixed";
    od.columns = common_columns_();
    md.outputs.push_back(std::move(od));

    md.params["observable"] = opt_.observable_name;
    md.params["frame_start"] = std::to_string(opt_.range.frame_start);
    md.params["frame_end"] = std::to_string(opt_.range.frame_end);
    md.params["waiting_window_frames"] = std::to_string(opt_.waiting_window_frames);
    md.params["origin_stride_frames"] = std::to_string(opt_.origin_stride_frames);
    md.params["waiting_frames"] = join_size_t_(opt_.waiting_frames_abs);
    md.params["lag_frames"] = join_size_t_(opt_.lag_frames);
    md.params["components"] = components_name_();
    md.params["remove_drift"] = opt_.remove_drift ? "true" : "false";
    if (!opt_.q_values.empty()) md.params["q_values"] = measure_ext::join_doubles(opt_.q_values);
    if (!opt_.a_values.empty()) md.params["a_values"] = measure_ext::join_doubles(opt_.a_values);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    double dx = 0.0, dy = 0.0, dz = 0.0;
    if (opt_.remove_drift) {
      for (const std::size_t i : drift_sel_.idx) {
        dx += xu[i];
        dy += yu[i];
        dz += zu[i];
      }
      const double inv = 1.0 / static_cast<double>(drift_sel_.idx.size());
      dx *= inv; dy *= inv; dz *= inv;
    }

    StoredFrame sf;
    sf.frame_index_abs = frame_index;
    sf.timestep = frame.timestep;
    sf.x.resize(sel_.idx.size());
    sf.y.resize(sel_.idx.size());
    sf.z.resize(sel_.idx.size());
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t i = sel_.idx[p];
      sf.x[p] = xu[i] - dx;
      sf.y[p] = yu[i] - dy;
      sf.z[p] = zu[i] - dz;
    }
    frames_.push_back(std::move(sf));
  }

  void flush_partial() override {}

  void finalize() override {
    if (!started_ || opt_.range.dry_run) return;
    if (frames_.empty()) {
      throw std::runtime_error(type_name_ + ": no frames were buffered in the requested analysis window");
    }
    write_output_();
  }

private:
  double displacement_sq_(const StoredFrame& cur, const StoredFrame& org, std::size_t p) const {
    double dr2 = 0.0;
    if ((opt_.diag_mask & 1) != 0) {
      const double dx = cur.x[p] - org.x[p];
      dr2 += dx * dx;
    }
    if ((opt_.diag_mask & 2) != 0) {
      const double dy = cur.y[p] - org.y[p];
      dr2 += dy * dy;
    }
    if ((opt_.diag_mask & 4) != 0) {
      const double dz = cur.z[p] - org.z[p];
      dr2 += dz * dz;
    }
    return dr2;
  }

  static double sem_from_sum_(double sum, double sumsq, std::size_t n) {
    if (n <= 1) return 0.0;
    const double mean = sum / static_cast<double>(n);
    const double var = std::max(0.0, (sumsq - static_cast<double>(n) * mean * mean) / static_cast<double>(n - 1));
    return std::sqrt(var / static_cast<double>(n));
  }

  std::vector<std::string> common_columns_() const {
    switch (mode_) {
      case TwoTimeMode::MSD:
        return {"waiting_frame", "waiting_timestep", "waiting_time", "lag_frames", "lag_dtimestep", "lag_time", "msd", "sem", "n_origins"};
      case TwoTimeMode::FS:
        return {"waiting_frame", "waiting_timestep", "waiting_time", "lag_frames", "lag_dtimestep", "lag_time", "q", "fs", "sem", "n_origins"};
      case TwoTimeMode::Q:
        return {"waiting_frame", "waiting_timestep", "waiting_time", "lag_frames", "lag_dtimestep", "lag_time", "a", "q", "sem", "n_origins"};
      case TwoTimeMode::CHI4:
        return {"waiting_frame", "waiting_timestep", "waiting_time", "lag_frames", "lag_dtimestep", "lag_time", "a", "q_mean", "q2_mean", "chi4", "q_sem", "n_origins"};
    }
    return {};
  }

  std::string components_name_() const {
    switch (opt_.diag_mask) {
      case 1: return "xx";
      case 2: return "yy";
      case 4: return "zz";
      case 3: return "xxyy";
      case 5: return "xxzz";
      case 6: return "yyzz";
      case 7: return "xxyyzz";
      default: return "custom";
    }
  }

  static std::string join_size_t_(const std::vector<std::size_t>& vals) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < vals.size(); ++i) {
      if (i) oss << ",";
      oss << vals[i];
    }
    return oss.str();
  }

  void write_output_() const {
    const std::int64_t t0 = frames_.front().timestep;
    const std::size_t nframes = frames_.size();
    const std::size_t nsel = sel_.idx.size();

    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: two-time / aging observable\n";
      ofs << "# type: " << type_name_ << ", observable: " << opt_.observable_name << "\n";
      ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
      ofs << "# frame_range: [" << opt_.range.frame_start << ", " << opt_.range.frame_end << "]\n";
      ofs << "# waiting_window_frames: " << opt_.waiting_window_frames << ", origin_stride_frames: " << opt_.origin_stride_frames << "\n";
      ofs << "# waiting_frames: " << join_size_t_(opt_.waiting_frames_abs) << "\n";
      ofs << "# lag_frames: " << join_size_t_(opt_.lag_frames) << "\n";
      ofs << "# components: " << components_name_() << " (n_dims=" << count_diag_dims(opt_.diag_mask) << ")\n";
      ofs << "# remove_drift: " << (opt_.remove_drift ? "true" : "false") << "\n";
      if (!opt_.q_values.empty()) ofs << "# q_values: " << dstr(opt_.q_values.front()) << ((opt_.q_values.size() > 1) ? (",..." ) : "") << "\n";
      if (!opt_.a_values.empty()) ofs << "# a_values: " << dstr(opt_.a_values.front()) << ((opt_.a_values.size() > 1) ? (",..." ) : "") << "\n";
      ofs << "# columns:";
      for (const auto& c : common_columns_()) ofs << ' ' << c;
      ofs << "\n";

      for (const std::size_t waiting_abs : opt_.waiting_frames_abs) {
        const std::size_t waiting_rel = waiting_abs - static_cast<std::size_t>(opt_.range.frame_start);
        if (waiting_rel >= nframes) continue;
        const std::int64_t waiting_ts = frames_[waiting_rel].timestep;
        const double waiting_time = static_cast<double>(waiting_ts - t0) * opt_.dt;

        for (const std::size_t lag : opt_.lag_frames) {
          if (waiting_rel + lag >= nframes) continue;

          if (mode_ == TwoTimeMode::MSD) {
            double sum = 0.0, sumsq = 0.0, sum_dt = 0.0;
            std::size_t norig = 0;
            for (std::size_t o = waiting_rel; o < std::min(nframes, waiting_rel + opt_.waiting_window_frames); o += opt_.origin_stride_frames) {
              if (o + lag >= nframes) break;
              double acc = 0.0;
              for (std::size_t p = 0; p < nsel; ++p) acc += displacement_sq_(frames_[o + lag], frames_[o], p);
              const double v = acc / static_cast<double>(nsel);
              sum += v;
              sumsq += v * v;
              sum_dt += static_cast<double>(frames_[o + lag].timestep - frames_[o].timestep);
              ++norig;
            }
            if (norig == 0) continue;
            const double mean = sum / static_cast<double>(norig);
            const double sem = sem_from_sum_(sum, sumsq, norig);
            const double mean_dt = sum_dt / static_cast<double>(norig);
            ofs << waiting_abs << ' ' << waiting_ts << ' ' << std::setprecision(17) << waiting_time << ' '
                << lag << ' ' << std::setprecision(17) << mean_dt << ' ' << std::setprecision(17) << (mean_dt * opt_.dt) << ' '
                << std::setprecision(17) << mean << ' ' << std::setprecision(17) << sem << ' ' << norig << '\n';
            continue;
          }

          if (mode_ == TwoTimeMode::FS) {
            for (const double q : opt_.q_values) {
              double sum = 0.0, sumsq = 0.0, sum_dt = 0.0;
              std::size_t norig = 0;
              for (std::size_t o = waiting_rel; o < std::min(nframes, waiting_rel + opt_.waiting_window_frames); o += opt_.origin_stride_frames) {
                if (o + lag >= nframes) break;
                double acc = 0.0;
                for (std::size_t p = 0; p < nsel; ++p) {
                  acc += spherical_bessel_j0(q * std::sqrt(displacement_sq_(frames_[o + lag], frames_[o], p)));
                }
                const double v = acc / static_cast<double>(nsel);
                sum += v;
                sumsq += v * v;
                sum_dt += static_cast<double>(frames_[o + lag].timestep - frames_[o].timestep);
                ++norig;
              }
              if (norig == 0) continue;
              const double mean = sum / static_cast<double>(norig);
              const double sem = sem_from_sum_(sum, sumsq, norig);
              const double mean_dt = sum_dt / static_cast<double>(norig);
              ofs << waiting_abs << ' ' << waiting_ts << ' ' << std::setprecision(17) << waiting_time << ' '
                  << lag << ' ' << std::setprecision(17) << mean_dt << ' ' << std::setprecision(17) << (mean_dt * opt_.dt) << ' '
                  << std::setprecision(17) << q << ' ' << std::setprecision(17) << mean << ' '
                  << std::setprecision(17) << sem << ' ' << norig << '\n';
            }
            continue;
          }

          if (mode_ == TwoTimeMode::Q || mode_ == TwoTimeMode::CHI4) {
            for (const double a : opt_.a_values) {
              const double a2 = a * a;
              double sum = 0.0, sumsq = 0.0, sum_dt = 0.0;
              std::size_t norig = 0;
              for (std::size_t o = waiting_rel; o < std::min(nframes, waiting_rel + opt_.waiting_window_frames); o += opt_.origin_stride_frames) {
                if (o + lag >= nframes) break;
                double acc = 0.0;
                for (std::size_t p = 0; p < nsel; ++p) {
                  acc += (displacement_sq_(frames_[o + lag], frames_[o], p) <= a2) ? 1.0 : 0.0;
                }
                const double qv = acc / static_cast<double>(nsel);
                sum += qv;
                sumsq += qv * qv;
                sum_dt += static_cast<double>(frames_[o + lag].timestep - frames_[o].timestep);
                ++norig;
              }
              if (norig == 0) continue;
              const double qmean = sum / static_cast<double>(norig);
              const double mean_dt = sum_dt / static_cast<double>(norig);
              if (mode_ == TwoTimeMode::Q) {
                const double sem = sem_from_sum_(sum, sumsq, norig);
                ofs << waiting_abs << ' ' << waiting_ts << ' ' << std::setprecision(17) << waiting_time << ' '
                    << lag << ' ' << std::setprecision(17) << mean_dt << ' ' << std::setprecision(17) << (mean_dt * opt_.dt) << ' '
                    << std::setprecision(17) << a << ' ' << std::setprecision(17) << qmean << ' '
                    << std::setprecision(17) << sem << ' ' << norig << '\n';
              } else {
                if (norig < 2) continue;
                const double q2mean = sumsq / static_cast<double>(norig);
                const double chi4 = static_cast<double>(nsel) * std::max(0.0, q2mean - qmean * qmean);
                const double qsem = sem_from_sum_(sum, sumsq, norig);
                ofs << waiting_abs << ' ' << waiting_ts << ' ' << std::setprecision(17) << waiting_time << ' '
                    << lag << ' ' << std::setprecision(17) << mean_dt << ' ' << std::setprecision(17) << (mean_dt * opt_.dt) << ' '
                    << std::setprecision(17) << a << ' ' << std::setprecision(17) << qmean << ' '
                    << std::setprecision(17) << q2mean << ' ' << std::setprecision(17) << chi4 << ' '
                    << std::setprecision(17) << qsem << ' ' << norig << '\n';
              }
            }
          }
        }
      }
    });
  }

  std::string type_name_;
  std::string instance_name_;
  std::string output_path_;
  TwoTimeMode mode_ = TwoTimeMode::MSD;
  TwoTimeOptions opt_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  bool started_ = false;
  std::vector<StoredFrame> frames_;
};

TwoTimeOptions parse_two_time_options(const IniConfig& cfg,
                                      const std::string& section,
                                      const MeasureBuildEnv& env,
                                      const std::string& type_name,
                                      TwoTimeMode mode) {
  TwoTimeOptions opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  const std::int64_t frame_end_cfg = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.frame_end = resolve_finite_frame_end(env.follow, opt.range.frame_start, frame_end_cfg, env.input_path, type_name);
  opt.range.dry_run = env.dry_run;
  opt.waiting_frames_abs = build_waiting_anchor_list(cfg, section, opt.range.frame_start, opt.range.frame_end);
  opt.waiting_window_frames = static_cast<std::size_t>(cfg.get_int64(section, "waiting_window_frames", std::optional<std::int64_t>((mode == TwoTimeMode::CHI4) ? 8 : 1)));
  opt.origin_stride_frames = static_cast<std::size_t>(cfg.get_int64(section, "origin_stride_frames", std::optional<std::int64_t>(1)));
  opt.lag_frames = build_lag_list(cfg, section, opt.range.frame_start, opt.range.frame_end);
  opt.diag_mask = parse_diag_mask(cfg.get_string(section, "components", std::optional<std::string>("xxyyzz")));
  opt.remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  opt.dt = env.dt;
  opt.observable_name = two_time_mode_name(mode);
  if (mode == TwoTimeMode::FS) opt.q_values = parse_double_list(cfg, section, "q_values");
  if (mode == TwoTimeMode::Q || mode == TwoTimeMode::CHI4) opt.a_values = parse_double_list(cfg, section, "a_values");
  return opt;
}

std::unique_ptr<IMeasure> create_two_time_measure(const IniConfig& cfg,
                                                  const std::string& section,
                                                  const std::string& instance,
                                                  const MeasureBuildEnv& env,
                                                  const SystemContext& sysctx,
                                                  std::string type_name,
                                                  TwoTimeMode mode) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) {
    throw std::runtime_error(type_name + " factory: missing first_frame or SelectionProvider");
  }
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);

  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
  SelectionView drift_sel = remove_drift
      ? get_static_group_view(*env.selection_provider, frame0, drift_group, type_name)
      : get_static_group_view(*env.selection_provider, frame0, "all", type_name);

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>(type_name + std::string(".dat")))).lexically_normal();

  TwoTimeOptions opt = parse_two_time_options(cfg, section, env, type_name, mode);
  return std::make_unique<TwoTimeMeasure>(std::move(type_name), instance, out.string(), sel, drift_sel, mode, std::move(opt));
}

MeasureCapabilities msd_tw_caps(const IniConfig& cfg,
                                const std::string& section,
                                const std::string& instance,
                                const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  append_static_caps(cfg, section, caps);
  return caps;
}

MeasureCapabilities fs_tw_caps(const IniConfig& cfg,
                               const std::string& section,
                               const std::string& instance,
                               const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  append_static_caps(cfg, section, caps);
  return caps;
}

MeasureCapabilities q_tw_caps(const IniConfig& cfg,
                              const std::string& section,
                              const std::string& instance,
                              const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  append_static_caps(cfg, section, caps);
  return caps;
}

MeasureCapabilities chi4_tw_caps(const IniConfig& cfg,
                                 const std::string& section,
                                 const std::string& instance,
                                 const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  append_static_caps(cfg, section, caps);
  return caps;
}

MeasureCapabilities waiting_time_scan_caps(const IniConfig& cfg,
                                           const std::string& section,
                                           const std::string& instance,
                                           const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  append_static_caps(cfg, section, caps);
  return caps;
}

std::unique_ptr<IMeasure> msd_tw_create(const IniConfig& cfg,
                                        const std::string& section,
                                        const std::string& instance,
                                        const MeasureBuildEnv& env,
                                        const SystemContext& sysctx) {
  return create_two_time_measure(cfg, section, instance, env, sysctx, "msd_tw", TwoTimeMode::MSD);
}

std::unique_ptr<IMeasure> fs_tw_create(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& instance,
                                       const MeasureBuildEnv& env,
                                       const SystemContext& sysctx) {
  return create_two_time_measure(cfg, section, instance, env, sysctx, "fs_tw", TwoTimeMode::FS);
}

std::unique_ptr<IMeasure> q_tw_create(const IniConfig& cfg,
                                      const std::string& section,
                                      const std::string& instance,
                                      const MeasureBuildEnv& env,
                                      const SystemContext& sysctx) {
  return create_two_time_measure(cfg, section, instance, env, sysctx, "q_tw", TwoTimeMode::Q);
}

std::unique_ptr<IMeasure> chi4_tw_create(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env,
                                         const SystemContext& sysctx) {
  return create_two_time_measure(cfg, section, instance, env, sysctx, "chi4_tw", TwoTimeMode::CHI4);
}

std::unique_ptr<IMeasure> waiting_time_scan_create(const IniConfig& cfg,
                                                   const std::string& section,
                                                   const std::string& instance,
                                                   const MeasureBuildEnv& env,
                                                   const SystemContext& sysctx) {
  const TwoTimeMode mode = parse_two_time_mode(cfg.get_string(section, "observable"));
  return create_two_time_measure(cfg, section, instance, env, sysctx, "waiting_time_scan", mode);
}

static MeasureRegistrar g_register_msd_tw("msd_tw", &msd_tw_caps, &msd_tw_create);
static MeasureRegistrar g_register_fs_tw("fs_tw", &fs_tw_caps, &fs_tw_create);
static MeasureRegistrar g_register_q_tw("q_tw", &q_tw_caps, &q_tw_create);
static MeasureRegistrar g_register_chi4_tw("chi4_tw", &chi4_tw_caps, &chi4_tw_create);
static MeasureRegistrar g_register_waiting_time_scan("waiting_time_scan", &waiting_time_scan_caps, &waiting_time_scan_create);

} // namespace
} // namespace pilots
