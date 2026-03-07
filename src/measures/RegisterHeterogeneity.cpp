#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
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

using measure_ext::Vec3;
using measure_ext::atom_vec3;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::get_static_group_view;
using measure_ext::integer_like_field_to_i64;
using measure_ext::norm;
using measure_ext::resolve_path;
using measure_ext::x_unit_for_axis;

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
  if (opt.max_lag_frames == 0 || opt.max_lag_frames < -1) {
    throw std::runtime_error(section + ": max_lag_frames must be -1 or >= 1");
  }
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

inline double j0(double x) {
  const double ax = std::abs(x);
  if (ax < 1e-12) return 1.0;
  return std::sin(x) / x;
}

inline double mic_delta(double d, double L) {
  if (!(L > 0.0)) return d;
  return d - L * std::round(d / L);
}

inline Vec3 mic_vec(const Vec3& a, const Vec3& b, double lx, double ly, double lz) {
  return Vec3{mic_delta(a.x - b.x, lx), mic_delta(a.y - b.y, ly), mic_delta(a.z - b.z, lz)};
}

struct StoredFrame {
  std::int64_t timestep = 0;
  double lx = 0.0, ly = 0.0, lz = 0.0;
  std::vector<Vec3> pos;
};

struct DSU {
  explicit DSU(std::size_t n) : p(n), r(n, 0) {
    std::iota(p.begin(), p.end(), std::size_t{0});
  }
  std::size_t find(std::size_t x) {
    while (p[x] != x) {
      p[x] = p[p[x]];
      x = p[x];
    }
    return x;
  }
  void unite(std::size_t a, std::size_t b) {
    a = find(a); b = find(b);
    if (a == b) return;
    if (r[a] < r[b]) std::swap(a, b);
    p[b] = a;
    if (r[a] == r[b]) ++r[a];
  }
  std::vector<std::size_t> p;
  std::vector<int> r;
};

class BufferedSelectedPositionsMeasure : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    bool remove_drift = false;
  };

  BufferedSelectedPositionsMeasure(std::string type_name,
                                   std::string instance_name,
                                   SelectionView sel,
                                   SelectionView drift_sel,
                                   Options opt)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error(type_name_ + ": selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) {
      throw std::runtime_error(type_name_ + ": remove_drift requested but drift selection is empty");
    }
  }

  std::string type() const override { return type_name_; }
  std::string instance_name() const override { return instance_name_; }

  void on_start(const Frame& first_frame) override { (void)first_frame; started_ = true; }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    double comx = 0.0, comy = 0.0, comz = 0.0;
    if (opt_.remove_drift) {
      for (const std::size_t i : drift_sel_.idx) {
        comx += xu[i];
        comy += yu[i];
        comz += zu[i];
      }
      const double inv = 1.0 / static_cast<double>(drift_sel_.idx.size());
      comx *= inv; comy *= inv; comz *= inv;
    }

    StoredFrame sf;
    sf.timestep = frame.timestep;
    sf.lx = frame.box.lx();
    sf.ly = frame.box.ly();
    sf.lz = frame.box.lz();
    sf.pos.resize(sel_.idx.size());
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t i = sel_.idx[p];
      sf.pos[p] = Vec3{xu[i] - comx, yu[i] - comy, zu[i] - comz};
    }
    frames_.push_back(std::move(sf));
  }

  void flush_partial() override {}
  void finalize() override { if (started_ && !opt_.range.dry_run) finalize_buffered(); }

protected:
  virtual void finalize_buffered() = 0;

  std::string type_name_;
  std::string instance_name_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  Options opt_;
  bool started_ = false;
  std::vector<StoredFrame> frames_;
};

enum class FourPointMode {
  G4_RT,
  S4_QT,
  XI4_T,
};

class FourPointMeasure final : public BufferedSelectedPositionsMeasure {
public:
  struct Options {
    BufferedSelectedPositionsMeasure::Options base;
    LagOptions lag;
    double overlap_a = 1.0;
    double r_max = 0.0;
    double dr = 0.1;
    std::vector<double> q_list;
    std::size_t fit_n_q = 3;
  };

  FourPointMeasure(std::string type_name,
                   std::string instance_name,
                   std::string output_path,
                   SelectionView sel,
                   SelectionView drift_sel,
                   FourPointMode mode,
                   Options opt)
      : BufferedSelectedPositionsMeasure(type_name, instance_name, sel, drift_sel, opt.base),
        output_path_(std::move(output_path)),
        mode_(mode),
        opt2_(std::move(opt)) {
    if (mode_ == FourPointMode::G4_RT) {
      if (!(opt2_.dr > 0.0)) throw std::runtime_error(type_name_ + ": dr must be > 0");
      if (!(opt2_.r_max > 0.0)) throw std::runtime_error(type_name_ + ": r_max must be > 0");
    }
    if ((mode_ == FourPointMode::S4_QT || mode_ == FourPointMode::XI4_T) && opt2_.q_list.empty()) {
      throw std::runtime_error(type_name_ + ": q_list must not be empty");
    }
  }

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
    if (mode_ == FourPointMode::G4_RT) {
      od.columns = {"lag", "time", "r_lo", "r_hi", "r_mid", "g4", "count_pairs", "Q_lag", "mean_dtimestep"};
    } else if (mode_ == FourPointMode::S4_QT) {
      od.columns = {"lag", "time", "q", "S4", "Q_lag", "mean_dtimestep"};
    } else {
      od.columns = {"lag", "time", "xi4", "S4_q0", "fit_intercept", "fit_slope", "n_q_fit", "Q_lag", "mean_dtimestep"};
    }
    md.outputs.push_back(std::move(od));
    return md;
  }

private:
  void finalize_buffered() override {
    if (frames_.empty()) return;
    const auto lags = build_lag_list(frames_.size(), opt2_.lag, true);
    std::vector<std::int64_t> timesteps(frames_.size(), 0);
    for (std::size_t i = 0; i < frames_.size(); ++i) timesteps[i] = frames_[i].timestep;

    std::vector<double> qlag(lags.size(), 0.0);
    std::vector<std::size_t> n_orig(lags.size(), 0);
    for (std::size_t li = 0; li < lags.size(); ++li) {
      const std::size_t lag = lags[li];
      for (std::size_t o = 0; o + lag < frames_.size(); o += static_cast<std::size_t>(opt2_.lag.origin_stride)) {
        const auto& f0 = frames_[o];
        const auto& f1 = frames_[o + lag];
        double qsum = 0.0;
        for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
          const Vec3 dr = f1.pos[i] - f0.pos[i];
          qsum += (norm(dr) <= opt2_.overlap_a) ? 1.0 : 0.0;
        }
        qlag[li] += qsum / static_cast<double>(sel_.idx.size());
        ++n_orig[li];
      }
      if (n_orig[li] > 0) qlag[li] /= static_cast<double>(n_orig[li]);
    }

    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      if (mode_ == FourPointMode::G4_RT) {
        const std::size_t nbins = static_cast<std::size_t>(std::ceil(opt2_.r_max / opt2_.dr));
        ofs << "# PILOTS: connected four-point radial correlation of overlap mobility\n";
        ofs << "# columns: lag time r_lo r_hi r_mid g4 count_pairs Q_lag mean_dtimestep\n";
        for (std::size_t li = 0; li < lags.size(); ++li) {
          const std::size_t lag = lags[li];
          std::vector<double> sum(nbins, 0.0);
          std::vector<std::size_t> cnt(nbins, 0);
          for (std::size_t o = 0; o + lag < frames_.size(); o += static_cast<std::size_t>(opt2_.lag.origin_stride)) {
            const auto& f0 = frames_[o];
            const auto& f1 = frames_[o + lag];
            std::vector<double> dw(sel_.idx.size(), 0.0);
            for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
              const Vec3 dr = f1.pos[i] - f0.pos[i];
              const double wi = (norm(dr) <= opt2_.overlap_a) ? 1.0 : 0.0;
              dw[i] = wi - qlag[li];
            }
            for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
              for (std::size_t j = i + 1; j < sel_.idx.size(); ++j) {
                const Vec3 rij = mic_vec(f0.pos[i], f0.pos[j], f0.lx, f0.ly, f0.lz);
                const double r = norm(rij);
                if (r >= opt2_.r_max) continue;
                const std::size_t b = static_cast<std::size_t>(std::floor(r / opt2_.dr));
                if (b < nbins) {
                  sum[b] += dw[i] * dw[j];
                  ++cnt[b];
                }
              }
            }
          }
          const double dt = mean_dtimestep_for_lag(timesteps, lag, static_cast<std::size_t>(opt2_.lag.origin_stride));
          for (std::size_t b = 0; b < nbins; ++b) {
            const double rlo = static_cast<double>(b) * opt2_.dr;
            const double rhi = rlo + opt2_.dr;
            const double rmid = 0.5 * (rlo + rhi);
            const double g4 = (cnt[b] > 0) ? (sum[b] / static_cast<double>(cnt[b])) : 0.0;
            ofs << lags[li] << ' '
                << std::setprecision(17) << dt << ' '
                << std::setprecision(17) << rlo << ' '
                << std::setprecision(17) << rhi << ' '
                << std::setprecision(17) << rmid << ' '
                << std::setprecision(17) << g4 << ' '
                << cnt[b] << ' '
                << std::setprecision(17) << qlag[li] << ' '
                << std::setprecision(17) << dt << '\n';
          }
        }
      } else if (mode_ == FourPointMode::S4_QT) {
        ofs << "# PILOTS: four-point structure factor of overlap mobility\n";
        ofs << "# columns: lag time q S4 Q_lag mean_dtimestep\n";
        for (std::size_t li = 0; li < lags.size(); ++li) {
          const std::size_t lag = lags[li];
          std::vector<double> s4(opt2_.q_list.size(), 0.0);
          std::size_t norig = 0;
          for (std::size_t o = 0; o + lag < frames_.size(); o += static_cast<std::size_t>(opt2_.lag.origin_stride)) {
            const auto& f0 = frames_[o];
            const auto& f1 = frames_[o + lag];
            std::vector<double> dw(sel_.idx.size(), 0.0);
            double self2 = 0.0;
            for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
              const Vec3 dr = f1.pos[i] - f0.pos[i];
              const double wi = (norm(dr) <= opt2_.overlap_a) ? 1.0 : 0.0;
              dw[i] = wi - qlag[li];
              self2 += dw[i] * dw[i];
            }
            for (std::size_t qi = 0; qi < opt2_.q_list.size(); ++qi) {
              double val = self2;
              for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
                for (std::size_t j = i + 1; j < sel_.idx.size(); ++j) {
                  const Vec3 rij = mic_vec(f0.pos[i], f0.pos[j], f0.lx, f0.ly, f0.lz);
                  val += 2.0 * dw[i] * dw[j] * j0(opt2_.q_list[qi] * norm(rij));
                }
              }
              s4[qi] += val / static_cast<double>(sel_.idx.size());
            }
            ++norig;
          }
          const double dt = mean_dtimestep_for_lag(timesteps, lag, static_cast<std::size_t>(opt2_.lag.origin_stride));
          for (std::size_t qi = 0; qi < opt2_.q_list.size(); ++qi) {
            const double sval = (norig > 0) ? (s4[qi] / static_cast<double>(norig)) : 0.0;
            ofs << lags[li] << ' '
                << std::setprecision(17) << dt << ' '
                << std::setprecision(17) << opt2_.q_list[qi] << ' '
                << std::setprecision(17) << sval << ' '
                << std::setprecision(17) << qlag[li] << ' '
                << std::setprecision(17) << dt << '\n';
          }
        }
      } else {
        ofs << "# PILOTS: dynamic four-point correlation length from low-q Ornstein-Zernike fit\n";
        ofs << "# columns: lag time xi4 S4_q0 fit_intercept fit_slope n_q_fit Q_lag mean_dtimestep\n";
        for (std::size_t li = 0; li < lags.size(); ++li) {
          const std::size_t lag = lags[li];
          std::vector<double> s4(opt2_.q_list.size(), 0.0);
          std::size_t norig = 0;
          for (std::size_t o = 0; o + lag < frames_.size(); o += static_cast<std::size_t>(opt2_.lag.origin_stride)) {
            const auto& f0 = frames_[o];
            const auto& f1 = frames_[o + lag];
            std::vector<double> dw(sel_.idx.size(), 0.0);
            double self2 = 0.0;
            for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
              const Vec3 dr = f1.pos[i] - f0.pos[i];
              const double wi = (norm(dr) <= opt2_.overlap_a) ? 1.0 : 0.0;
              dw[i] = wi - qlag[li];
              self2 += dw[i] * dw[i];
            }
            for (std::size_t qi = 0; qi < opt2_.q_list.size(); ++qi) {
              double val = self2;
              for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
                for (std::size_t j = i + 1; j < sel_.idx.size(); ++j) {
                  const Vec3 rij = mic_vec(f0.pos[i], f0.pos[j], f0.lx, f0.ly, f0.lz);
                  val += 2.0 * dw[i] * dw[j] * j0(opt2_.q_list[qi] * norm(rij));
                }
              }
              s4[qi] += val / static_cast<double>(sel_.idx.size());
            }
            ++norig;
          }
          for (double& v : s4) v = (norig > 0) ? (v / static_cast<double>(norig)) : 0.0;
          std::vector<std::pair<double,double>> fitpts;
          fitpts.reserve(opt2_.q_list.size());
          for (std::size_t qi = 0; qi < opt2_.q_list.size(); ++qi) {
            const double q = opt2_.q_list[qi];
            const double s = s4[qi];
            if (s > 0.0) fitpts.push_back({q * q, 1.0 / s});
          }
          std::sort(fitpts.begin(), fitpts.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
          if (fitpts.size() > opt2_.fit_n_q) fitpts.resize(opt2_.fit_n_q);
          double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
          for (const auto& p : fitpts) {
            sx += p.first; sy += p.second; sxx += p.first * p.first; sxy += p.first * p.second;
          }
          double a = 0.0, b = 0.0;
          const double n = static_cast<double>(fitpts.size());
          if (fitpts.size() >= 2) {
            const double den = n * sxx - sx * sx;
            if (std::abs(den) > 0.0) {
              b = (n * sxy - sx * sy) / den;
              a = (sy - b * sx) / n;
            }
          } else if (fitpts.size() == 1) {
            a = fitpts[0].second;
            b = 0.0;
          }
          double xi = 0.0;
          if (a > 0.0 && b > 0.0) xi = std::sqrt(b / a);
          const double dt = mean_dtimestep_for_lag(timesteps, lag, static_cast<std::size_t>(opt2_.lag.origin_stride));
          const double s40 = !s4.empty() ? s4.front() : 0.0;
          ofs << lags[li] << ' '
              << std::setprecision(17) << dt << ' '
              << std::setprecision(17) << xi << ' '
              << std::setprecision(17) << s40 << ' '
              << std::setprecision(17) << a << ' '
              << std::setprecision(17) << b << ' '
              << fitpts.size() << ' '
              << std::setprecision(17) << qlag[li] << ' '
              << std::setprecision(17) << dt << '\n';
        }
      }
    });
  }

  std::string output_path_;
  FourPointMode mode_;
  Options opt2_;
};

class StringMeasure final : public BufferedSelectedPositionsMeasure {
public:
  enum class Mode { Summary, Distribution };
  struct Options {
    BufferedSelectedPositionsMeasure::Options base;
    LagOptions lag;
    double mobile_cutoff = 1.0;
    double link_cutoff = 0.3;
  };

  StringMeasure(std::string type_name,
                std::string instance_name,
                std::string output_path,
                SelectionView sel,
                SelectionView drift_sel,
                Mode mode,
                Options opt)
      : BufferedSelectedPositionsMeasure(type_name, instance_name, sel, drift_sel, opt.base),
        output_path_(std::move(output_path)),
        mode_(mode),
        opt2_(std::move(opt)) {
    if (!(opt2_.mobile_cutoff > 0.0)) throw std::runtime_error(type_name_ + ": mobile_cutoff must be > 0");
    if (!(opt2_.link_cutoff > 0.0)) throw std::runtime_error(type_name_ + ": link_cutoff must be > 0");
  }

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
    if (mode_ == Mode::Summary) {
      od.columns = {"lag", "time", "mean_string_size", "max_string_size", "n_strings", "n_mobile", "mobile_fraction", "mean_dtimestep"};
    } else {
      od.columns = {"lag", "time", "string_size", "probability", "count", "mean_dtimestep"};
    }
    md.outputs.push_back(std::move(od));
    return md;
  }

private:
  void finalize_buffered() override {
    if (frames_.empty()) return;
    const auto lags = build_lag_list(frames_.size(), opt2_.lag, false);
    std::vector<std::int64_t> timesteps(frames_.size(), 0);
    for (std::size_t i = 0; i < frames_.size(); ++i) timesteps[i] = frames_[i].timestep;

    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      if (mode_ == Mode::Summary) {
        ofs << "# PILOTS: string-like cooperative motion summary\n";
        ofs << "# columns: lag time mean_string_size max_string_size n_strings n_mobile mobile_fraction mean_dtimestep\n";
      } else {
        ofs << "# PILOTS: string length distribution\n";
        ofs << "# columns: lag time string_size probability count mean_dtimestep\n";
      }
      for (const std::size_t lag : lags) {
        double sum_mean = 0.0, sum_max = 0.0, sum_nstr = 0.0, sum_nmob = 0.0;
        std::size_t norig = 0;
        std::unordered_map<std::size_t, std::size_t> hist;
        std::size_t hist_total = 0;
        for (std::size_t o = 0; o + lag < frames_.size(); o += static_cast<std::size_t>(opt2_.lag.origin_stride)) {
          const auto& f0 = frames_[o];
          const auto& f1 = frames_[o + lag];
          std::vector<std::size_t> mobile;
          mobile.reserve(sel_.idx.size());
          for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
            if (norm(f1.pos[i] - f0.pos[i]) >= opt2_.mobile_cutoff) mobile.push_back(i);
          }
          const std::size_t nm = mobile.size();
          if (nm == 0) {
            ++norig;
            continue;
          }
          DSU dsu(nm);
          for (std::size_t a = 0; a < nm; ++a) {
            for (std::size_t b = a + 1; b < nm; ++b) {
              const std::size_t i = mobile[a];
              const std::size_t j = mobile[b];
              const double d1 = norm(mic_vec(f0.pos[i], f1.pos[j], f0.lx, f0.ly, f0.lz));
              const double d2 = norm(mic_vec(f0.pos[j], f1.pos[i], f0.lx, f0.ly, f0.lz));
              if (std::min(d1, d2) <= opt2_.link_cutoff) dsu.unite(a, b);
            }
          }
          std::unordered_map<std::size_t, std::size_t> sizes;
          for (std::size_t a = 0; a < nm; ++a) ++sizes[dsu.find(a)];
          double mean_sz = 0.0;
          std::size_t max_sz = 0;
          for (const auto& kv : sizes) {
            mean_sz += static_cast<double>(kv.second);
            max_sz = std::max(max_sz, kv.second);
            ++hist[kv.second];
            ++hist_total;
          }
          mean_sz /= static_cast<double>(sizes.size());
          sum_mean += mean_sz;
          sum_max += static_cast<double>(max_sz);
          sum_nstr += static_cast<double>(sizes.size());
          sum_nmob += static_cast<double>(nm);
          ++norig;
        }
        const double dt = mean_dtimestep_for_lag(timesteps, lag, static_cast<std::size_t>(opt2_.lag.origin_stride));
        if (mode_ == Mode::Summary) {
          const double denom = (norig > 0) ? static_cast<double>(norig) : 1.0;
          ofs << lag << ' '
              << std::setprecision(17) << dt << ' '
              << std::setprecision(17) << (sum_mean / denom) << ' '
              << std::setprecision(17) << (sum_max / denom) << ' '
              << std::setprecision(17) << (sum_nstr / denom) << ' '
              << std::setprecision(17) << (sum_nmob / denom) << ' '
              << std::setprecision(17) << ((sum_nmob / denom) / static_cast<double>(sel_.idx.size())) << ' '
              << std::setprecision(17) << dt << '\n';
        } else {
          std::vector<std::pair<std::size_t,std::size_t>> items(hist.begin(), hist.end());
          std::sort(items.begin(), items.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
          for (const auto& kv : items) {
            const double p = (hist_total > 0) ? (static_cast<double>(kv.second) / static_cast<double>(hist_total)) : 0.0;
            ofs << lag << ' '
                << std::setprecision(17) << dt << ' '
                << kv.first << ' '
                << std::setprecision(17) << p << ' '
                << kv.second << ' '
                << std::setprecision(17) << dt << '\n';
          }
        }
      }
    });
  }

  std::string output_path_;
  Mode mode_;
  Options opt2_;
};

class MobileClusterMeasure final : public BufferedSelectedPositionsMeasure {
public:
  struct Options {
    BufferedSelectedPositionsMeasure::Options base;
    LagOptions lag;
    double mobile_cutoff = 1.0;
    double cluster_cutoff = 1.0;
  };

  MobileClusterMeasure(std::string instance_name,
                       std::string output_path,
                       SelectionView sel,
                       SelectionView drift_sel,
                       Options opt)
      : BufferedSelectedPositionsMeasure("mobile_cluster", instance_name, sel, drift_sel, opt.base),
        output_path_(std::move(output_path)),
        opt2_(std::move(opt)) {
    if (!(opt2_.mobile_cutoff > 0.0)) throw std::runtime_error("mobile_cluster: mobile_cutoff must be > 0");
    if (!(opt2_.cluster_cutoff > 0.0)) throw std::runtime_error("mobile_cluster: cluster_cutoff must be > 0");
  }

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
    od.columns = {"lag", "time", "mean_cluster_size", "max_cluster_size", "n_clusters", "n_mobile", "mobile_fraction", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));
    return md;
  }

private:
  void finalize_buffered() override {
    if (frames_.empty()) return;
    const auto lags = build_lag_list(frames_.size(), opt2_.lag, false);
    std::vector<std::int64_t> timesteps(frames_.size(), 0);
    for (std::size_t i = 0; i < frames_.size(); ++i) timesteps[i] = frames_[i].timestep;

    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: clusters of mobile particles\n";
      ofs << "# columns: lag time mean_cluster_size max_cluster_size n_clusters n_mobile mobile_fraction mean_dtimestep\n";
      for (const std::size_t lag : lags) {
        double sum_mean = 0.0, sum_max = 0.0, sum_ncl = 0.0, sum_nmob = 0.0;
        std::size_t norig = 0;
        for (std::size_t o = 0; o + lag < frames_.size(); o += static_cast<std::size_t>(opt2_.lag.origin_stride)) {
          const auto& f0 = frames_[o];
          const auto& f1 = frames_[o + lag];
          std::vector<std::size_t> mobile;
          for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
            if (norm(f1.pos[i] - f0.pos[i]) >= opt2_.mobile_cutoff) mobile.push_back(i);
          }
          const std::size_t nm = mobile.size();
          if (nm == 0) { ++norig; continue; }
          DSU dsu(nm);
          for (std::size_t a = 0; a < nm; ++a) {
            for (std::size_t b = a + 1; b < nm; ++b) {
              const Vec3 rij = mic_vec(f0.pos[mobile[a]], f0.pos[mobile[b]], f0.lx, f0.ly, f0.lz);
              if (norm(rij) <= opt2_.cluster_cutoff) dsu.unite(a, b);
            }
          }
          std::unordered_map<std::size_t, std::size_t> sizes;
          for (std::size_t a = 0; a < nm; ++a) ++sizes[dsu.find(a)];
          double mean_sz = 0.0;
          std::size_t max_sz = 0;
          for (const auto& kv : sizes) {
            mean_sz += static_cast<double>(kv.second);
            max_sz = std::max(max_sz, kv.second);
          }
          mean_sz /= static_cast<double>(sizes.size());
          sum_mean += mean_sz;
          sum_max += static_cast<double>(max_sz);
          sum_ncl += static_cast<double>(sizes.size());
          sum_nmob += static_cast<double>(nm);
          ++norig;
        }
        const double dt = mean_dtimestep_for_lag(timesteps, lag, static_cast<std::size_t>(opt2_.lag.origin_stride));
        const double denom = (norig > 0) ? static_cast<double>(norig) : 1.0;
        ofs << lag << ' '
            << std::setprecision(17) << dt << ' '
            << std::setprecision(17) << (sum_mean / denom) << ' '
            << std::setprecision(17) << (sum_max / denom) << ' '
            << std::setprecision(17) << (sum_ncl / denom) << ' '
            << std::setprecision(17) << (sum_nmob / denom) << ' '
            << std::setprecision(17) << ((sum_nmob / denom) / static_cast<double>(sel_.idx.size())) << ' '
            << std::setprecision(17) << dt << '\n';
      }
    });
  }

  std::string output_path_;
  Options opt2_;
};

class ExcitationMapMeasure final : public BufferedSelectedPositionsMeasure {
public:
  struct Options {
    BufferedSelectedPositionsMeasure::Options base;
    std::size_t lag_frames = 1;
    std::size_t origin_stride = 1;
    double mobile_cutoff = 1.0;
    std::string id_field = "id";
  };

  ExcitationMapMeasure(std::string instance_name,
                       std::string output_path,
                       SelectionView sel,
                       SelectionView drift_sel,
                       std::vector<std::int64_t> id_by_atom,
                       Options opt)
      : BufferedSelectedPositionsMeasure("excitation_map", instance_name, sel, drift_sel, opt.base),
        output_path_(std::move(output_path)),
        id_by_atom_(std::move(id_by_atom)),
        opt2_(std::move(opt)) {
    if (opt2_.lag_frames == 0) throw std::runtime_error("excitation_map: lag_frames must be >= 1");
    if (opt2_.origin_stride == 0) throw std::runtime_error("excitation_map: origin_stride must be >= 1");
  }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();
    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.columns = {"atom", "id", "p_exc", "n_exc", "n_samples", "mean_x0", "mean_y0", "mean_z0"};
    md.outputs.push_back(std::move(od));
    return md;
  }

private:
  void finalize_buffered() override {
    if (frames_.size() <= opt2_.lag_frames) return;
    std::vector<std::size_t> n_exc(sel_.idx.size(), 0);
    std::vector<std::size_t> n_samp(sel_.idx.size(), 0);
    std::vector<double> mx(sel_.idx.size(), 0.0), my(sel_.idx.size(), 0.0), mz(sel_.idx.size(), 0.0);
    for (std::size_t o = 0; o + opt2_.lag_frames < frames_.size(); o += opt2_.origin_stride) {
      const auto& f0 = frames_[o];
      const auto& f1 = frames_[o + opt2_.lag_frames];
      for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
        if (norm(f1.pos[i] - f0.pos[i]) >= opt2_.mobile_cutoff) ++n_exc[i];
        ++n_samp[i];
        mx[i] += f0.pos[i].x;
        my[i] += f0.pos[i].y;
        mz[i] += f0.pos[i].z;
      }
    }
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: excitation propensity map for a fixed lag\n";
      ofs << "# columns: atom id p_exc n_exc n_samples mean_x0 mean_y0 mean_z0\n";
      for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
        const std::size_t atom = sel_.idx[p];
        const std::size_t ns = n_samp[p];
        const double inv = (ns > 0) ? (1.0 / static_cast<double>(ns)) : 0.0;
        ofs << atom << ' '
            << ((atom < id_by_atom_.size()) ? id_by_atom_[atom] : static_cast<std::int64_t>(atom)) << ' '
            << std::setprecision(17) << ((ns > 0) ? (static_cast<double>(n_exc[p]) * inv) : 0.0) << ' '
            << n_exc[p] << ' '
            << ns << ' '
            << std::setprecision(17) << (mx[p] * inv) << ' '
            << std::setprecision(17) << (my[p] * inv) << ' '
            << std::setprecision(17) << (mz[p] * inv) << '\n';
      }
    });
  }

  std::string output_path_;
  std::vector<std::int64_t> id_by_atom_;
  Options opt2_;
};

class FacilitationACFMeasure final : public BufferedSelectedPositionsMeasure {
public:
  struct Options {
    BufferedSelectedPositionsMeasure::Options base;
    std::size_t window_lag_frames = 1;
    std::size_t max_sep_frames = 1;
    std::size_t sep_stride = 1;
    std::size_t origin_stride = 1;
    double mobile_cutoff = 1.0;
    double neighbor_cutoff = 1.0;
  };

  FacilitationACFMeasure(std::string instance_name,
                         std::string output_path,
                         SelectionView sel,
                         SelectionView drift_sel,
                         Options opt)
      : BufferedSelectedPositionsMeasure("facilitation_acf", instance_name, sel, drift_sel, opt.base),
        output_path_(std::move(output_path)),
        opt2_(std::move(opt)) {
    if (opt2_.window_lag_frames == 0) throw std::runtime_error("facilitation_acf: window_lag_frames must be >= 1");
    if (opt2_.sep_stride == 0 || opt2_.origin_stride == 0) {
      throw std::runtime_error("facilitation_acf: sep_stride and origin_stride must be >= 1");
    }
  }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();
    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.x_axis = "sep";
    od.x_unit = "frames";
    od.columns = {"sep", "time", "p_neighbor_excited", "p_later_excited", "joint_prob", "covariance", "facilitation_ratio", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));
    return md;
  }

private:
  void finalize_buffered() override {
    if (frames_.size() <= opt2_.window_lag_frames) return;
    std::vector<std::int64_t> timesteps(frames_.size(), 0);
    for (std::size_t i = 0; i < frames_.size(); ++i) timesteps[i] = frames_[i].timestep;

    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: facilitation cross-correlation between earlier excited neighborhoods and later excitations\n";
      ofs << "# columns: sep time p_neighbor_excited p_later_excited joint_prob covariance facilitation_ratio mean_dtimestep\n";
      for (std::size_t sep = 0; sep <= opt2_.max_sep_frames; sep += opt2_.sep_stride) {
        double sum_n = 0.0, sum_e = 0.0, sum_ne = 0.0;
        std::size_t count = 0;
        for (std::size_t o = 0; o + opt2_.window_lag_frames + sep + opt2_.window_lag_frames < frames_.size(); o += opt2_.origin_stride) {
          const auto& f0 = frames_[o];
          const auto& f1 = frames_[o + opt2_.window_lag_frames];
          const auto& g0 = frames_[o + sep];
          const auto& g1 = frames_[o + sep + opt2_.window_lag_frames];
          std::vector<unsigned char> e0(sel_.idx.size(), 0), e1(sel_.idx.size(), 0), nbh(sel_.idx.size(), 0);
          for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
            e0[i] = (norm(f1.pos[i] - f0.pos[i]) >= opt2_.mobile_cutoff) ? 1 : 0;
            e1[i] = (norm(g1.pos[i] - g0.pos[i]) >= opt2_.mobile_cutoff) ? 1 : 0;
          }
          for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
            bool has_excited_neighbor = false;
            for (std::size_t j = 0; j < sel_.idx.size(); ++j) {
              if (i == j || e0[j] == 0) continue;
              const Vec3 rij = mic_vec(f0.pos[i], f0.pos[j], f0.lx, f0.ly, f0.lz);
              if (norm(rij) <= opt2_.neighbor_cutoff) {
                has_excited_neighbor = true;
                break;
              }
            }
            nbh[i] = has_excited_neighbor ? 1 : 0;
            sum_n += static_cast<double>(nbh[i]);
            sum_e += static_cast<double>(e1[i]);
            sum_ne += static_cast<double>(nbh[i] * e1[i]);
            ++count;
          }
        }
        const double inv = (count > 0) ? (1.0 / static_cast<double>(count)) : 0.0;
        const double pn = sum_n * inv;
        const double pe = sum_e * inv;
        const double jne = sum_ne * inv;
        const double cov = jne - pn * pe;
        const double ratio = (pn > 0.0 && pe > 0.0) ? (jne / (pn * pe)) : 0.0;
        const double dt = mean_dtimestep_for_lag(timesteps, sep, opt2_.origin_stride);
        ofs << sep << ' '
            << std::setprecision(17) << dt << ' '
            << std::setprecision(17) << pn << ' '
            << std::setprecision(17) << pe << ' '
            << std::setprecision(17) << jne << ' '
            << std::setprecision(17) << cov << ' '
            << std::setprecision(17) << ratio << ' '
            << std::setprecision(17) << dt << '\n';
      }
    });
  }

  std::string output_path_;
  Options opt2_;
};

MeasureCapabilities buffered_pos_caps(const IniConfig& cfg,
                                      const std::string& section,
                                      const std::string& instance,
                                      const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(false));
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.requires_identity_consistent = true;
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  if (cfg.has_key(section, "group_b")) caps.group_refs.push_back(cfg.get_string(section, "group_b"));
  if (remove_drift) caps.group_refs.push_back(cfg.get_string(section, "drift_group", std::optional<std::string>("all")));
  return caps;
}

std::unique_ptr<IMeasure> create_four_point_like(const IniConfig& cfg,
                                                 const std::string& section,
                                                 const std::string& instance,
                                                 const MeasureBuildEnv& env,
                                                 const SystemContext& sysctx,
                                                 std::string type_name,
                                                 FourPointMode mode) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error(type_name + " factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(false));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
  SelectionView drift_sel = remove_drift ? get_static_group_view(*env.selection_provider, frame0, drift_group, type_name)
                                         : get_static_group_view(*env.selection_provider, frame0, "all", type_name);
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>(type_name + std::string(".dat")))).lexically_normal();

  FourPointMeasure::Options opt;
  opt.base.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.base.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.base.range.dry_run = env.dry_run;
  opt.base.remove_drift = remove_drift;
  opt.lag = lag_options_from_config(cfg, section);
  opt.overlap_a = cfg.get_double(section, "overlap_a", std::optional<double>(1.0));
  opt.r_max = cfg.get_double(section, "r_max", std::optional<double>(0.0));
  opt.dr = cfg.get_double(section, "dr", std::optional<double>(0.1));
  if (cfg.has_key(section, "q_list")) opt.q_list = measure_ext::parse_double_list(cfg, section, "q_list");
  opt.fit_n_q = static_cast<std::size_t>(cfg.get_int64(section, "fit_n_q", std::optional<std::int64_t>(3)));
  return std::make_unique<FourPointMeasure>(type_name, instance, out.string(), sel, drift_sel, mode, opt);
}

std::unique_ptr<IMeasure> g4_rt_create(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& instance,
                                       const MeasureBuildEnv& env,
                                       const SystemContext& sysctx) {
  return create_four_point_like(cfg, section, instance, env, sysctx, "g4_rt", FourPointMode::G4_RT);
}

std::unique_ptr<IMeasure> s4_qt_create(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& instance,
                                       const MeasureBuildEnv& env,
                                       const SystemContext& sysctx) {
  return create_four_point_like(cfg, section, instance, env, sysctx, "s4_qt", FourPointMode::S4_QT);
}

std::unique_ptr<IMeasure> xi4_t_create(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& instance,
                                       const MeasureBuildEnv& env,
                                       const SystemContext& sysctx) {
  return create_four_point_like(cfg, section, instance, env, sysctx, "xi4_t", FourPointMode::XI4_T);
}

std::unique_ptr<IMeasure> create_string_like(const IniConfig& cfg,
                                             const std::string& section,
                                             const std::string& instance,
                                             const MeasureBuildEnv& env,
                                             const SystemContext& sysctx,
                                             std::string type_name,
                                             StringMeasure::Mode mode) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error(type_name + " factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(false));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
  SelectionView drift_sel = remove_drift ? get_static_group_view(*env.selection_provider, frame0, drift_group, type_name)
                                         : get_static_group_view(*env.selection_provider, frame0, "all", type_name);
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>(type_name + std::string(".dat")))).lexically_normal();

  StringMeasure::Options opt;
  opt.base.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.base.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.base.range.dry_run = env.dry_run;
  opt.base.remove_drift = remove_drift;
  opt.lag = lag_options_from_config(cfg, section);
  opt.mobile_cutoff = cfg.get_double(section, "mobile_cutoff", std::optional<double>(1.0));
  opt.link_cutoff = cfg.get_double(section, "link_cutoff", std::optional<double>(0.3));
  return std::make_unique<StringMeasure>(type_name, instance, out.string(), sel, drift_sel, mode, opt);
}

std::unique_ptr<IMeasure> string_motion_create(const IniConfig& cfg,
                                               const std::string& section,
                                               const std::string& instance,
                                               const MeasureBuildEnv& env,
                                               const SystemContext& sysctx) {
  return create_string_like(cfg, section, instance, env, sysctx, "string_motion", StringMeasure::Mode::Summary);
}

std::unique_ptr<IMeasure> string_length_dist_create(const IniConfig& cfg,
                                                    const std::string& section,
                                                    const std::string& instance,
                                                    const MeasureBuildEnv& env,
                                                    const SystemContext& sysctx) {
  return create_string_like(cfg, section, instance, env, sysctx, "string_length_dist", StringMeasure::Mode::Distribution);
}

std::unique_ptr<IMeasure> mobile_cluster_create(const IniConfig& cfg,
                                                const std::string& section,
                                                const std::string& instance,
                                                const MeasureBuildEnv& env,
                                                const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("mobile_cluster factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "mobile_cluster");
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(false));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
  SelectionView drift_sel = remove_drift ? get_static_group_view(*env.selection_provider, frame0, drift_group, "mobile_cluster")
                                         : get_static_group_view(*env.selection_provider, frame0, "all", "mobile_cluster");
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("mobile_cluster.dat"))).lexically_normal();

  MobileClusterMeasure::Options opt;
  opt.base.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.base.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.base.range.dry_run = env.dry_run;
  opt.base.remove_drift = remove_drift;
  opt.lag = lag_options_from_config(cfg, section);
  opt.mobile_cutoff = cfg.get_double(section, "mobile_cutoff", std::optional<double>(1.0));
  opt.cluster_cutoff = cfg.get_double(section, "cluster_cutoff", std::optional<double>(1.0));
  return std::make_unique<MobileClusterMeasure>(instance, out.string(), sel, drift_sel, opt);
}

MeasureCapabilities excitation_map_caps(const IniConfig& cfg,
                                        const std::string& section,
                                        const std::string& instance,
                                        const MeasureBuildEnv& env) {
  MeasureCapabilities caps = buffered_pos_caps(cfg, section, instance, env);
  caps.requires_intfields.push_back(cfg.get_string(section, "id_field", std::optional<std::string>("id")));
  return caps;
}

std::unique_ptr<IMeasure> excitation_map_create(const IniConfig& cfg,
                                                const std::string& section,
                                                const std::string& instance,
                                                const MeasureBuildEnv& env,
                                                const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("excitation_map factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "excitation_map");
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(false));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
  SelectionView drift_sel = remove_drift ? get_static_group_view(*env.selection_provider, frame0, drift_group, "excitation_map")
                                         : get_static_group_view(*env.selection_provider, frame0, "all", "excitation_map");
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("excitation_map.dat"))).lexically_normal();
  const std::string id_field = cfg.get_string(section, "id_field", std::optional<std::string>("id"));
  const auto ids = integer_like_field_to_i64(frame0, id_field, true);

  ExcitationMapMeasure::Options opt;
  opt.base.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.base.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.base.range.dry_run = env.dry_run;
  opt.base.remove_drift = remove_drift;
  opt.lag_frames = static_cast<std::size_t>(cfg.get_int64(section, "lag_frames", std::optional<std::int64_t>(1)));
  opt.origin_stride = static_cast<std::size_t>(cfg.get_int64(section, "origin_stride", std::optional<std::int64_t>(1)));
  opt.mobile_cutoff = cfg.get_double(section, "mobile_cutoff", std::optional<double>(1.0));
  opt.id_field = id_field;
  return std::make_unique<ExcitationMapMeasure>(instance, out.string(), sel, drift_sel, ids, opt);
}

std::unique_ptr<IMeasure> facilitation_acf_create(const IniConfig& cfg,
                                                  const std::string& section,
                                                  const std::string& instance,
                                                  const MeasureBuildEnv& env,
                                                  const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("facilitation_acf factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "facilitation_acf");
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(false));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
  SelectionView drift_sel = remove_drift ? get_static_group_view(*env.selection_provider, frame0, drift_group, "facilitation_acf")
                                         : get_static_group_view(*env.selection_provider, frame0, "all", "facilitation_acf");
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("facilitation_acf.dat"))).lexically_normal();

  FacilitationACFMeasure::Options opt;
  opt.base.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.base.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.base.range.dry_run = env.dry_run;
  opt.base.remove_drift = remove_drift;
  opt.window_lag_frames = static_cast<std::size_t>(cfg.get_int64(section, "window_lag_frames", std::optional<std::int64_t>(1)));
  opt.max_sep_frames = static_cast<std::size_t>(cfg.get_int64(section, "max_sep_frames", std::optional<std::int64_t>(10)));
  opt.sep_stride = static_cast<std::size_t>(cfg.get_int64(section, "sep_stride", std::optional<std::int64_t>(1)));
  opt.origin_stride = static_cast<std::size_t>(cfg.get_int64(section, "origin_stride", std::optional<std::int64_t>(1)));
  opt.mobile_cutoff = cfg.get_double(section, "mobile_cutoff", std::optional<double>(1.0));
  opt.neighbor_cutoff = cfg.get_double(section, "neighbor_cutoff", std::optional<double>(1.0));
  return std::make_unique<FacilitationACFMeasure>(instance, out.string(), sel, drift_sel, opt);
}

static MeasureRegistrar g_reg_g4_rt("g4_rt", &buffered_pos_caps, &g4_rt_create);
static MeasureRegistrar g_reg_s4_qt("s4_qt", &buffered_pos_caps, &s4_qt_create);
static MeasureRegistrar g_reg_xi4_t("xi4_t", &buffered_pos_caps, &xi4_t_create);
static MeasureRegistrar g_reg_string_motion("string_motion", &buffered_pos_caps, &string_motion_create);
static MeasureRegistrar g_reg_string_length_dist("string_length_dist", &buffered_pos_caps, &string_length_dist_create);
static MeasureRegistrar g_reg_mobile_cluster("mobile_cluster", &buffered_pos_caps, &mobile_cluster_create);
static MeasureRegistrar g_reg_excitation_map("excitation_map", &excitation_map_caps, &excitation_map_create);
static MeasureRegistrar g_reg_facilitation_acf("facilitation_acf", &buffered_pos_caps, &facilitation_acf_create);

} // namespace
} // namespace pilots
