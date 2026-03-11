#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
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
using measure_ext::SlabFrameTransform;
using measure_ext::SlabTransformOptions;
using measure_ext::Vec3;
using measure_ext::atom_vec3;
using measure_ext::build_optional_anchor_selection;
using measure_ext::axis1d_name;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::make_slab_frame_transform;
using measure_ext::min_image_difference;
using measure_ext::norm2;
using measure_ext::parse_axis1d;

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

inline double sem_from_sums(double sum, double sumsq, std::size_t n) {
  if (n < 2) return 0.0;
  const double mean = sum / static_cast<double>(n);
  const double sample_var = std::max(0.0, (sumsq - static_cast<double>(n) * mean * mean)
                                              / static_cast<double>(n - 1));
  return std::sqrt(sample_var / static_cast<double>(n));
}

enum class PairEnergyModel {
  None,
  LJ126,
};

inline PairEnergyModel parse_pair_energy_model(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  if (s == "none" || s == "hard" || s == "hard_sphere") return PairEnergyModel::None;
  if (s == "lj" || s == "lj126" || s == "lennard-jones" || s == "lennard_jones") {
    return PairEnergyModel::LJ126;
  }
  throw std::runtime_error("invalid energy_model: '" + s + "' (use none|lj126)");
}

inline const char* pair_energy_model_name(PairEnergyModel m) {
  switch (m) {
    case PairEnergyModel::None: return "none";
    case PairEnergyModel::LJ126: return "lj126";
  }
  return "none";
}

struct ReservoirReferenceOptions {
  bool enabled = true;
  double left_lo_frac = 0.0;
  double left_hi_frac = 0.1;
  bool use_right_window = false;
  double right_lo_frac = 0.9;
  double right_hi_frac = 1.0;
};

struct InsertionOptions {
  RangeOptions range;
  SlabTransformOptions slab;
  ReservoirReferenceOptions reservoir;
  std::size_t n_bins = 0;
  std::size_t samples_plane = 16;
  std::size_t samples_axis = 2;
  std::size_t convergence_refine_factor = 2;
  double probe_radius = 0.0;
  std::string radius_field;
  double occluder_radius = 0.0;
  std::string sigma_field;
  double occluder_sigma = 0.0;
  double probe_sigma = 0.0;
  PairEnergyModel energy_model = PairEnergyModel::None;
  double beta = 1.0;
  std::string epsilon_field;
  double occluder_epsilon = 1.0;
  double probe_epsilon = 1.0;
  double soft_cutoff = 0.0;
};

struct OccluderData {
  std::vector<Vec3> pos;
  std::vector<double> radius;
  std::vector<double> sigma;
  std::vector<double> epsilon;
};

struct FrameInsertionStats {
  double domain_length = 0.0;
  std::size_t samples_per_bin = 0;
  std::vector<std::size_t> n_total;
  std::vector<std::size_t> n_access;
  std::vector<double> weight_sum;
  std::vector<double> h;
  std::vector<double> total_mean;
  std::vector<double> cond_mean;
};

struct WindowedMeans {
  bool valid = false;
  double h_ref = 1.0;
  double total_ref = 1.0;
  double cond_ref = 1.0;
};

inline double clamp_unit(double x) {
  return std::clamp(x, 0.0, 1.0);
}

inline void validate_window(double lo, double hi, const std::string& what) {
  if (!(lo >= 0.0 && lo <= 1.0 && hi >= 0.0 && hi <= 1.0 && hi > lo)) {
    throw std::runtime_error(what + " must satisfy 0 <= lo < hi <= 1");
  }
}

inline std::vector<std::size_t> reservoir_bins(std::size_t n_bins,
                                               const ReservoirReferenceOptions& opt) {
  std::vector<std::size_t> bins;
  bins.reserve(n_bins);
  for (std::size_t b = 0; b < n_bins; ++b) {
    const double center = (static_cast<double>(b) + 0.5) / static_cast<double>(n_bins);
    const bool in_left = (center >= opt.left_lo_frac && center <= opt.left_hi_frac);
    const bool in_right = opt.use_right_window && (center >= opt.right_lo_frac && center <= opt.right_hi_frac);
    if (in_left || in_right) bins.push_back(b);
  }
  if (bins.empty()) {
    throw std::runtime_error("reservoir reference window does not contain any bins; adjust n_bins or reservoir window fractions");
  }
  std::sort(bins.begin(), bins.end());
  bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
  return bins;
}

inline double average_over_bins(const std::vector<double>& vals,
                                const std::vector<std::size_t>& bins,
                                const std::vector<unsigned char>* valid = nullptr) {
  double sum = 0.0;
  std::size_t n = 0;
  for (const std::size_t b : bins) {
    if (b >= vals.size()) continue;
    if (valid && (b >= valid->size() || (*valid)[b] == 0)) continue;
    sum += vals[b];
    ++n;
  }
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  return sum / static_cast<double>(n);
}

inline double sigma_from_radius(double r) {
  if (r <= 0.0) return 0.0;
  return (2.0 * r) / std::pow(2.0, 1.0 / 6.0);
}

inline double sigma_fallback(double contact_distance) {
  if (contact_distance <= 0.0) return 0.0;
  return contact_distance / std::pow(2.0, 1.0 / 6.0);
}

inline double pair_soft_energy(double r,
                               double sigma_eff,
                               double eps_eff,
                               double soft_cutoff,
                               PairEnergyModel model) {
  if (model == PairEnergyModel::None) return 0.0;
  if (!(sigma_eff > 0.0) || !(eps_eff > 0.0)) return 0.0;
  if (r <= 0.0) return std::numeric_limits<double>::infinity();
  if (soft_cutoff > 0.0 && r > soft_cutoff) return 0.0;
  const double sr = sigma_eff / r;
  const double sr2 = sr * sr;
  const double sr6 = sr2 * sr2 * sr2;
  const double sr12 = sr6 * sr6;
  return 4.0 * eps_eff * (sr12 - sr6);
}

OccluderData build_occluder_data(const Frame& frame,
                                 const SelectionView& occluders,
                                 const InsertionOptions& opt) {
  if (occluders.idx.empty()) {
    throw std::runtime_error("insertion-profile evaluation requires a non-empty occluder selection");
  }
  const auto rad_field = (!opt.radius_field.empty())
      ? std::optional<std::span<const double>>(frame.require_dfield(opt.radius_field))
      : std::nullopt;
  const auto sigma_field = (!opt.sigma_field.empty())
      ? std::optional<std::span<const double>>(frame.require_dfield(opt.sigma_field))
      : std::nullopt;
  const auto eps_field = (!opt.epsilon_field.empty())
      ? std::optional<std::span<const double>>(frame.require_dfield(opt.epsilon_field))
      : std::nullopt;

  OccluderData out;
  out.pos.reserve(occluders.idx.size());
  out.radius.reserve(occluders.idx.size());
  out.sigma.reserve(occluders.idx.size());
  out.epsilon.reserve(occluders.idx.size());
  for (const std::size_t i : occluders.idx) {
    out.pos.push_back(atom_vec3(frame, i));
    const double ri = rad_field ? (*rad_field)[i] : opt.occluder_radius;
    if (ri < 0.0) throw std::runtime_error("occluder radius must be >= 0");
    out.radius.push_back(ri);
    double sigma_i = sigma_field ? (*sigma_field)[i] : 0.0;
    if (!(sigma_i > 0.0)) {
      sigma_i = (opt.occluder_sigma > 0.0) ? opt.occluder_sigma : sigma_from_radius(ri);
    }
    if (sigma_i < 0.0) throw std::runtime_error("occluder sigma must be >= 0");
    out.sigma.push_back(sigma_i);
    const double eps_i = eps_field ? (*eps_field)[i] : opt.occluder_epsilon;
    if (eps_i < 0.0) throw std::runtime_error("occluder epsilon must be >= 0");
    out.epsilon.push_back(eps_i);
  }
  return out;
}

FrameInsertionStats evaluate_insertion_frame(const Frame& frame,
                                             const OccluderData& occ,
                                             const SelectionView* anchor_sel,
                                             const InsertionOptions& opt,
                                             std::size_t samples_plane,
                                             std::size_t samples_axis) {
  if (opt.n_bins == 0) throw std::runtime_error("insertion-profile evaluation requires n_bins >= 1");
  if (samples_plane == 0 || samples_axis == 0) {
    throw std::runtime_error("insertion-profile evaluation requires samples_plane >= 1 and samples_axis >= 1");
  }
  if (opt.probe_radius < 0.0) throw std::runtime_error("probe_radius must be >= 0");
  if (opt.beta <= 0.0) throw std::runtime_error("beta must be > 0");

  const SlabFrameTransform tf = make_slab_frame_transform(frame, anchor_sel, opt.slab);
  const std::size_t branches = opt.slab.halfcell_fold ? 2 : 1;
  const std::size_t per_bin = samples_plane * samples_plane * samples_axis * branches;

  FrameInsertionStats out;
  out.domain_length = tf.domain_length();
  out.samples_per_bin = per_bin;
  out.n_total.assign(opt.n_bins, 0);
  out.n_access.assign(opt.n_bins, 0);
  out.weight_sum.assign(opt.n_bins, 0.0);
  out.h.assign(opt.n_bins, 0.0);
  out.total_mean.assign(opt.n_bins, 0.0);
  out.cond_mean.assign(opt.n_bins, 0.0);

  const double probe_sigma_default = (opt.probe_sigma > 0.0) ? opt.probe_sigma : sigma_from_radius(opt.probe_radius);

  for (std::size_t b = 0; b < opt.n_bins; ++b) {
    std::size_t n_total = 0;
    std::size_t n_access = 0;
    double weight_sum = 0.0;
    for (std::size_t az = 0; az < samples_axis; ++az) {
      const double domain_frac = (static_cast<double>(b)
                                + (static_cast<double>(az) + 0.5) / static_cast<double>(samples_axis))
                               / static_cast<double>(opt.n_bins);
      for (std::size_t iu = 0; iu < samples_plane; ++iu) {
        const double fu = (static_cast<double>(iu) + 0.5) / static_cast<double>(samples_plane);
        for (std::size_t iv = 0; iv < samples_plane; ++iv) {
          const double fv = (static_cast<double>(iv) + 0.5) / static_cast<double>(samples_plane);
          for (std::size_t branch = 0; branch < branches; ++branch) {
            const auto lam = tf.lambda_for_domain_sample(domain_frac, fu, fv, static_cast<int>(branch));
            const auto cart = frame.box.from_lambda(lam[0], lam[1], lam[2]);
            const Vec3 probe{cart[0], cart[1], cart[2]};
            bool accessible = true;
            double u_soft = 0.0;
            for (std::size_t k = 0; k < occ.pos.size(); ++k) {
              const Vec3 dr = min_image_difference(frame.box, occ.pos[k], probe);
              const double r2 = norm2(dr);
              const double contact = opt.probe_radius + occ.radius[k];
              if (r2 < contact * contact) {
                accessible = false;
                break;
              }
              if (opt.energy_model != PairEnergyModel::None) {
                const double r = std::sqrt(r2);
                double sigma_eff = 0.0;
                if (probe_sigma_default > 0.0 || occ.sigma[k] > 0.0) {
                  sigma_eff = 0.5 * (probe_sigma_default + occ.sigma[k]);
                }
                if (!(sigma_eff > 0.0)) sigma_eff = sigma_fallback(contact);
                const double eps_eff = std::sqrt(std::max(0.0, opt.probe_epsilon) * std::max(0.0, occ.epsilon[k]));
                u_soft += pair_soft_energy(r, sigma_eff, eps_eff, opt.soft_cutoff, opt.energy_model);
              }
            }
            ++n_total;
            if (accessible) {
              ++n_access;
              weight_sum += std::exp(-opt.beta * u_soft);
            }
          }
        }
      }
    }
    out.n_total[b] = n_total;
    out.n_access[b] = n_access;
    out.weight_sum[b] = weight_sum;
    out.h[b] = (n_total > 0) ? (static_cast<double>(n_access) / static_cast<double>(n_total)) : 0.0;
    out.total_mean[b] = (n_total > 0) ? (weight_sum / static_cast<double>(n_total)) : 0.0;
    out.cond_mean[b] = (n_access > 0) ? (weight_sum / static_cast<double>(n_access)) : 0.0;
  }
  return out;
}

struct ReservoirFrameValues {
  bool valid_h = false;
  bool valid_total = false;
  bool valid_cond = false;
  double h = 1.0;
  double total = 1.0;
  double cond = 1.0;
};

ReservoirFrameValues frame_reservoir_values(const FrameInsertionStats& s,
                                            const std::vector<std::size_t>& res_bins,
                                            const ReservoirReferenceOptions& opt) {
  ReservoirFrameValues out;
  if (!opt.enabled) {
    out.valid_h = true;
    out.valid_total = true;
    out.valid_cond = true;
    out.h = 1.0;
    out.total = 1.0;
    out.cond = 1.0;
    return out;
  }
  out.h = average_over_bins(s.h, res_bins, nullptr);
  out.valid_h = std::isfinite(out.h) && out.h > 0.0;
  out.total = average_over_bins(s.total_mean, res_bins, nullptr);
  out.valid_total = std::isfinite(out.total) && out.total > 0.0;
  std::vector<unsigned char> cond_valid(s.cond_mean.size(), 0);
  for (std::size_t b = 0; b < s.cond_mean.size(); ++b) cond_valid[b] = (s.n_access[b] > 0) ? 1u : 0u;
  out.cond = average_over_bins(s.cond_mean, res_bins, &cond_valid);
  out.valid_cond = std::isfinite(out.cond) && out.cond > 0.0;
  return out;
}

struct StatsAccumulator {
  std::vector<double> mean_sum;
  std::vector<double> mean_sumsq;
  std::vector<std::size_t> mean_n;

  void resize(std::size_t n) {
    mean_sum.assign(n, 0.0);
    mean_sumsq.assign(n, 0.0);
    mean_n.assign(n, 0);
  }

  void add(std::size_t i, double x) {
    if (!std::isfinite(x)) return;
    mean_sum[i] += x;
    mean_sumsq[i] += x * x;
    ++mean_n[i];
  }

  double mean(std::size_t i) const {
    return (mean_n[i] > 0) ? (mean_sum[i] / static_cast<double>(mean_n[i])) : 0.0;
  }

  double sem(std::size_t i) const {
    return sem_from_sums(mean_sum[i], mean_sumsq[i], mean_n[i]);
  }
};

class AccessibilityProfileMeasure final : public IMeasure {
public:
  AccessibilityProfileMeasure(std::string type_name,
                              std::string instance_name,
                              std::string output_path,
                              SelectionView occluders,
                              std::optional<SelectionView> anchor_sel,
                              InsertionOptions opt)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)) {
    occluder_name_owned_ = std::string(occluders.name);
    occluders_ = SelectionView{occluder_name_owned_, occluders.idx};
    if (anchor_sel.has_value()) {
      have_anchor_ = true;
      anchor_name_owned_ = std::string(anchor_sel->name);
      anchor_sel_ = SelectionView{anchor_name_owned_, anchor_sel->idx};
    }
    if (occluders_.idx.empty()) throw std::runtime_error(type_name_ + ": occluder selection is empty");
    if (opt_.n_bins == 0) throw std::runtime_error(type_name_ + ": n_bins must be >= 1");
    h_frame_.resize(opt_.n_bins);
    h_rel_frame_.resize(opt_.n_bins);
    minus_log_h_rel_frame_.resize(opt_.n_bins);
    conv_abs_frame_.resize(opt_.n_bins);
    conv_rel_frame_.resize(opt_.n_bins);
    pooled_total_trials_.assign(opt_.n_bins, 0);
    pooled_access_trials_.assign(opt_.n_bins, 0);
    reservoir_bins_ = opt_.reservoir.enabled ? reservoir_bins(opt_.n_bins, opt_.reservoir) : std::vector<std::size_t>{};
  }

  std::string type() const override { return type_name_; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(occluders_.name);
    md.n_selected = occluders_.idx.size();
    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.x_axis = "bin";
    od.x_unit = "fractional_coordinate";
    od.columns = {
        "bin", "coord_lo_frac", "coord_hi_frac", "coord_center_frac", "coord_center_mean",
        "h_mean", "h_sem", "h_res", "h_rel", "h_rel_sem",
        "minus_log_h_rel", "minus_log_h_rel_sem",
        "conv_absdiff", "conv_absdiff_sem", "conv_reldiff", "conv_reldiff_sem",
        "n_trials_total", "n_access_total", "n_samples_per_frame", "n_frames"};
    md.outputs.push_back(std::move(od));
    md.params["axis"] = axis1d_name(opt_.slab.axis);
    md.params["probe_radius"] = dstr(opt_.probe_radius);
    md.params["samples_plane"] = std::to_string(opt_.samples_plane);
    md.params["samples_axis"] = std::to_string(opt_.samples_axis);
    md.params["convergence_refine_factor"] = std::to_string(opt_.convergence_refine_factor);
    md.params["slab_align_recenter"] = opt_.slab.slab_align_recenter ? "true" : "false";
    md.params["halfcell_fold"] = opt_.slab.halfcell_fold ? "true" : "false";
    md.params["target_center_frac"] = dstr(opt_.slab.target_center_frac);
    md.params["reservoir_left_lo_frac"] = dstr(opt_.reservoir.left_lo_frac);
    md.params["reservoir_left_hi_frac"] = dstr(opt_.reservoir.left_hi_frac);
    md.params["reservoir_use_right"] = opt_.reservoir.use_right_window ? "true" : "false";
    md.params["reservoir_right_lo_frac"] = dstr(opt_.reservoir.right_lo_frac);
    md.params["reservoir_right_hi_frac"] = dstr(opt_.reservoir.right_hi_frac);
    md.params["reservoir_normalize"] = opt_.reservoir.enabled ? "true" : "false";
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const SelectionView* anchor_ptr = have_anchor_ ? &anchor_sel_ : nullptr;
    const OccluderData occ = build_occluder_data(frame, occluders_, opt_);
    const FrameInsertionStats coarse = evaluate_insertion_frame(frame, occ, anchor_ptr, opt_, opt_.samples_plane, opt_.samples_axis);
    std::optional<FrameInsertionStats> refined;
    if (opt_.convergence_refine_factor > 1) {
      const std::size_t plane_ref = opt_.samples_plane * opt_.convergence_refine_factor;
      const std::size_t axis_ref = opt_.samples_axis * opt_.convergence_refine_factor;
      refined = evaluate_insertion_frame(frame, occ, anchor_ptr, opt_, plane_ref, axis_ref);
    }
    mean_domain_length_ += coarse.domain_length;
    samples_per_bin_ = coarse.samples_per_bin;
    for (std::size_t b = 0; b < opt_.n_bins; ++b) {
      pooled_total_trials_[b] += coarse.n_total[b];
      pooled_access_trials_[b] += coarse.n_access[b];
      h_frame_.add(b, coarse.h[b]);
      if (refined.has_value()) {
        const double absdiff = std::abs(refined->h[b] - coarse.h[b]);
        conv_abs_frame_.add(b, absdiff);
        const double denom = std::max(1e-12, std::abs(refined->h[b]));
        conv_rel_frame_.add(b, absdiff / denom);
      }
    }

    const ReservoirFrameValues refvals = frame_reservoir_values(coarse, reservoir_bins_, opt_.reservoir);
    for (std::size_t b = 0; b < opt_.n_bins; ++b) {
      if (refvals.valid_h && coarse.h[b] > 0.0) {
        const double h_rel = coarse.h[b] / refvals.h;
        h_rel_frame_.add(b, h_rel);
        minus_log_h_rel_frame_.add(b, -std::log(h_rel));
      }
    }
    ++n_frames_;
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: slab accessibility profile from probe insertion geometry\n";
      ofs << "# occluders: " << occluders_.name << " (n=" << occluders_.idx.size() << ")\n";
      if (have_anchor_) ofs << "# anchor: " << anchor_sel_.name << "\n";
      ofs << "# axis: " << axis1d_name(opt_.slab.axis)
          << ", probe_radius: " << std::setprecision(17) << opt_.probe_radius
          << ", samples_plane: " << opt_.samples_plane
          << ", samples_axis: " << opt_.samples_axis
          << ", convergence_refine_factor: " << opt_.convergence_refine_factor
          << ", slab_align_recenter: " << (opt_.slab.slab_align_recenter ? "true" : "false")
          << ", halfcell_fold: " << (opt_.slab.halfcell_fold ? "true" : "false") << "\n";
      if (opt_.reservoir.enabled) {
        ofs << "# reservoir windows: left=[" << std::setprecision(17) << opt_.reservoir.left_lo_frac << ", "
            << opt_.reservoir.left_hi_frac << "]";
        if (opt_.reservoir.use_right_window) {
          ofs << ", right=[" << opt_.reservoir.right_lo_frac << ", " << opt_.reservoir.right_hi_frac << "]";
        }
        ofs << "\n";
      } else {
        ofs << "# reservoir normalization: disabled (reported *_rel and mu_* use the absolute z-resolved means; *_res = 1)\n";
      }
      ofs << "# columns: bin coord_lo_frac coord_hi_frac coord_center_frac coord_center_mean h_mean h_sem h_res h_rel h_rel_sem minus_log_h_rel minus_log_h_rel_sem conv_absdiff conv_absdiff_sem conv_reldiff conv_reldiff_sem n_trials_total n_access_total n_samples_per_frame n_frames\n";
      const double mean_D = (n_frames_ > 0) ? (mean_domain_length_ / static_cast<double>(n_frames_)) : 0.0;
      std::vector<double> h_pooled(opt_.n_bins, 0.0);
      for (std::size_t b = 0; b < opt_.n_bins; ++b) {
        h_pooled[b] = (pooled_total_trials_[b] > 0)
            ? (static_cast<double>(pooled_access_trials_[b]) / static_cast<double>(pooled_total_trials_[b]))
            : 0.0;
      }
      const double h_res = opt_.reservoir.enabled ? average_over_bins(h_pooled, reservoir_bins_, nullptr) : 1.0;
      for (std::size_t b = 0; b < opt_.n_bins; ++b) {
        const double lo_frac = static_cast<double>(b) / static_cast<double>(opt_.n_bins);
        const double hi_frac = static_cast<double>(b + 1) / static_cast<double>(opt_.n_bins);
        const double center_frac = 0.5 * (lo_frac + hi_frac);
        const double h_mean = h_pooled[b];
        const double h_sem = h_frame_.sem(b);
        const double h_rel = (std::isfinite(h_res) && h_res > 0.0) ? (h_mean / h_res) : std::numeric_limits<double>::quiet_NaN();
        const double minus_log_h_rel = (std::isfinite(h_rel) && h_rel > 0.0) ? (-std::log(h_rel)) : std::numeric_limits<double>::infinity();
        ofs << b << ' '
            << std::setprecision(17) << lo_frac << ' '
            << std::setprecision(17) << hi_frac << ' '
            << std::setprecision(17) << center_frac << ' '
            << std::setprecision(17) << center_frac * mean_D << ' '
            << std::setprecision(17) << h_mean << ' '
            << std::setprecision(17) << h_sem << ' '
            << std::setprecision(17) << h_res << ' '
            << std::setprecision(17) << h_rel << ' '
            << std::setprecision(17) << h_rel_frame_.sem(b) << ' '
            << std::setprecision(17) << minus_log_h_rel << ' '
            << std::setprecision(17) << minus_log_h_rel_frame_.sem(b) << ' '
            << std::setprecision(17) << conv_abs_frame_.mean(b) << ' '
            << std::setprecision(17) << conv_abs_frame_.sem(b) << ' '
            << std::setprecision(17) << conv_rel_frame_.mean(b) << ' '
            << std::setprecision(17) << conv_rel_frame_.sem(b) << ' '
            << pooled_total_trials_[b] << ' '
            << pooled_access_trials_[b] << ' '
            << samples_per_bin_ << ' '
            << n_frames_ << '\n';
      }
    });
  }

  void finalize() override {
    if (!opt_.range.dry_run) flush_partial();
  }

private:
  std::string type_name_;
  std::string instance_name_;
  std::string output_path_;
  std::string occluder_name_owned_;
  std::string anchor_name_owned_;
  SelectionView occluders_{};
  SelectionView anchor_sel_{};
  bool have_anchor_ = false;
  InsertionOptions opt_;
  bool started_ = false;
  StatsAccumulator h_frame_;
  StatsAccumulator h_rel_frame_;
  StatsAccumulator minus_log_h_rel_frame_;
  StatsAccumulator conv_abs_frame_;
  StatsAccumulator conv_rel_frame_;
  std::vector<std::size_t> reservoir_bins_;
  std::vector<std::size_t> pooled_total_trials_;
  std::vector<std::size_t> pooled_access_trials_;
  double mean_domain_length_ = 0.0;
  std::size_t samples_per_bin_ = 0;
  std::size_t n_frames_ = 0;
};

class WidomProfileMeasure final : public IMeasure {
public:
  enum class Mode { Total, Conditional };

  WidomProfileMeasure(std::string type_name,
                      std::string instance_name,
                      std::string output_path,
                      SelectionView occluders,
                      std::optional<SelectionView> anchor_sel,
                      InsertionOptions opt,
                      Mode mode)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        mode_(mode) {
    occluder_name_owned_ = std::string(occluders.name);
    occluders_ = SelectionView{occluder_name_owned_, occluders.idx};
    if (anchor_sel.has_value()) {
      have_anchor_ = true;
      anchor_name_owned_ = std::string(anchor_sel->name);
      anchor_sel_ = SelectionView{anchor_name_owned_, anchor_sel->idx};
    }
    if (occluders_.idx.empty()) throw std::runtime_error(type_name_ + ": occluder selection is empty");
    if (opt_.n_bins == 0) throw std::runtime_error(type_name_ + ": n_bins must be >= 1");
    reservoir_bins_ = opt_.reservoir.enabled ? reservoir_bins(opt_.n_bins, opt_.reservoir) : std::vector<std::size_t>{};

    h_frame_.resize(opt_.n_bins);
    total_frame_.resize(opt_.n_bins);
    total_rel_frame_.resize(opt_.n_bins);
    mu_total_frame_.resize(opt_.n_bins);
    cond_frame_.resize(opt_.n_bins);
    cond_rel_frame_.resize(opt_.n_bins);
    mu_cond_frame_.resize(opt_.n_bins);
    h_rel_frame_.resize(opt_.n_bins);
    mu_hard_frame_.resize(opt_.n_bins);
    audit_frame_.resize(opt_.n_bins);
    conv_abs_frame_.resize(opt_.n_bins);
    conv_rel_frame_.resize(opt_.n_bins);

    pooled_total_trials_.assign(opt_.n_bins, 0);
    pooled_access_trials_.assign(opt_.n_bins, 0);
    pooled_weight_sum_.assign(opt_.n_bins, 0.0);
  }

  std::string type() const override { return type_name_; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(occluders_.name);
    md.n_selected = occluders_.idx.size();
    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.x_axis = "bin";
    od.x_unit = "fractional_coordinate";
    if (mode_ == Mode::Conditional) {
      od.columns = {
          "bin", "coord_lo_frac", "coord_hi_frac", "coord_center_frac", "coord_center_mean",
          "h_mean", "h_sem", "h_res", "h_rel", "h_rel_sem", "mu_hard_rel", "mu_hard_rel_sem",
          "cond_boltz_mean", "cond_boltz_sem", "cond_boltz_res", "cond_rel", "cond_rel_sem", "mu_cond_rel", "mu_cond_rel_sem",
          "full_boltz_mean", "full_boltz_sem", "full_boltz_res", "full_rel", "full_rel_sem", "mu_full_rel", "mu_full_rel_sem",
          "identity_residual", "identity_residual_sem",
          "conv_absdiff", "conv_absdiff_sem", "conv_reldiff", "conv_reldiff_sem",
          "n_trials_total", "n_access_total", "n_samples_per_frame", "n_frames"};
    } else {
      od.columns = {
          "bin", "coord_lo_frac", "coord_hi_frac", "coord_center_frac", "coord_center_mean",
          "h_mean", "h_sem",
          "boltz_mean", "boltz_sem", "boltz_res", "boltz_rel", "boltz_rel_sem", "mu_rel", "mu_rel_sem",
          "conv_absdiff", "conv_absdiff_sem", "conv_reldiff", "conv_reldiff_sem",
          "n_trials_total", "n_access_total", "n_samples_per_frame", "n_frames"};
    }
    md.outputs.push_back(std::move(od));
    md.params["axis"] = axis1d_name(opt_.slab.axis);
    md.params["probe_radius"] = dstr(opt_.probe_radius);
    md.params["energy_model"] = pair_energy_model_name(opt_.energy_model);
    md.params["beta"] = dstr(opt_.beta);
    md.params["samples_plane"] = std::to_string(opt_.samples_plane);
    md.params["samples_axis"] = std::to_string(opt_.samples_axis);
    md.params["convergence_refine_factor"] = std::to_string(opt_.convergence_refine_factor);
    md.params["slab_align_recenter"] = opt_.slab.slab_align_recenter ? "true" : "false";
    md.params["halfcell_fold"] = opt_.slab.halfcell_fold ? "true" : "false";
    md.params["target_center_frac"] = dstr(opt_.slab.target_center_frac);
    md.params["reservoir_left_lo_frac"] = dstr(opt_.reservoir.left_lo_frac);
    md.params["reservoir_left_hi_frac"] = dstr(opt_.reservoir.left_hi_frac);
    md.params["reservoir_use_right"] = opt_.reservoir.use_right_window ? "true" : "false";
    md.params["reservoir_right_lo_frac"] = dstr(opt_.reservoir.right_lo_frac);
    md.params["reservoir_right_hi_frac"] = dstr(opt_.reservoir.right_hi_frac);
    md.params["reservoir_normalize"] = opt_.reservoir.enabled ? "true" : "false";
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const SelectionView* anchor_ptr = have_anchor_ ? &anchor_sel_ : nullptr;
    const OccluderData occ = build_occluder_data(frame, occluders_, opt_);
    const FrameInsertionStats coarse = evaluate_insertion_frame(frame, occ, anchor_ptr, opt_, opt_.samples_plane, opt_.samples_axis);
    std::optional<FrameInsertionStats> refined;
    if (opt_.convergence_refine_factor > 1) {
      const std::size_t plane_ref = opt_.samples_plane * opt_.convergence_refine_factor;
      const std::size_t axis_ref = opt_.samples_axis * opt_.convergence_refine_factor;
      refined = evaluate_insertion_frame(frame, occ, anchor_ptr, opt_, plane_ref, axis_ref);
    }
    mean_domain_length_ += coarse.domain_length;
    samples_per_bin_ = coarse.samples_per_bin;
    for (std::size_t b = 0; b < opt_.n_bins; ++b) {
      pooled_total_trials_[b] += coarse.n_total[b];
      pooled_access_trials_[b] += coarse.n_access[b];
      pooled_weight_sum_[b] += coarse.weight_sum[b];
      h_frame_.add(b, coarse.h[b]);
      total_frame_.add(b, coarse.total_mean[b]);
      if (coarse.n_access[b] > 0) cond_frame_.add(b, coarse.cond_mean[b]);
      if (refined.has_value()) {
        const double absdiff = std::abs(refined->h[b] - coarse.h[b]);
        conv_abs_frame_.add(b, absdiff);
        const double denom = std::max(1e-12, std::abs(refined->h[b]));
        conv_rel_frame_.add(b, absdiff / denom);
      }
    }

    const ReservoirFrameValues refvals = frame_reservoir_values(coarse, reservoir_bins_, opt_.reservoir);
    for (std::size_t b = 0; b < opt_.n_bins; ++b) {
      if (refvals.valid_h && coarse.h[b] > 0.0) {
        const double h_rel = coarse.h[b] / refvals.h;
        h_rel_frame_.add(b, h_rel);
        mu_hard_frame_.add(b, -std::log(h_rel) / opt_.beta);
      }
      if (refvals.valid_total && coarse.total_mean[b] > 0.0) {
        const double total_rel = coarse.total_mean[b] / refvals.total;
        total_rel_frame_.add(b, total_rel);
        mu_total_frame_.add(b, -std::log(total_rel) / opt_.beta);
      }
      if (refvals.valid_cond && coarse.n_access[b] > 0 && coarse.cond_mean[b] > 0.0) {
        const double cond_rel = coarse.cond_mean[b] / refvals.cond;
        cond_rel_frame_.add(b, cond_rel);
        mu_cond_frame_.add(b, -std::log(cond_rel) / opt_.beta);
      }
      if (coarse.h[b] > 0.0 && coarse.total_mean[b] > 0.0 && coarse.n_access[b] > 0 && coarse.cond_mean[b] > 0.0) {
        audit_frame_.add(b, -std::log(coarse.total_mean[b] / (coarse.h[b] * coarse.cond_mean[b])) / opt_.beta);
      }
    }
    ++n_frames_;
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: z-resolved charge-off Widom insertion profile\n";
      ofs << "# mode: " << (mode_ == Mode::Conditional ? "conditional" : "total")
          << ", energy_model: " << pair_energy_model_name(opt_.energy_model)
          << ", beta: " << std::setprecision(17) << opt_.beta << "\n";
      ofs << "# occluders: " << occluders_.name << " (n=" << occluders_.idx.size() << ")\n";
      if (have_anchor_) ofs << "# anchor: " << anchor_sel_.name << "\n";
      if (opt_.reservoir.enabled) {
        ofs << "# reservoir windows: left=[" << std::setprecision(17) << opt_.reservoir.left_lo_frac << ", "
            << opt_.reservoir.left_hi_frac << "]";
        if (opt_.reservoir.use_right_window) {
          ofs << ", right=[" << opt_.reservoir.right_lo_frac << ", " << opt_.reservoir.right_hi_frac << "]";
        }
        ofs << "\n";
      } else {
        ofs << "# reservoir normalization: disabled (reported *_rel and mu_* use the absolute z-resolved means; *_res = 1)\n";
      }
      if (mode_ == Mode::Conditional) {
        ofs << "# columns: bin coord_lo_frac coord_hi_frac coord_center_frac coord_center_mean h_mean h_sem h_res h_rel h_rel_sem mu_hard_rel mu_hard_rel_sem cond_boltz_mean cond_boltz_sem cond_boltz_res cond_rel cond_rel_sem mu_cond_rel mu_cond_rel_sem full_boltz_mean full_boltz_sem full_boltz_res full_rel full_rel_sem mu_full_rel mu_full_rel_sem identity_residual identity_residual_sem conv_absdiff conv_absdiff_sem conv_reldiff conv_reldiff_sem n_trials_total n_access_total n_samples_per_frame n_frames\n";
      } else {
        ofs << "# columns: bin coord_lo_frac coord_hi_frac coord_center_frac coord_center_mean h_mean h_sem boltz_mean boltz_sem boltz_res boltz_rel boltz_rel_sem mu_rel mu_rel_sem conv_absdiff conv_absdiff_sem conv_reldiff conv_reldiff_sem n_trials_total n_access_total n_samples_per_frame n_frames\n";
      }
      const double mean_D = (n_frames_ > 0) ? (mean_domain_length_ / static_cast<double>(n_frames_)) : 0.0;

      std::vector<double> h_pooled(opt_.n_bins, 0.0);
      std::vector<double> total_pooled(opt_.n_bins, 0.0);
      std::vector<double> cond_pooled(opt_.n_bins, 0.0);
      for (std::size_t b = 0; b < opt_.n_bins; ++b) {
        h_pooled[b] = (pooled_total_trials_[b] > 0)
            ? (static_cast<double>(pooled_access_trials_[b]) / static_cast<double>(pooled_total_trials_[b]))
            : 0.0;
        total_pooled[b] = (pooled_total_trials_[b] > 0)
            ? (pooled_weight_sum_[b] / static_cast<double>(pooled_total_trials_[b]))
            : 0.0;
        cond_pooled[b] = (pooled_access_trials_[b] > 0)
            ? (pooled_weight_sum_[b] / static_cast<double>(pooled_access_trials_[b]))
            : 0.0;
      }
      std::vector<unsigned char> cond_valid(opt_.n_bins, 0);
      for (std::size_t b = 0; b < opt_.n_bins; ++b) cond_valid[b] = (pooled_access_trials_[b] > 0) ? 1u : 0u;
      const double h_res = opt_.reservoir.enabled ? average_over_bins(h_pooled, reservoir_bins_, nullptr) : 1.0;
      const double total_res = opt_.reservoir.enabled ? average_over_bins(total_pooled, reservoir_bins_, nullptr) : 1.0;
      const double cond_res = opt_.reservoir.enabled ? average_over_bins(cond_pooled, reservoir_bins_, &cond_valid) : 1.0;

      for (std::size_t b = 0; b < opt_.n_bins; ++b) {
        const double lo_frac = static_cast<double>(b) / static_cast<double>(opt_.n_bins);
        const double hi_frac = static_cast<double>(b + 1) / static_cast<double>(opt_.n_bins);
        const double center_frac = 0.5 * (lo_frac + hi_frac);
        const double h_mean = h_pooled[b];
        const double total_mean = total_pooled[b];
        const double cond_mean = cond_pooled[b];
        if (mode_ == Mode::Conditional) {
          const double h_rel = (std::isfinite(h_res) && h_res > 0.0) ? (h_mean / h_res) : std::numeric_limits<double>::quiet_NaN();
          const double total_rel = (std::isfinite(total_res) && total_res > 0.0) ? (total_mean / total_res) : std::numeric_limits<double>::quiet_NaN();
          const double cond_rel = (std::isfinite(cond_res) && cond_res > 0.0) ? (cond_mean / cond_res) : std::numeric_limits<double>::quiet_NaN();
          const double mu_hard = (std::isfinite(h_rel) && h_rel > 0.0) ? (-std::log(h_rel) / opt_.beta) : std::numeric_limits<double>::infinity();
          const double mu_total = (std::isfinite(total_rel) && total_rel > 0.0) ? (-std::log(total_rel) / opt_.beta) : std::numeric_limits<double>::infinity();
          const double mu_cond = (std::isfinite(cond_rel) && cond_rel > 0.0) ? (-std::log(cond_rel) / opt_.beta) : std::numeric_limits<double>::infinity();
          const double identity = (h_mean > 0.0 && total_mean > 0.0 && cond_mean > 0.0)
              ? (-std::log(total_mean / (h_mean * cond_mean)) / opt_.beta)
              : std::numeric_limits<double>::quiet_NaN();
          ofs << b << ' '
              << std::setprecision(17) << lo_frac << ' '
              << std::setprecision(17) << hi_frac << ' '
              << std::setprecision(17) << center_frac << ' '
              << std::setprecision(17) << center_frac * mean_D << ' '
              << std::setprecision(17) << h_mean << ' '
              << std::setprecision(17) << h_frame_.sem(b) << ' '
              << std::setprecision(17) << h_res << ' '
              << std::setprecision(17) << h_rel << ' '
              << std::setprecision(17) << h_rel_frame_.sem(b) << ' '
              << std::setprecision(17) << mu_hard << ' '
              << std::setprecision(17) << mu_hard_frame_.sem(b) << ' '
              << std::setprecision(17) << cond_mean << ' '
              << std::setprecision(17) << cond_frame_.sem(b) << ' '
              << std::setprecision(17) << cond_res << ' '
              << std::setprecision(17) << cond_rel << ' '
              << std::setprecision(17) << cond_rel_frame_.sem(b) << ' '
              << std::setprecision(17) << mu_cond << ' '
              << std::setprecision(17) << mu_cond_frame_.sem(b) << ' '
              << std::setprecision(17) << total_mean << ' '
              << std::setprecision(17) << total_frame_.sem(b) << ' '
              << std::setprecision(17) << total_res << ' '
              << std::setprecision(17) << total_rel << ' '
              << std::setprecision(17) << total_rel_frame_.sem(b) << ' '
              << std::setprecision(17) << mu_total << ' '
              << std::setprecision(17) << mu_total_frame_.sem(b) << ' '
              << std::setprecision(17) << identity << ' '
              << std::setprecision(17) << audit_frame_.sem(b) << ' '
              << std::setprecision(17) << conv_abs_frame_.mean(b) << ' '
              << std::setprecision(17) << conv_abs_frame_.sem(b) << ' '
              << std::setprecision(17) << conv_rel_frame_.mean(b) << ' '
              << std::setprecision(17) << conv_rel_frame_.sem(b) << ' '
              << pooled_total_trials_[b] << ' '
              << pooled_access_trials_[b] << ' '
              << samples_per_bin_ << ' '
              << n_frames_ << '\n';
        } else {
          const double total_rel = (std::isfinite(total_res) && total_res > 0.0) ? (total_mean / total_res) : std::numeric_limits<double>::quiet_NaN();
          const double mu_total = (std::isfinite(total_rel) && total_rel > 0.0) ? (-std::log(total_rel) / opt_.beta) : std::numeric_limits<double>::infinity();
          ofs << b << ' '
              << std::setprecision(17) << lo_frac << ' '
              << std::setprecision(17) << hi_frac << ' '
              << std::setprecision(17) << center_frac << ' '
              << std::setprecision(17) << center_frac * mean_D << ' '
              << std::setprecision(17) << h_mean << ' '
              << std::setprecision(17) << h_frame_.sem(b) << ' '
              << std::setprecision(17) << total_mean << ' '
              << std::setprecision(17) << total_frame_.sem(b) << ' '
              << std::setprecision(17) << total_res << ' '
              << std::setprecision(17) << total_rel << ' '
              << std::setprecision(17) << total_rel_frame_.sem(b) << ' '
              << std::setprecision(17) << mu_total << ' '
              << std::setprecision(17) << mu_total_frame_.sem(b) << ' '
              << std::setprecision(17) << conv_abs_frame_.mean(b) << ' '
              << std::setprecision(17) << conv_abs_frame_.sem(b) << ' '
              << std::setprecision(17) << conv_rel_frame_.mean(b) << ' '
              << std::setprecision(17) << conv_rel_frame_.sem(b) << ' '
              << pooled_total_trials_[b] << ' '
              << pooled_access_trials_[b] << ' '
              << samples_per_bin_ << ' '
              << n_frames_ << '\n';
        }
      }
    });
  }

  void finalize() override {
    if (!opt_.range.dry_run) flush_partial();
  }

private:
  std::string type_name_;
  std::string instance_name_;
  std::string output_path_;
  std::string occluder_name_owned_;
  std::string anchor_name_owned_;
  SelectionView occluders_{};
  SelectionView anchor_sel_{};
  bool have_anchor_ = false;
  InsertionOptions opt_;
  Mode mode_ = Mode::Total;
  bool started_ = false;
  std::vector<std::size_t> reservoir_bins_;
  StatsAccumulator h_frame_;
  StatsAccumulator total_frame_;
  StatsAccumulator total_rel_frame_;
  StatsAccumulator mu_total_frame_;
  StatsAccumulator cond_frame_;
  StatsAccumulator cond_rel_frame_;
  StatsAccumulator mu_cond_frame_;
  StatsAccumulator h_rel_frame_;
  StatsAccumulator mu_hard_frame_;
  StatsAccumulator audit_frame_;
  StatsAccumulator conv_abs_frame_;
  StatsAccumulator conv_rel_frame_;
  std::vector<std::size_t> pooled_total_trials_;
  std::vector<std::size_t> pooled_access_trials_;
  std::vector<double> pooled_weight_sum_;
  double mean_domain_length_ = 0.0;
  std::size_t samples_per_bin_ = 0;
  std::size_t n_frames_ = 0;
};

void append_insertion_caps(const IniConfig& cfg,
                           const std::string& section,
                           MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  if (cfg.has_key(section, "anchor_group")) caps.group_refs.push_back(cfg.get_string(section, "anchor_group"));
  if (cfg.has_key(section, "radius_field")) caps.requires_dfields.push_back(cfg.get_string(section, "radius_field"));
  if (cfg.has_key(section, "sigma_field")) caps.requires_dfields.push_back(cfg.get_string(section, "sigma_field"));
  const auto model = parse_pair_energy_model(cfg.get_string(section, "energy_model", std::optional<std::string>("none")));
  if (model != PairEnergyModel::None && cfg.has_key(section, "epsilon_field")) {
    caps.requires_dfields.push_back(cfg.get_string(section, "epsilon_field"));
  }
}

std::optional<SelectionView> build_anchor_selection(const IniConfig& cfg,
                                                    const std::string& section,
                                                    SelectionProvider& sp,
                                                    const Frame& frame0,
                                                    const SelectionView& fallback_sel,
                                                    bool needed) {
  return build_optional_anchor_selection(cfg, section, sp, frame0, fallback_sel, needed, "slab");
}

inline ReservoirReferenceOptions reservoir_options_from_config(const IniConfig& cfg,
                                                               const std::string& section,
                                                               bool halfcell_fold) {
  ReservoirReferenceOptions opt;
  opt.enabled = cfg.get_bool(section, "reservoir_normalize", std::optional<bool>(true));
  opt.left_lo_frac = cfg.get_double(section, "reservoir_lo_frac", std::optional<double>(0.0));
  opt.left_hi_frac = cfg.get_double(section, "reservoir_hi_frac", std::optional<double>(0.1));
  if (cfg.has_key(section, "reservoir_left_lo_frac")) {
    opt.left_lo_frac = cfg.get_double(section, "reservoir_left_lo_frac");
  }
  if (cfg.has_key(section, "reservoir_left_hi_frac")) {
    opt.left_hi_frac = cfg.get_double(section, "reservoir_left_hi_frac");
  }
  validate_window(opt.left_lo_frac, opt.left_hi_frac, section + ": reservoir left window");
  if (cfg.has_key(section, "reservoir_use_right")) {
    opt.use_right_window = cfg.get_bool(section, "reservoir_use_right");
  } else {
    opt.use_right_window = !halfcell_fold;
  }
  opt.right_lo_frac = cfg.get_double(section, "reservoir_right_lo_frac", std::optional<double>(0.9));
  opt.right_hi_frac = cfg.get_double(section, "reservoir_right_hi_frac", std::optional<double>(1.0));
  if (opt.use_right_window) validate_window(opt.right_lo_frac, opt.right_hi_frac, section + ": reservoir right window");
  return opt;
}

InsertionOptions insertion_options_from_config(const IniConfig& cfg,
                                               const std::string& section,
                                               bool require_beta) {
  InsertionOptions opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.slab.axis = parse_axis1d(cfg.get_string(section, "axis", std::optional<std::string>("z")));
  opt.slab.slab_align_recenter = cfg.get_bool(section, "slab_align_recenter", std::optional<bool>(false));
  opt.slab.halfcell_fold = cfg.get_bool(section, "halfcell_fold", std::optional<bool>(false));
  opt.slab.target_center_frac = cfg.get_double(section, "target_center_frac", std::optional<double>(0.5));
  opt.reservoir = reservoir_options_from_config(cfg, section, opt.slab.halfcell_fold);
  opt.n_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_bins", std::optional<std::int64_t>(100)));
  opt.samples_plane = static_cast<std::size_t>(cfg.get_int64(section, "samples_plane", std::optional<std::int64_t>(16)));
  opt.samples_axis = static_cast<std::size_t>(cfg.get_int64(section, "samples_axis", std::optional<std::int64_t>(2)));
  opt.convergence_refine_factor = static_cast<std::size_t>(cfg.get_int64(section, "convergence_refine_factor", std::optional<std::int64_t>(2)));
  opt.probe_radius = cfg.get_double(section, "probe_radius");
  opt.radius_field = cfg.get_string(section, "radius_field", std::optional<std::string>(""));
  opt.occluder_radius = cfg.get_double(section, "occluder_radius", std::optional<double>(0.0));
  opt.sigma_field = cfg.get_string(section, "sigma_field", std::optional<std::string>(""));
  opt.occluder_sigma = cfg.get_double(section, "occluder_sigma", std::optional<double>(0.0));
  opt.probe_sigma = cfg.get_double(section, "probe_sigma", std::optional<double>(0.0));
  opt.energy_model = parse_pair_energy_model(cfg.get_string(section, "energy_model", std::optional<std::string>("none")));
  opt.epsilon_field = cfg.get_string(section, "epsilon_field", std::optional<std::string>(""));
  opt.occluder_epsilon = cfg.get_double(section, "occluder_epsilon", std::optional<double>(1.0));
  opt.probe_epsilon = cfg.get_double(section, "probe_epsilon", std::optional<double>(1.0));
  opt.soft_cutoff = cfg.get_double(section, "soft_cutoff", std::optional<double>(0.0));
  if (require_beta || cfg.has_key(section, "beta")) {
    opt.beta = cfg.get_double(section, "beta", std::optional<double>(1.0));
  }
  if (opt.range.frame_start < 0) throw std::runtime_error(section + ": frame_start must be >= 0");
  if (opt.range.frame_end >= 0 && opt.range.frame_end < opt.range.frame_start) {
    throw std::runtime_error(section + ": frame_end must be -1 or >= frame_start");
  }
  if (opt.n_bins == 0) throw std::runtime_error(section + ": n_bins must be >= 1");
  if (opt.samples_plane == 0 || opt.samples_axis == 0) {
    throw std::runtime_error(section + ": samples_plane and samples_axis must be >= 1");
  }
  if (opt.convergence_refine_factor == 0) throw std::runtime_error(section + ": convergence_refine_factor must be >= 1");
  if (opt.probe_radius < 0.0) throw std::runtime_error(section + ": probe_radius must be >= 0");
  if (opt.occluder_radius < 0.0) throw std::runtime_error(section + ": occluder_radius must be >= 0");
  if (opt.probe_sigma < 0.0) throw std::runtime_error(section + ": probe_sigma must be >= 0");
  if (opt.occluder_sigma < 0.0) throw std::runtime_error(section + ": occluder_sigma must be >= 0");
  if (opt.probe_epsilon < 0.0) throw std::runtime_error(section + ": probe_epsilon must be >= 0");
  if (opt.occluder_epsilon < 0.0) throw std::runtime_error(section + ": occluder_epsilon must be >= 0");
  if (opt.soft_cutoff < 0.0) throw std::runtime_error(section + ": soft_cutoff must be >= 0");
  if (opt.beta <= 0.0) throw std::runtime_error(section + ": beta must be > 0");
  return opt;
}

MeasureCapabilities accessibility_caps(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& instance,
                                       const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_insertion_caps(cfg, section, caps);
  return caps;
}

std::unique_ptr<IMeasure> accessibility_create(const IniConfig& cfg,
                                               const std::string& section,
                                               const std::string& instance,
                                               const MeasureBuildEnv& env,
                                               const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) {
    throw std::runtime_error("accessibility_profile factory: missing first_frame or SelectionProvider");
  }
  const Frame& frame0 = *env.first_frame;
  const std::string type_name = cfg.get_string(section, "type", std::optional<std::string>("accessibility_profile"));
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView occluders = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);
  InsertionOptions opt = insertion_options_from_config(cfg, section, false);
  opt.range.dry_run = env.dry_run;
  auto anchor = build_anchor_selection(cfg, section, *env.selection_provider, frame0, occluders,
                                       opt.slab.slab_align_recenter || opt.slab.halfcell_fold);
  const fs::path out = measure_ext::resolve_measure_output_path(cfg, section, env, "output", "accessibility_profile.dat");
  return std::make_unique<AccessibilityProfileMeasure>(type_name, instance, out.string(), occluders, anchor, std::move(opt));
}

MeasureCapabilities widom_caps(const IniConfig& cfg,
                               const std::string& section,
                               const std::string& instance,
                               const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_insertion_caps(cfg, section, caps);
  return caps;
}

std::unique_ptr<IMeasure> widom_create(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& instance,
                                       const MeasureBuildEnv& env,
                                       const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) {
    throw std::runtime_error("widom_profile factory: missing first_frame or SelectionProvider");
  }
  const Frame& frame0 = *env.first_frame;
  const std::string type_name = cfg.get_string(section, "type", std::optional<std::string>("widom_profile"));
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView occluders = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);
  InsertionOptions opt = insertion_options_from_config(cfg, section, true);
  opt.range.dry_run = env.dry_run;
  auto anchor = build_anchor_selection(cfg, section, *env.selection_provider, frame0, occluders,
                                       opt.slab.slab_align_recenter || opt.slab.halfcell_fold);
  const fs::path out = measure_ext::resolve_measure_output_path(cfg, section, env, "output", type_name + std::string(".dat"));
  const auto mode = (type_name == "conditional_widom_profile") ? WidomProfileMeasure::Mode::Conditional
                                                                : WidomProfileMeasure::Mode::Total;
  return std::make_unique<WidomProfileMeasure>(type_name, instance, out.string(), occluders, anchor, std::move(opt), mode);
}

static MeasureRegistrar g_reg_accessibility("accessibility_profile", &accessibility_caps, &accessibility_create);
static MeasureRegistrar g_reg_probe_accessibility("probe_accessibility_profile", &accessibility_caps, &accessibility_create);
static MeasureRegistrar g_reg_widom("widom_profile", &widom_caps, &widom_create);
static MeasureRegistrar g_reg_conditional_widom("conditional_widom_profile", &widom_caps, &widom_create);

} // namespace
} // namespace pilots
