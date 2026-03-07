#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <memory>
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

using measure_ext::integer_like_field_to_i64;
using measure_ext::get_static_combined_view;
using measure_ext::get_static_group_view;
using measure_ext::resolve_path;
using measure_ext::same_index_set;

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

inline double mic_delta(double d, double L) {
  if (!(L > 0.0)) return d;
  return d - L * std::round(d / L);
}

struct Vec3 { double x = 0.0, y = 0.0, z = 0.0; };
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3{a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator*(double s, const Vec3& a) { return Vec3{s * a.x, s * a.y, s * a.z}; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline double norm(const Vec3& a) { return std::sqrt(dot(a,a)); }
inline Vec3 mic_vec(const Vec3& a, const Vec3& b, double lx, double ly, double lz) {
  return Vec3{mic_delta(a.x - b.x, lx), mic_delta(a.y - b.y, ly), mic_delta(a.z - b.z, lz)};
}

struct DSU {
  explicit DSU(std::size_t n) : p(n), r(n, 0) { for (std::size_t i = 0; i < n; ++i) p[i] = i; }
  std::size_t find(std::size_t x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
  void unite(std::size_t a, std::size_t b) { a = find(a); b = find(b); if (a == b) return; if (r[a] < r[b]) std::swap(a,b); p[b] = a; if (r[a] == r[b]) ++r[a]; }
  std::vector<std::size_t> p;
  std::vector<int> r;
};

inline long double fact_int(int n) {
  if (n < 0) return 0.0L;
  long double v = 1.0L;
  for (int i = 2; i <= n; ++i) v *= static_cast<long double>(i);
  return v;
}

inline long double double_fact_odd(int n) {
  if (n <= 0) return 1.0L;
  long double v = 1.0L;
  for (int k = n; k >= 1; k -= 2) v *= static_cast<long double>(k);
  return v;
}

inline double assoc_legendre(int l, int m, double x) {
  m = std::abs(m);
  x = std::clamp(x, -1.0, 1.0);
  if (m > l) return 0.0;
  const double sx = std::sqrt(std::max(0.0, 1.0 - x*x));
  double pmm = ((m % 2) ? -1.0 : 1.0) * static_cast<double>(double_fact_odd(2*m - 1)) * std::pow(sx, m);
  if (l == m) return pmm;
  double pm1m = x * static_cast<double>(2*m + 1) * pmm;
  if (l == m + 1) return pm1m;
  double pkm2 = pmm;
  double pkm1 = pm1m;
  double pk = 0.0;
  for (int k = m + 2; k <= l; ++k) {
    pk = ((2.0 * k - 1.0) * x * pkm1 - (k + m - 1.0) * pkm2) / static_cast<double>(k - m);
    pkm2 = pkm1;
    pkm1 = pk;
  }
  return pk;
}

inline std::complex<double> spherical_harmonic(int l, int m, const Vec3& u) {
  const double r = norm(u);
  if (r <= 0.0) return {0.0, 0.0};
  const double cth = std::clamp(u.z / r, -1.0, 1.0);
  const double phi = std::atan2(u.y, u.x);
  if (m < 0) {
    const int mp = -m;
    const auto y = spherical_harmonic(l, mp, u);
    return ((mp % 2) ? -1.0 : 1.0) * std::conj(y);
  }
  const long double pref = std::sqrt(((2.0L * l + 1.0L) / (4.0L * std::numbers::pi_v<long double>))
                                     * fact_int(l - m) / fact_int(l + m));
  const double plm = assoc_legendre(l, m, cth);
  const std::complex<double> eimphi(std::cos(static_cast<double>(m) * phi), std::sin(static_cast<double>(m) * phi));
  return static_cast<double>(pref) * plm * eimphi;
}

inline long double wigner_3j(int j1, int j2, int j3, int m1, int m2, int m3) {
  if (m1 + m2 + m3 != 0) return 0.0L;
  if (std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(m3) > j3) return 0.0L;
  if (j3 > j1 + j2 || j3 < std::abs(j1 - j2)) return 0.0L;
  const int t1 = j2 - j3 - m1;
  const int t2 = j1 - j3 + m2;
  const int tmin = std::max(0, std::max(t1, t2));
  const int tmax = std::min({j1 + j2 - j3, j1 - m1, j2 + m2});
  if (tmin > tmax) return 0.0L;

  const long double delta = std::sqrt(fact_int(j1 + j2 - j3)
                                      * fact_int(j1 - j2 + j3)
                                      * fact_int(-j1 + j2 + j3)
                                      / fact_int(j1 + j2 + j3 + 1));
  const long double normf = std::sqrt(fact_int(j1 + m1) * fact_int(j1 - m1)
                                      * fact_int(j2 + m2) * fact_int(j2 - m2)
                                      * fact_int(j3 + m3) * fact_int(j3 - m3));
  long double sum = 0.0L;
  for (int t = tmin; t <= tmax; ++t) {
    const long double den = fact_int(t)
                           * fact_int(j1 + j2 - j3 - t)
                           * fact_int(j1 - m1 - t)
                           * fact_int(j2 + m2 - t)
                           * fact_int(j3 - j2 + m1 + t)
                           * fact_int(j3 - j1 - m2 + t);
    const long double term = (((t % 2) == 0) ? 1.0L : -1.0L) / den;
    sum += term;
  }
  const long double phase = (((j1 - j2 - m3) % 2) == 0) ? 1.0L : -1.0L;
  return phase * delta * normf * sum;
}

struct FrameDescriptors {
  std::vector<double> q4;
  std::vector<double> q6;
  std::vector<double> w6;
  std::vector<double> w6_hat;
  std::vector<int> cn;
};

struct DescriptorOptions {
  double neighbor_cutoff = 0.0;
};

FrameDescriptors compute_descriptors(const Frame& frame,
                                     const SelectionView& sel,
                                     const SelectionView& nbr_sel,
                                     const DescriptorOptions& opt) {
  if (!(opt.neighbor_cutoff > 0.0)) throw std::runtime_error("local structure descriptors: neighbor_cutoff must be > 0");
  const auto xu = frame.require_dfield("xu");
  const auto yu = frame.require_dfield("yu");
  const auto zu = frame.require_dfield("zu");

  std::vector<Vec3> pos(frame.natoms);
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < frame.natoms; ++i) pos[i] = Vec3{xu[i], yu[i], zu[i]};

  FrameDescriptors out;
  out.q4.assign(sel.idx.size(), 0.0);
  out.q6.assign(sel.idx.size(), 0.0);
  out.w6.assign(sel.idx.size(), 0.0);
  out.w6_hat.assign(sel.idx.size(), 0.0);
  out.cn.assign(sel.idx.size(), 0);

  const bool same = same_index_set(sel, nbr_sel);
  const double rc = opt.neighbor_cutoff;
  const double rc2 = rc * rc;

#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
  for (std::size_t p = 0; p < sel.idx.size(); ++p) {
    const std::size_t i = sel.idx[p];
    std::vector<std::complex<double>> q4m(9, {0.0, 0.0});
    std::vector<std::complex<double>> q6m(13, {0.0, 0.0});
    int cn = 0;
    for (const std::size_t j : nbr_sel.idx) {
      if (same && j == i) continue;
      const Vec3 rij = mic_vec(pos[j], pos[i], frame.box.lx(), frame.box.ly(), frame.box.lz());
      const double r2 = dot(rij, rij);
      if (r2 <= 0.0 || r2 > rc2) continue;
      ++cn;
      for (int m = -4; m <= 4; ++m) q4m[m + 4] += spherical_harmonic(4, m, rij);
      for (int m = -6; m <= 6; ++m) q6m[m + 6] += spherical_harmonic(6, m, rij);
    }
    out.cn[p] = cn;
    if (cn == 0) continue;
    const double inv = 1.0 / static_cast<double>(cn);
    double s4 = 0.0, s6 = 0.0;
    for (auto& z : q4m) { z *= inv; s4 += std::norm(z); }
    for (auto& z : q6m) { z *= inv; s6 += std::norm(z); }
    out.q4[p] = std::sqrt((4.0 * std::numbers::pi) / 9.0 * s4);
    out.q6[p] = std::sqrt((4.0 * std::numbers::pi) / 13.0 * s6);

    std::complex<double> w6 = {0.0, 0.0};
    for (int m1 = -6; m1 <= 6; ++m1) {
      for (int m2 = -6; m2 <= 6; ++m2) {
        const int m3 = -(m1 + m2);
        if (m3 < -6 || m3 > 6) continue;
        const long double c3j = wigner_3j(6, 6, 6, m1, m2, m3);
        if (std::abs(static_cast<double>(c3j)) < 1e-16) continue;
        w6 += static_cast<double>(c3j) * q6m[m1 + 6] * q6m[m2 + 6] * q6m[m3 + 6];
      }
    }
    out.w6[p] = w6.real();
    const double denom = std::pow(s6, 1.5);
    out.w6_hat[p] = (denom > 0.0) ? (w6.real() / denom) : 0.0;
  }
  return out;
}

struct SummaryRecord {
  std::size_t frame_index = 0;
  std::int64_t timestep = 0;
  double a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0;
  std::size_t n = 0;
};

struct LFSRecord {
  std::size_t frame_index = 0;
  std::int64_t timestep = 0;
  double frac = 0.0;
  double n_lfs = 0.0;
  double mean_cluster = 0.0;
  double max_cluster = 0.0;
  double n_clusters = 0.0;
  double mean_q6 = 0.0;
  double mean_w6hat = 0.0;
};

struct SoftnessRecord {
  std::size_t frame_index = 0;
  std::int64_t timestep = 0;
  double mean_soft = 0.0;
  double std_soft = 0.0;
  double mean_q6 = 0.0;
  double mean_w6hat = 0.0;
  double mean_cn = 0.0;
  double mean_free = 0.0;
  std::size_t n = 0;
};

struct Q4PerRecord {
  std::size_t frame_index = 0;
  std::int64_t timestep = 0;
  std::size_t local_pos = 0;
  double q4 = 0.0;
  double q6 = 0.0;
  double w6 = 0.0;
  double w6_hat = 0.0;
  int cn = 0;
};

struct SoftnessPerRecord {
  std::size_t frame_index = 0;
  std::int64_t timestep = 0;
  std::size_t local_pos = 0;
  double softness = 0.0;
  double q6 = 0.0;
  double w6_hat = 0.0;
  int cn = 0;
  double free_volume = 0.0;
};

class Q4Q6W6Measure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    DescriptorOptions desc;
    bool write_per_particle = false;
  };

  Q4Q6W6Measure(std::string instance_name,
                std::string summary_path,
                std::string per_particle_path,
                SelectionView sel,
                SelectionView nbr_sel,
                std::vector<std::int64_t> ids,
                Options opt)
      : instance_name_(std::move(instance_name)),
        summary_path_(std::move(summary_path)),
        per_particle_path_(std::move(per_particle_path)),
        ids_(std::move(ids)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    nbr_name_owned_ = std::string(nbr_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    nbr_sel_ = SelectionView{nbr_name_owned_, nbr_sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error("q4_q6_w6: selection is empty");
    if (nbr_sel_.idx.empty()) throw std::runtime_error("q4_q6_w6: neighbor selection is empty");
  }

  std::string type() const override { return "q4_q6_w6"; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();
    output::OutputFileDescriptor s;
    s.path = summary_path_;
    s.format = "text";
    s.columns = {"frame", "timestep", "mean_q4", "mean_q6", "mean_w6", "mean_w6_hat", "mean_cn", "n_particles"};
    md.outputs.push_back(std::move(s));
    if (opt_.write_per_particle) {
      output::OutputFileDescriptor p;
      p.path = per_particle_path_;
      p.format = "text";
      p.columns = {"frame", "timestep", "atom", "id", "q4", "q6", "w6", "w6_hat", "cn"};
      md.outputs.push_back(std::move(p));
    }
    return md;
  }

  void on_start(const Frame& first_frame) override { (void)first_frame; }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto d = compute_descriptors(frame, sel_, nbr_sel_, opt_.desc);
    SummaryRecord r;
    r.frame_index = frame_index;
    r.timestep = frame.timestep;
    r.n = sel_.idx.size();
    for (std::size_t i = 0; i < sel_.idx.size(); ++i) {
      r.a += d.q4[i];
      r.b += d.q6[i];
      r.c += d.w6[i];
      r.d += d.w6_hat[i];
      r.e += static_cast<double>(d.cn[i]);
      if (opt_.write_per_particle) per_.push_back(Q4PerRecord{frame_index, frame.timestep, i, d.q4[i], d.q6[i], d.w6[i], d.w6_hat[i], d.cn[i]});
    }
    const double inv = 1.0 / static_cast<double>(sel_.idx.size());
    r.a *= inv; r.b *= inv; r.c *= inv; r.d *= inv; r.e *= inv;
    rec_.push_back(r);
  }

  void flush_partial() override {}
  void finalize() override {
    if (opt_.range.dry_run) return;
    util::atomic_write_text(summary_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: Steinhardt local bond-orientational order descriptors\n";
      ofs << "# columns: frame timestep mean_q4 mean_q6 mean_w6 mean_w6_hat mean_cn n_particles\n";
      for (const auto& r : rec_) {
        ofs << r.frame_index << ' ' << r.timestep << ' '
            << std::setprecision(17) << r.a << ' '
            << std::setprecision(17) << r.b << ' '
            << std::setprecision(17) << r.c << ' '
            << std::setprecision(17) << r.d << ' '
            << std::setprecision(17) << r.e << ' '
            << r.n << '\n';
      }
    });
    if (opt_.write_per_particle) {
      util::atomic_write_text(per_particle_path_, [&](std::ostream& ofs) {
        ofs << "# PILOTS: per-particle q4 q6 w6 w6_hat coordination number\n";
        ofs << "# columns: frame timestep atom id q4 q6 w6 w6_hat cn\n";
        for (const auto& v : per_) {
          const std::size_t atom = sel_.idx[v.local_pos];
          const std::int64_t id = (atom < ids_.size()) ? ids_[atom] : static_cast<std::int64_t>(atom);
          ofs << v.frame_index << ' ' << v.timestep << ' ' << atom << ' ' << id << ' '
              << std::setprecision(17) << v.q4 << ' '
              << std::setprecision(17) << v.q6 << ' '
              << std::setprecision(17) << v.w6 << ' '
              << std::setprecision(17) << v.w6_hat << ' '
              << v.cn << '\n';
        }
      });
    }
  }

private:
  std::string instance_name_;
  std::string summary_path_;
  std::string per_particle_path_;
  std::string sel_name_owned_;
  std::string nbr_name_owned_;
  SelectionView sel_;
  SelectionView nbr_sel_;
  std::vector<std::int64_t> ids_;
  Options opt_;
  std::vector<SummaryRecord> rec_;
  std::vector<Q4PerRecord> per_;
};

class LFSMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    DescriptorOptions desc;
    double q4_min = -std::numeric_limits<double>::infinity();
    double q4_max =  std::numeric_limits<double>::infinity();
    double q6_min = -std::numeric_limits<double>::infinity();
    double q6_max =  std::numeric_limits<double>::infinity();
    double w6h_min = -std::numeric_limits<double>::infinity();
    double w6h_max =  std::numeric_limits<double>::infinity();
    int cn_min = std::numeric_limits<int>::min();
    int cn_max = std::numeric_limits<int>::max();
    double cluster_cutoff = 0.0;
  };

  LFSMeasure(std::string instance_name,
             std::string output_path,
             SelectionView sel,
             SelectionView nbr_sel,
             Options opt)
      : instance_name_(std::move(instance_name)), output_path_(std::move(output_path)), opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    nbr_name_owned_ = std::string(nbr_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    nbr_sel_ = SelectionView{nbr_name_owned_, nbr_sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error("locally_favored_structure: selection is empty");
  }

  std::string type() const override { return "locally_favored_structure"; }
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
    od.columns = {"frame", "timestep", "fraction_lfs", "n_lfs", "mean_lfs_cluster_size", "max_lfs_cluster_size", "n_lfs_clusters", "mean_q6_lfs", "mean_w6hat_lfs"};
    md.outputs.push_back(std::move(od));
    return md;
  }

  void on_start(const Frame& first_frame) override { (void)first_frame; }
  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto d = compute_descriptors(frame, sel_, nbr_sel_, opt_.desc);
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    std::vector<Vec3> pos(frame.natoms);
    for (std::size_t i = 0; i < frame.natoms; ++i) pos[i] = Vec3{xu[i], yu[i], zu[i]};

    std::vector<unsigned char> is_lfs(sel_.idx.size(), 0);
    std::vector<std::size_t> lfs_pos;
    double q6_sum = 0.0, w6h_sum = 0.0;
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const bool ok = d.q4[p] >= opt_.q4_min && d.q4[p] <= opt_.q4_max
                   && d.q6[p] >= opt_.q6_min && d.q6[p] <= opt_.q6_max
                   && d.w6_hat[p] >= opt_.w6h_min && d.w6_hat[p] <= opt_.w6h_max
                   && d.cn[p] >= opt_.cn_min && d.cn[p] <= opt_.cn_max;
      if (ok) {
        is_lfs[p] = 1;
        lfs_pos.push_back(p);
        q6_sum += d.q6[p];
        w6h_sum += d.w6_hat[p];
      }
    }
    double mean_cluster = 0.0;
    std::size_t max_cluster = 0;
    std::size_t n_clusters = 0;
    if (opt_.cluster_cutoff > 0.0 && !lfs_pos.empty()) {
      DSU dsu(lfs_pos.size());
      for (std::size_t a = 0; a < lfs_pos.size(); ++a) {
        for (std::size_t b = a + 1; b < lfs_pos.size(); ++b) {
          const std::size_t ia = sel_.idx[lfs_pos[a]];
          const std::size_t ib = sel_.idx[lfs_pos[b]];
          const Vec3 ra = pos[ia], rb = pos[ib];
          if (norm(mic_vec(ra, rb, frame.box.lx(), frame.box.ly(), frame.box.lz())) <= opt_.cluster_cutoff) dsu.unite(a,b);
        }
      }
      std::unordered_map<std::size_t, std::size_t> sizes;
      for (std::size_t a = 0; a < lfs_pos.size(); ++a) ++sizes[dsu.find(a)];
      n_clusters = sizes.size();
      for (const auto& kv : sizes) {
        mean_cluster += static_cast<double>(kv.second);
        max_cluster = std::max(max_cluster, kv.second);
      }
      mean_cluster /= static_cast<double>(n_clusters);
    }
    LFSRecord r;
    r.frame_index = frame_index;
    r.timestep = frame.timestep;
    r.frac = static_cast<double>(lfs_pos.size()) / static_cast<double>(sel_.idx.size());
    r.n_lfs = static_cast<double>(lfs_pos.size());
    r.mean_cluster = mean_cluster;
    r.max_cluster = static_cast<double>(max_cluster);
    r.n_clusters = static_cast<double>(n_clusters);
    r.mean_q6 = lfs_pos.empty() ? 0.0 : q6_sum / static_cast<double>(lfs_pos.size());
    r.mean_w6hat = lfs_pos.empty() ? 0.0 : w6h_sum / static_cast<double>(lfs_pos.size());
    rec_.push_back(r);
  }
  void flush_partial() override {}
  void finalize() override {
    if (opt_.range.dry_run) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: generic locally favored structure classifier from structural descriptors\n";
      ofs << "# columns: frame timestep fraction_lfs n_lfs mean_lfs_cluster_size max_lfs_cluster_size n_lfs_clusters mean_q6_lfs mean_w6hat_lfs\n";
      for (const auto& r : rec_) {
        ofs << r.frame_index << ' ' << r.timestep << ' '
            << std::setprecision(17) << r.frac << ' '
            << std::setprecision(17) << r.n_lfs << ' '
            << std::setprecision(17) << r.mean_cluster << ' '
            << std::setprecision(17) << r.max_cluster << ' '
            << std::setprecision(17) << r.n_clusters << ' '
            << std::setprecision(17) << r.mean_q6 << ' '
            << std::setprecision(17) << r.mean_w6hat << '\n';
      }
    });
  }

private:
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::string nbr_name_owned_;
  SelectionView sel_;
  SelectionView nbr_sel_;
  Options opt_;
  std::vector<LFSRecord> rec_;
};

struct FieldStatsRow {
  std::size_t frame_index = 0;
  std::int64_t timestep = 0;
  double mean = 0.0, stddev = 0.0, vmin = 0.0, vmax = 0.0;
  std::size_t n = 0;
};

class FieldStatsMeasure final : public IMeasure {
public:
  enum class Kind { Voronoi, FreeVolume };
  struct Options {
    RangeOptions range;
    Kind kind = Kind::Voronoi;
    std::string voronoi_field;
    std::string free_volume_field;
    std::string excluded_volume_field;
    std::string radius_field;
    double particle_radius = std::numeric_limits<double>::quiet_NaN();
  };

  FieldStatsMeasure(std::string type_name,
                    std::string instance_name,
                    std::string output_path,
                    SelectionView sel,
                    Options opt)
      : type_name_(std::move(type_name)), instance_name_(std::move(instance_name)), output_path_(std::move(output_path)), opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error(type_name_ + ": selection is empty");
  }

  std::string type() const override { return type_name_; }
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
    od.columns = {"frame", "timestep", "mean", "stddev", "min", "max", "n_particles"};
    md.outputs.push_back(std::move(od));
    return md;
  }
  void on_start(const Frame& first_frame) override { (void)first_frame; }
  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    std::vector<double> values(sel_.idx.size(), 0.0);
    if (opt_.kind == Kind::Voronoi) {
      if (opt_.voronoi_field.empty()) throw std::runtime_error("voronoi_volume: voronoi_volume_field must be provided");
      const auto v = frame.require_dfield(opt_.voronoi_field);
      for (std::size_t p = 0; p < sel_.idx.size(); ++p) values[p] = v[sel_.idx[p]];
    } else {
      if (!opt_.free_volume_field.empty()) {
        const auto v = frame.require_dfield(opt_.free_volume_field);
        for (std::size_t p = 0; p < sel_.idx.size(); ++p) values[p] = v[sel_.idx[p]];
      } else {
        if (opt_.voronoi_field.empty()) throw std::runtime_error("free_volume: need free_volume_field or voronoi_volume_field");
        const auto vv = frame.require_dfield(opt_.voronoi_field);
        std::optional<std::span<const double>> excl_field;
        std::optional<std::span<const double>> rad_field;
        if (!opt_.excluded_volume_field.empty()) excl_field = frame.require_dfield(opt_.excluded_volume_field);
        if (!opt_.radius_field.empty()) rad_field = frame.require_dfield(opt_.radius_field);
        for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
          const std::size_t i = sel_.idx[p];
          double excl = 0.0;
          if (excl_field) excl = (*excl_field)[i];
          else if (rad_field) { const double r = (*rad_field)[i]; excl = (4.0 / 3.0) * std::numbers::pi * r * r * r; }
          else if (std::isfinite(opt_.particle_radius)) { const double r = opt_.particle_radius; excl = (4.0 / 3.0) * std::numbers::pi * r * r * r; }
          values[p] = vv[i] - excl;
        }
      }
    }
    double mean = 0.0, mean2 = 0.0;
    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    for (const double x : values) {
      mean += x; mean2 += x * x; vmin = std::min(vmin, x); vmax = std::max(vmax, x);
    }
    const double inv = 1.0 / static_cast<double>(values.size());
    mean *= inv; mean2 *= inv;
    const double stddev = std::sqrt(std::max(0.0, mean2 - mean * mean));
    rec_.push_back({frame_index, frame.timestep, mean, stddev, vmin, vmax, values.size()});
  }
  void flush_partial() override {}
  void finalize() override {
    if (opt_.range.dry_run) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: field statistics for " << type_name_ << "\n";
      ofs << "# columns: frame timestep mean stddev min max n_particles\n";
      for (const auto& r : rec_) {
        ofs << r.frame_index << ' ' << r.timestep << ' '
            << std::setprecision(17) << r.mean << ' '
            << std::setprecision(17) << r.stddev << ' '
            << std::setprecision(17) << r.vmin << ' '
            << std::setprecision(17) << r.vmax << ' '
            << r.n << '\n';
      }
    });
  }
private:
  std::string type_name_;
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  SelectionView sel_;
  Options opt_;
  std::vector<FieldStatsRow> rec_;
};

class SoftnessProxyMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    DescriptorOptions desc;
    std::string voronoi_field;
    std::string free_volume_field;
    std::string excluded_volume_field;
    std::string radius_field;
    double particle_radius = std::numeric_limits<double>::quiet_NaN();
    double c0 = 0.0, c_q6 = 1.0, c_w6h = 0.0, c_cn = 0.0, c_free = 0.0;
    bool write_per_particle = false;
  };

  SoftnessProxyMeasure(std::string instance_name,
                       std::string summary_path,
                       std::string per_particle_path,
                       SelectionView sel,
                       SelectionView nbr_sel,
                       std::vector<std::int64_t> ids,
                       Options opt)
      : instance_name_(std::move(instance_name)), summary_path_(std::move(summary_path)), per_particle_path_(std::move(per_particle_path)), ids_(std::move(ids)), opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    nbr_name_owned_ = std::string(nbr_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    nbr_sel_ = SelectionView{nbr_name_owned_, nbr_sel.idx};
  }

  std::string type() const override { return "softness_proxy"; }
  std::string instance_name() const override { return instance_name_; }
  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();
    output::OutputFileDescriptor s;
    s.path = summary_path_;
    s.format = "text";
    s.columns = {"frame", "timestep", "mean_softness", "std_softness", "mean_q6", "mean_w6hat", "mean_cn", "mean_free_volume", "n_particles"};
    md.outputs.push_back(std::move(s));
    if (opt_.write_per_particle) {
      output::OutputFileDescriptor p;
      p.path = per_particle_path_;
      p.format = "text";
      p.columns = {"frame", "timestep", "atom", "id", "softness", "q6", "w6_hat", "cn", "free_volume"};
      md.outputs.push_back(std::move(p));
    }
    return md;
  }
  void on_start(const Frame& first_frame) override { (void)first_frame; }
  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto d = compute_descriptors(frame, sel_, nbr_sel_, opt_.desc);
    std::vector<double> fv(sel_.idx.size(), 0.0);
    if (!opt_.free_volume_field.empty()) {
      const auto v = frame.require_dfield(opt_.free_volume_field);
      for (std::size_t p = 0; p < sel_.idx.size(); ++p) fv[p] = v[sel_.idx[p]];
    } else if (!opt_.voronoi_field.empty()) {
      const auto vv = frame.require_dfield(opt_.voronoi_field);
      std::optional<std::span<const double>> excl_field;
      std::optional<std::span<const double>> rad_field;
      if (!opt_.excluded_volume_field.empty()) excl_field = frame.require_dfield(opt_.excluded_volume_field);
      if (!opt_.radius_field.empty()) rad_field = frame.require_dfield(opt_.radius_field);
      for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
        const std::size_t i = sel_.idx[p];
        double excl = 0.0;
        if (excl_field) excl = (*excl_field)[i];
        else if (rad_field) { const double r = (*rad_field)[i]; excl = (4.0 / 3.0) * std::numbers::pi * r * r * r; }
        else if (std::isfinite(opt_.particle_radius)) { const double r = opt_.particle_radius; excl = (4.0 / 3.0) * std::numbers::pi * r * r * r; }
        fv[p] = vv[i] - excl;
      }
    }
    double ms = 0.0, ms2 = 0.0, mq6 = 0.0, mw = 0.0, mcn = 0.0, mfv = 0.0;
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const double s = opt_.c0 + opt_.c_q6 * d.q6[p] + opt_.c_w6h * d.w6_hat[p] + opt_.c_cn * static_cast<double>(d.cn[p]) + opt_.c_free * fv[p];
      ms += s; ms2 += s*s; mq6 += d.q6[p]; mw += d.w6_hat[p]; mcn += static_cast<double>(d.cn[p]); mfv += fv[p];
      if (opt_.write_per_particle) per_.push_back(SoftnessPerRecord{frame_index, frame.timestep, p, s, d.q6[p], d.w6_hat[p], d.cn[p], fv[p]});
    }
    const double inv = 1.0 / static_cast<double>(sel_.idx.size());
    SoftnessRecord r;
    r.frame_index = frame_index;
    r.timestep = frame.timestep;
    r.mean_soft = ms * inv;
    r.std_soft = std::sqrt(std::max(0.0, ms2 * inv - (ms * inv) * (ms * inv)));
    r.mean_q6 = mq6 * inv;
    r.mean_w6hat = mw * inv;
    r.mean_cn = mcn * inv;
    r.mean_free = mfv * inv;
    r.n = sel_.idx.size();
    rec_.push_back(r);
  }
  void flush_partial() override {}
  void finalize() override {
    if (opt_.range.dry_run) return;
    util::atomic_write_text(summary_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: linear softness proxy from structural descriptors\n";
      ofs << "# columns: frame timestep mean_softness std_softness mean_q6 mean_w6hat mean_cn mean_free_volume n_particles\n";
      for (const auto& r : rec_) {
        ofs << r.frame_index << ' ' << r.timestep << ' '
            << std::setprecision(17) << r.mean_soft << ' '
            << std::setprecision(17) << r.std_soft << ' '
            << std::setprecision(17) << r.mean_q6 << ' '
            << std::setprecision(17) << r.mean_w6hat << ' '
            << std::setprecision(17) << r.mean_cn << ' '
            << std::setprecision(17) << r.mean_free << ' '
            << r.n << '\n';
      }
    });
    if (opt_.write_per_particle) {
      util::atomic_write_text(per_particle_path_, [&](std::ostream& ofs) {
        ofs << "# PILOTS: per-particle linear softness proxy\n";
        ofs << "# columns: frame timestep atom id softness q6 w6_hat cn free_volume\n";
        for (const auto& v : per_) {
          const std::size_t atom = sel_.idx[v.local_pos];
          const std::int64_t id = (atom < ids_.size()) ? ids_[atom] : static_cast<std::int64_t>(atom);
          ofs << v.frame_index << ' ' << v.timestep << ' ' << atom << ' ' << id << ' '
              << std::setprecision(17) << v.softness << ' '
              << std::setprecision(17) << v.q6 << ' '
              << std::setprecision(17) << v.w6_hat << ' '
              << v.cn << ' '
              << std::setprecision(17) << v.free_volume << '\n';
        }
      });
    }
  }

private:
  std::string instance_name_;
  std::string summary_path_;
  std::string per_particle_path_;
  std::string sel_name_owned_;
  std::string nbr_name_owned_;
  SelectionView sel_;
  SelectionView nbr_sel_;
  std::vector<std::int64_t> ids_;
  Options opt_;
  std::vector<SoftnessRecord> rec_;
  std::vector<SoftnessPerRecord> per_;
};

void append_local_caps(const IniConfig& cfg,
                       const std::string& section,
                       MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.requires_identity_consistent = true;
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  if (cfg.has_key(section, "neighbor_group")) caps.group_refs.push_back(cfg.get_string(section, "neighbor_group"));
}

std::unique_ptr<IMeasure> q4_q6_w6_create(const IniConfig& cfg,
                                          const std::string& section,
                                          const std::string& instance,
                                          const MeasureBuildEnv& env,
                                          const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("q4_q6_w6 factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "q4_q6_w6");
  SelectionView nbr_sel = cfg.has_key(section, "neighbor_group")
      ? get_static_group_view(*env.selection_provider, frame0, cfg.get_string(section, "neighbor_group"), "q4_q6_w6")
      : sel;
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path sum = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("q4_q6_w6.dat"))).lexically_normal();
  const fs::path per = (output_dir / cfg.get_string(section, "output_per_particle", std::optional<std::string>("q4_q6_w6_per_particle.dat"))).lexically_normal();
  const auto ids = integer_like_field_to_i64(frame0, cfg.get_string(section, "id_field", std::optional<std::string>("id")), true);

  Q4Q6W6Measure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.desc.neighbor_cutoff = cfg.get_double(section, "neighbor_cutoff");
  opt.write_per_particle = cfg.get_bool(section, "write_per_particle", std::optional<bool>(false));
  return std::make_unique<Q4Q6W6Measure>(instance, sum.string(), per.string(), sel, nbr_sel, ids, opt);
}

MeasureCapabilities q4_q6_w6_caps(const IniConfig& cfg,
                                  const std::string& section,
                                  const std::string& instance,
                                  const MeasureBuildEnv& env) {
  (void)instance; (void)env; MeasureCapabilities caps; append_local_caps(cfg, section, caps); caps.requires_intfields.push_back(cfg.get_string(section, "id_field", std::optional<std::string>("id"))); return caps; }

std::unique_ptr<IMeasure> lfs_create(const IniConfig& cfg,
                                     const std::string& section,
                                     const std::string& instance,
                                     const MeasureBuildEnv& env,
                                     const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("locally_favored_structure factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "locally_favored_structure");
  SelectionView nbr_sel = cfg.has_key(section, "neighbor_group")
      ? get_static_group_view(*env.selection_provider, frame0, cfg.get_string(section, "neighbor_group"), "locally_favored_structure")
      : sel;
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("locally_favored_structure.dat"))).lexically_normal();
  LFSMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.desc.neighbor_cutoff = cfg.get_double(section, "neighbor_cutoff");
  opt.q4_min = cfg.get_double(section, "q4_min", std::optional<double>(opt.q4_min));
  opt.q4_max = cfg.get_double(section, "q4_max", std::optional<double>(opt.q4_max));
  opt.q6_min = cfg.get_double(section, "q6_min", std::optional<double>(opt.q6_min));
  opt.q6_max = cfg.get_double(section, "q6_max", std::optional<double>(opt.q6_max));
  opt.w6h_min = cfg.get_double(section, "w6_hat_min", std::optional<double>(opt.w6h_min));
  opt.w6h_max = cfg.get_double(section, "w6_hat_max", std::optional<double>(opt.w6h_max));
  opt.cn_min = static_cast<int>(cfg.get_int64(section, "cn_min", std::optional<std::int64_t>(opt.cn_min)));
  opt.cn_max = static_cast<int>(cfg.get_int64(section, "cn_max", std::optional<std::int64_t>(opt.cn_max)));
  opt.cluster_cutoff = cfg.get_double(section, "cluster_cutoff", std::optional<double>(0.0));
  return std::make_unique<LFSMeasure>(instance, out.string(), sel, nbr_sel, opt);
}

MeasureCapabilities lfs_caps(const IniConfig& cfg,
                             const std::string& section,
                             const std::string& instance,
                             const MeasureBuildEnv& env) {
  (void)instance; (void)env; MeasureCapabilities caps; append_local_caps(cfg, section, caps); return caps; }

MeasureCapabilities voronoi_caps(const IniConfig& cfg,
                                 const std::string& section,
                                 const std::string& instance,
                                 const MeasureBuildEnv& env) {
  (void)instance; (void)env; MeasureCapabilities caps; caps.selection_policy = SelectionPolicy::RequireStatic; caps.requires_identity_consistent = true; caps.scale = ScaleCompatibility{true, true, true}; caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all"))); caps.requires_dfields.push_back(cfg.get_string(section, "voronoi_volume_field")); return caps; }

std::unique_ptr<IMeasure> voronoi_create(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env,
                                         const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("voronoi_volume factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "voronoi_volume");
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("voronoi_volume.dat"))).lexically_normal();
  FieldStatsMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.kind = FieldStatsMeasure::Kind::Voronoi;
  opt.voronoi_field = cfg.get_string(section, "voronoi_volume_field");
  return std::make_unique<FieldStatsMeasure>("voronoi_volume", instance, out.string(), sel, opt);
}

MeasureCapabilities free_volume_caps(const IniConfig& cfg,
                                     const std::string& section,
                                     const std::string& instance,
                                     const MeasureBuildEnv& env) {
  (void)instance; (void)env; MeasureCapabilities caps; caps.selection_policy = SelectionPolicy::RequireStatic; caps.requires_identity_consistent = true; caps.scale = ScaleCompatibility{true, true, true}; caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all"))); if (cfg.has_key(section, "free_volume_field")) caps.requires_dfields.push_back(cfg.get_string(section, "free_volume_field")); else caps.requires_dfields.push_back(cfg.get_string(section, "voronoi_volume_field")); if (cfg.has_key(section, "excluded_volume_field")) caps.requires_dfields.push_back(cfg.get_string(section, "excluded_volume_field")); if (cfg.has_key(section, "radius_field")) caps.requires_dfields.push_back(cfg.get_string(section, "radius_field")); return caps; }

std::unique_ptr<IMeasure> free_volume_create(const IniConfig& cfg,
                                             const std::string& section,
                                             const std::string& instance,
                                             const MeasureBuildEnv& env,
                                             const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("free_volume factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "free_volume");
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("free_volume.dat"))).lexically_normal();
  FieldStatsMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.kind = FieldStatsMeasure::Kind::FreeVolume;
  opt.voronoi_field = cfg.get_string(section, "voronoi_volume_field", std::optional<std::string>(""));
  opt.free_volume_field = cfg.get_string(section, "free_volume_field", std::optional<std::string>(""));
  opt.excluded_volume_field = cfg.get_string(section, "excluded_volume_field", std::optional<std::string>(""));
  opt.radius_field = cfg.get_string(section, "radius_field", std::optional<std::string>(""));
  opt.particle_radius = cfg.get_double(section, "particle_radius", std::optional<double>(std::numeric_limits<double>::quiet_NaN()));
  return std::make_unique<FieldStatsMeasure>("free_volume", instance, out.string(), sel, opt);
}

MeasureCapabilities softness_caps(const IniConfig& cfg,
                                  const std::string& section,
                                  const std::string& instance,
                                  const MeasureBuildEnv& env) {
  (void)instance; (void)env; MeasureCapabilities caps; append_local_caps(cfg, section, caps); caps.requires_intfields.push_back(cfg.get_string(section, "id_field", std::optional<std::string>("id"))); if (cfg.has_key(section, "free_volume_field")) caps.requires_dfields.push_back(cfg.get_string(section, "free_volume_field")); if (cfg.has_key(section, "voronoi_volume_field")) caps.requires_dfields.push_back(cfg.get_string(section, "voronoi_volume_field")); if (cfg.has_key(section, "excluded_volume_field")) caps.requires_dfields.push_back(cfg.get_string(section, "excluded_volume_field")); if (cfg.has_key(section, "radius_field")) caps.requires_dfields.push_back(cfg.get_string(section, "radius_field")); return caps; }

std::unique_ptr<IMeasure> softness_create(const IniConfig& cfg,
                                          const std::string& section,
                                          const std::string& instance,
                                          const MeasureBuildEnv& env,
                                          const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("softness_proxy factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "softness_proxy");
  SelectionView nbr_sel = cfg.has_key(section, "neighbor_group")
      ? get_static_group_view(*env.selection_provider, frame0, cfg.get_string(section, "neighbor_group"), "softness_proxy")
      : sel;
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path sum = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("softness_proxy.dat"))).lexically_normal();
  const fs::path per = (output_dir / cfg.get_string(section, "output_per_particle", std::optional<std::string>("softness_proxy_per_particle.dat"))).lexically_normal();
  const auto ids = integer_like_field_to_i64(frame0, cfg.get_string(section, "id_field", std::optional<std::string>("id")), true);

  SoftnessProxyMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.desc.neighbor_cutoff = cfg.get_double(section, "neighbor_cutoff");
  opt.voronoi_field = cfg.get_string(section, "voronoi_volume_field", std::optional<std::string>(""));
  opt.free_volume_field = cfg.get_string(section, "free_volume_field", std::optional<std::string>(""));
  opt.excluded_volume_field = cfg.get_string(section, "excluded_volume_field", std::optional<std::string>(""));
  opt.radius_field = cfg.get_string(section, "radius_field", std::optional<std::string>(""));
  opt.particle_radius = cfg.get_double(section, "particle_radius", std::optional<double>(std::numeric_limits<double>::quiet_NaN()));
  opt.c0 = cfg.get_double(section, "c0", std::optional<double>(0.0));
  opt.c_q6 = cfg.get_double(section, "c_q6", std::optional<double>(1.0));
  opt.c_w6h = cfg.get_double(section, "c_w6_hat", std::optional<double>(0.0));
  opt.c_cn = cfg.get_double(section, "c_cn", std::optional<double>(0.0));
  opt.c_free = cfg.get_double(section, "c_free", std::optional<double>(0.0));
  opt.write_per_particle = cfg.get_bool(section, "write_per_particle", std::optional<bool>(false));
  return std::make_unique<SoftnessProxyMeasure>(instance, sum.string(), per.string(), sel, nbr_sel, ids, opt);
}

static MeasureRegistrar g_reg_q4_q6_w6("q4_q6_w6", &q4_q6_w6_caps, &q4_q6_w6_create);
static MeasureRegistrar g_reg_lfs("locally_favored_structure", &lfs_caps, &lfs_create);
static MeasureRegistrar g_reg_voronoi("voronoi_volume", &voronoi_caps, &voronoi_create);
static MeasureRegistrar g_reg_free_volume("free_volume", &free_volume_caps, &free_volume_create);
static MeasureRegistrar g_reg_softness("softness_proxy", &softness_caps, &softness_create);

} // namespace
} // namespace pilots
