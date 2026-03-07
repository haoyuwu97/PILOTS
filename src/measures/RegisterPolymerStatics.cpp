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

using measure_ext::SelectedChains;
using measure_ext::OrderedChains;
using measure_ext::build_selected_chains;
using measure_ext::chain_id_per_atom_from_config;
using measure_ext::count_diag_dims;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::ordered_chains_from_config;
using measure_ext::parse_diag_mask;
using measure_ext::resolve_path;

struct RGReeOptions {
  std::int64_t frame_start = 0;
  std::int64_t frame_end = -1;
  int diag_mask = 7;
  double dt = 0.0;
  bool dry_run = false;
};

class RGReeMeasure final : public IMeasure {
public:
  struct Row {
    std::size_t frame = 0;
    std::int64_t timestep = 0;
    double time = 0.0;
    std::size_t n_chains = 0;
    double mean_rg2 = 0.0;
    double mean_rg = 0.0;
    double mean_ree2 = 0.0;
    double mean_ree = 0.0;
  };

  RGReeMeasure(std::string instance_name,
               std::string output_path,
               SelectionView sel,
               OrderedChains chains,
               RGReeOptions opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        chains_(std::move(chains)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};

    if (opt_.frame_start < 0) throw std::runtime_error("rg_ree: frame_start must be >= 0");
    if (opt_.frame_end >= 0 && opt_.frame_end < opt_.frame_start) {
      throw std::runtime_error("rg_ree: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.dt > 0.0)) throw std::runtime_error("rg_ree: dt must be > 0");
    if (opt_.diag_mask <= 0 || opt_.diag_mask > 7) throw std::runtime_error("rg_ree: invalid diag_mask");
    if (sel_.idx.empty()) throw std::runtime_error("rg_ree: selection is empty");
    if (chains_.n_chains() == 0) throw std::runtime_error("rg_ree: no ordered chains are available in the selected atoms");
  }

  std::string type() const override { return "rg_ree"; }
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
    od.x_axis = "frame";
    od.x_unit = "frames";
    od.columns = {"frame", "timestep", "time", "n_chains", "mean_rg2", "mean_rg", "mean_ree2", "mean_ree"};
    md.outputs.push_back(std::move(od));

    md.params["dt"] = dstr(opt_.dt);
    md.params["frame_start"] = std::to_string(opt_.frame_start);
    md.params["frame_end"] = std::to_string(opt_.frame_end);
    md.params["diag_mask"] = std::to_string(opt_.diag_mask);
    md.params["dimensions"] = std::to_string(count_diag_dims(opt_.diag_mask));
    md.params["n_ordered_chains"] = std::to_string(chains_.n_chains());
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.frame_start) return;
    if (opt_.frame_end >= 0 && fi > opt_.frame_end) return;

    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");

    const bool use_x = (opt_.diag_mask & 1) != 0;
    const bool use_y = (opt_.diag_mask & 2) != 0;
    const bool use_z = (opt_.diag_mask & 4) != 0;

    double sum_rg2 = 0.0;
    double sum_ree2 = 0.0;

    for (std::size_t c = 0; c < chains_.n_chains(); ++c) {
      const auto& chain = chains_.atom_by_chain[c];
      const double invn = 1.0 / static_cast<double>(chain.size());

      double cx = 0.0, cy = 0.0, cz = 0.0;
      for (const std::size_t i : chain) {
        cx += xu[i];
        cy += yu[i];
        cz += zu[i];
      }
      cx *= invn;
      cy *= invn;
      cz *= invn;

      double rg2 = 0.0;
      for (const std::size_t i : chain) {
        const double dx = xu[i] - cx;
        const double dy = yu[i] - cy;
        const double dz = zu[i] - cz;
        if (use_x) rg2 += dx * dx;
        if (use_y) rg2 += dy * dy;
        if (use_z) rg2 += dz * dz;
      }
      rg2 *= invn;
      sum_rg2 += rg2;

      const std::size_t first = chain.front();
      const std::size_t last = chain.back();
      const double dx = xu[last] - xu[first];
      const double dy = yu[last] - yu[first];
      const double dz = zu[last] - zu[first];
      double ree2 = 0.0;
      if (use_x) ree2 += dx * dx;
      if (use_y) ree2 += dy * dy;
      if (use_z) ree2 += dz * dz;
      sum_ree2 += ree2;
    }

    const double invc = 1.0 / static_cast<double>(chains_.n_chains());
    const double mean_rg2 = sum_rg2 * invc;
    const double mean_ree2 = sum_ree2 * invc;

    rows_.push_back(Row{
        frame_index,
        frame.timestep,
        static_cast<double>(frame.timestep) * opt_.dt,
        chains_.n_chains(),
        mean_rg2,
        (mean_rg2 > 0.0) ? std::sqrt(mean_rg2) : 0.0,
        mean_ree2,
        (mean_ree2 > 0.0) ? std::sqrt(mean_ree2) : 0.0,
    });
  }

  void flush_partial() override {
    if (opt_.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      for (const auto& r : rows_) {
        ofs << r.frame << ' '
            << r.timestep << ' '
            << std::setprecision(17) << r.time << ' '
            << r.n_chains << ' '
            << std::setprecision(17) << r.mean_rg2 << ' '
            << std::setprecision(17) << r.mean_rg << ' '
            << std::setprecision(17) << r.mean_ree2 << ' '
            << std::setprecision(17) << r.mean_ree << '\n';
      }
    });
  }

  void finalize() override {
    if (opt_.dry_run) return;
    flush_partial();
  }

private:
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  SelectionView sel_;
  RGReeOptions opt_;
  OrderedChains chains_;
  std::vector<Row> rows_;
  bool started_ = false;

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: polymer size observables <Rg^2>, <Rg>, <Ree^2>, <Ree>\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# n_ordered_chains: " << chains_.n_chains() << "\n";
    ofs << "# dt: " << std::setprecision(17) << opt_.dt << " (seconds per timestep)\n";
    ofs << "# frame_range: [" << opt_.frame_start << ", " << opt_.frame_end << "]\n";
    ofs << "# components_mask: " << opt_.diag_mask << " (xx=1,yy=2,zz=4)\n";
    ofs << "# dimensions: " << count_diag_dims(opt_.diag_mask) << "\n";
    ofs << "# columns: frame  timestep  time  n_chains  mean_rg2  mean_rg  mean_ree2  mean_ree\n";
  }
};

void append_common_static_caps(const IniConfig& cfg,
                               const std::string& section,
                               MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.scale = ScaleCompatibility{true, true, true};
  const std::string group_ref = cfg.get_string(section, "group", std::optional<std::string>("all"));
  caps.group_refs.push_back(group_ref);
}

void append_chain_membership_caps(const IniConfig& cfg,
                                  const std::string& section,
                                  const MeasureBuildEnv& env,
                                  MeasureCapabilities& caps) {
  if (cfg.has_key(section, "chain_id_field")) {
    caps.requires_dfields.push_back(cfg.get_string(section, "chain_id_field"));
    return;
  }
  if (env.first_frame && env.first_frame->has_mol) {
    caps.requires_i64fields.push_back("mol");
    return;
  }
  caps.requires_topology_sections.push_back("bonds");
}

void append_chain_order_caps(const IniConfig& cfg,
                             const std::string& section,
                             const MeasureBuildEnv& env,
                             MeasureCapabilities& caps) {
  append_chain_membership_caps(cfg, section, env, caps);
  if (cfg.has_key(section, "chain_pos_field")) {
    caps.requires_dfields.push_back(cfg.get_string(section, "chain_pos_field"));
  } else {
    caps.requires_topology_sections.push_back("bonds");
  }
}

MeasureCapabilities rg_ree_caps(const IniConfig& cfg,
                                const std::string& section,
                                const std::string& instance,
                                const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  append_common_static_caps(cfg, section, caps);
  append_chain_order_caps(cfg, section, env, caps);
  return caps;
}

std::unique_ptr<IMeasure> rg_ree_create(const IniConfig& cfg,
                                        const std::string& section,
                                        const std::string& instance,
                                        const MeasureBuildEnv& env,
                                        const SystemContext& sysctx) {
  if (!env.first_frame) throw std::runtime_error("rg_ree factory: first_frame is null");
  if (!env.selection_provider) throw std::runtime_error("rg_ree factory: SelectionProvider is missing");
  const Frame& frame0 = *env.first_frame;

  const std::string group_ref = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo_group_ref = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string combine_expr = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group_ref, topo_group_ref, combine_expr, "rg_ree");

  const auto chain_id = chain_id_per_atom_from_config(cfg, section, frame0, sysctx);
  const SelectedChains groups = build_selected_chains(sel, std::span<const std::int64_t>(chain_id.data(), chain_id.size()));
  const OrderedChains chains = ordered_chains_from_config(cfg, section, frame0, sysctx, sel, groups);

  const std::int64_t frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  const std::int64_t frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  if (frame_start < 0) throw std::runtime_error("frame_start must be >= 0");
  if (frame_end >= 0 && frame_end < frame_start) {
    throw std::runtime_error("frame_end must be -1 or >= frame_start");
  }

  const int diag_mask = parse_diag_mask(cfg.get_string(section, "components", std::optional<std::string>("xxyyzz")));
  const std::string out_file = cfg.get_string(section, "output", std::optional<std::string>("rg_ree.dat"));
  const std::string outdir_s = cfg.get_string(section, "output_dir", std::optional<std::string>(""));
  const fs::path output_dir = outdir_s.empty() ? env.output_dir_general : resolve_path(env.cfg_dir, outdir_s);
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out_path = (output_dir / out_file).lexically_normal();

  RGReeOptions opt;
  opt.frame_start = frame_start;
  opt.frame_end = frame_end;
  opt.diag_mask = diag_mask;
  opt.dt = env.dt;
  opt.dry_run = env.dry_run;

  return std::make_unique<RGReeMeasure>(instance, out_path.string(), sel, chains, opt);
}

static MeasureRegistrar g_register_rg_ree("rg_ree", &rg_ree_caps, &rg_ree_create);

} // namespace
} // namespace pilots
