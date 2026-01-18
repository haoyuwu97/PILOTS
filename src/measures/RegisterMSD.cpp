#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

#include "pilots/correlate/CorrelatorFactory.hpp"
#include "pilots/measures/MSDMeasure.hpp"
#include "pilots/measures/MeasureRegistry.hpp"

namespace fs = std::filesystem;

namespace pilots {
namespace {

fs::path resolve_path(const fs::path& base_dir, const std::string& p) {
  fs::path path(p);
  if (path.is_absolute()) return path;
  return (base_dir / path).lexically_normal();
}

int parse_diag_mask(const std::string& s_raw) {
  std::string s = s_raw;
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  // Tensor-style names
  if (s == "xxyyzz" || s == "all") return 7;
  if (s == "xx") return 1;
  if (s == "yy") return 2;
  if (s == "zz") return 4;
  if (s == "xxyy") return 3;
  if (s == "xxzz") return 5;
  if (s == "yyzz") return 6;
  // Backward-compatible aliases
  if (s == "xyz") return 7;
  if (s == "x") return 1;
  if (s == "y") return 2;
  if (s == "z") return 4;
  if (s == "xy") return 3;
  if (s == "xz") return 5;
  if (s == "yz") return 6;
  throw std::runtime_error("invalid components string: '" + s_raw + "' (use xxyyzz, xx, yy, zz, xxyy, xxzz, yyzz)");
}

std::size_t count_lammps_frames(const std::string& path) {
  std::ifstream ifs(path);
  if (!ifs) {
    throw std::runtime_error("failed to open input for frame count: " + path);
  }
  std::string line;
  std::size_t n = 0;
  while (std::getline(ifs, line)) {
    if (line.rfind("ITEM: TIMESTEP", 0) == 0) ++n;
  }
  return n;
}

MeasureCapabilities msd_caps(const IniConfig& cfg,
                            const std::string& section,
                            const std::string& instance,
                            const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;

  MeasureCapabilities caps;
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.scale = ScaleCompatibility{true, true, true};

  const std::string group_ref = cfg.get_string(section, "group", std::optional<std::string>("all"));
  caps.group_refs.push_back(group_ref);

  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  if (remove_drift) {
    const std::string drift_group_ref = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
    caps.group_refs.push_back(drift_group_ref);
  }
  return caps;
}

std::unique_ptr<IMeasure> msd_create(const IniConfig& cfg,
                                    const std::string& section,
                                    const std::string& instance,
                                    const MeasureBuildEnv& env,
                                    const SystemContext& sysctx) {
  (void)sysctx;

  if (!env.first_frame) throw std::runtime_error("MSD factory: first_frame is null");
  if (!env.selection_provider) {
    throw std::runtime_error("MSD factory: SelectionProvider is missing (Runner must build it)");
  }
  const Frame& frame0 = *env.first_frame;

  const std::string group_ref = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo_group_ref = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string combine_expr = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  const std::string drift_group_ref = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));

  const CorrelatorSpec corr_spec = parse_correlator_spec(cfg, section, env.dt);

  const std::int64_t frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  const std::int64_t frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  const std::string out_file = cfg.get_string(section, "output", std::optional<std::string>("msd.dat"));
  const std::string outdir_s = cfg.get_string(section, "output_dir", std::optional<std::string>(""));
  const std::string comp_s = cfg.get_string(section, "components", std::optional<std::string>("xxyyzz"));

  if (frame_start < 0) throw std::runtime_error("frame_start must be >= 0");
  if (frame_end >= 0 && frame_end < frame_start) throw std::runtime_error("frame_end must be -1 or >= frame_start");

  const fs::path output_dir = outdir_s.empty() ? env.output_dir_general : resolve_path(env.cfg_dir, outdir_s);
  if (!env.dry_run) {
    fs::create_directories(output_dir);
  }
  const fs::path out_path = (output_dir / out_file).lexically_normal();

  // Resolve effective frame_end for correlator=exact (finite window required).
  std::int64_t frame_end_eff = frame_end;
  if (corr_spec.type == "exact") {
    if (env.follow) {
      throw std::runtime_error("correlator=exact is not supported in follow mode (frame_end is unbounded). Use correlator=multitau.");
    }
    if (frame_end_eff < 0) {
      const std::size_t total_frames = count_lammps_frames(env.input_path.string());
      if (total_frames == 0) throw std::runtime_error("input contains 0 frames");
      frame_end_eff = static_cast<std::int64_t>(total_frames - 1);
    }
  }

  const int diag_mask = parse_diag_mask(comp_s);

  // Selections (AtomGroup + TopoGroup + combine).
  // MSD is scientifically meaningful only for identity-consistent, static selections.
  if (env.selection_provider->is_dynamic_spec(group_ref, topo_group_ref)) {
    throw std::runtime_error("MSD requires a static selection; the selection spec depends on a dynamic group/topo_group");
  }
  SelectionView sel = env.selection_provider->get_combined_view(frame0, 0, group_ref, topo_group_ref, combine_expr);

  SelectionView drift_sel;
  if (remove_drift) {
    if (env.selection_provider->is_dynamic_spec(drift_group_ref, "all")) {
      throw std::runtime_error("MSD drift_group requires static selection; drift_group depends on a dynamic group");
    }
    drift_sel = env.selection_provider->get_combined_view(frame0, 0, drift_group_ref, "all", "A");
  } else {
    drift_sel = env.selection_provider->get_combined_view(frame0, 0, "all", "all", "A");
  }

  MSDMeasure::Options opt;
  opt.frame_start = frame_start;
  opt.frame_end = (corr_spec.type == "exact") ? frame_end_eff : frame_end;
  opt.diag_mask = diag_mask;
  opt.remove_drift = remove_drift;
  opt.corr = corr_spec;
  opt.dry_run = env.dry_run;

  return std::make_unique<MSDMeasure>(instance, out_path.string(), sel, drift_sel, opt);
}

static MeasureRegistrar g_register_msd("msd", &msd_caps, &msd_create);

} // namespace
} // namespace pilots
