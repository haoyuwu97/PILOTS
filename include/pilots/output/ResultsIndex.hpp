#pragma once

#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

#include "pilots/util/AtomicFile.hpp"

namespace pilots::output {
namespace fs = std::filesystem;

// Results schema version for downstream tooling (Python, dashboards, etc.).
// Bump this when changing JSON structure in non-backward-compatible ways.
inline constexpr const char* RESULTS_SCHEMA_VERSION = "1.0";

struct FileFingerprint {
  std::string path; // as recorded in config / resolved path
  std::uint64_t size_bytes = 0;
  std::int64_t mtime_epoch_s = 0; // seconds since epoch (best-effort)

  // Hashing is optional because very large, actively-growing follow inputs can be expensive
  // or misleading to hash. When hash_computed=false, hash_fnv1a64_hex is empty.
  bool hash_computed = false;
  std::string hash_fnv1a64_hex;
  std::string hash_kind = "fnv1a64"; // "fnv1a64" or "none"
};

struct MappingAudit {
  bool enabled = false;
  std::string mode = "none";                 // none|identity|by_mol|file
  std::string position = "com_geom";         // com_mass|com_geom|representative_atom
  std::string source_group;                   // optional
  std::string file;                           // optional
  std::string file_hash_fnv1a64_hex;           // optional (mode=file)
  std::string spec_hash_fnv1a64_hex;           // hash of the mapping spec (stable string canonicalization)
};

struct OutputFileDescriptor {
  std::string path;                 // relative or absolute
  std::string format = "text";      // "text" for now
  // Dataset semantics (v1): the logical x-axis for downstream tooling.
  // Example: "frame" (unit: "frames"), "timestep" (unit: "timesteps"), "timebin" (unit: "simulation_time").
  std::string x_axis;
  std::string x_unit;
  std::vector<std::string> columns; // column names in order
  // Optional per-column metadata for downstream tooling.
  // If present, these vectors must be the same length as `columns`.
  std::vector<std::string> column_units;
  std::vector<std::string> column_descriptions;
};

struct MeasureDescriptor {
  std::string instance;                 // instance name, e.g. section [measure.foo] -> "foo"
  std::string type;                     // measure type, e.g. "msd"
  std::string selection;                // selection name/expr
  std::size_t n_selected = 0;           // selection size
  std::vector<OutputFileDescriptor> outputs;
  std::map<std::string, std::string> params; // key->string value (JSON as string)
};

struct ProfilingEntry {
  double seconds = 0.0;
};

struct MeasureProfiling {
  double on_start_s = 0.0;
  double on_frame_s = 0.0;
  double finalize_s = 0.0;
  std::size_t frames = 0;
};

struct GroupAudit {
  std::string name;
  std::string expr;
  bool is_dynamic = false;
  std::size_t samples = 0; // number of times (frames) this group was evaluated (dynamic)
  std::size_t size_min = 0;
  std::size_t size_max = 0;
  double size_mean = 0.0;
  double size_var = 0.0;
};

struct CombinedSelectionAudit {
  std::string key; // A=...;T=...;combine=...
  bool is_dynamic = false;
  std::size_t samples = 0;
  std::size_t size_min = 0;
  std::size_t size_max = 0;
  double size_mean = 0.0;
  double size_var = 0.0;
};

struct BondGraphAudit {
  std::string mode = "all_bonds";
  std::vector<int> bond_types; // filter
  std::string atom_group;      // atom-group filter
  std::size_t edges = 0;       // number of undirected edges retained
};

struct PolymerClassifierComponentAudit {
  std::size_t id = 0;
  std::size_t size = 0;
  std::size_t n_ends = 0;
  std::size_t n_branch_points = 0;
  std::int64_t cycle_rank = 0;
  std::string tag; // linear|ring|branched|network|other
};

struct PolymerClassifierAudit {
  std::string graph_source; // topology_graph|bond_graph|bead_graph
  std::vector<std::string> tags;
  std::string hint;
  std::size_t n_components = 0;
  std::size_t largest_component_size = 0;
  std::vector<PolymerClassifierComponentAudit> components;
};

struct ResultsIndex {
  std::string schema_version = RESULTS_SCHEMA_VERSION;
  std::string pilots_version;

  // Model/scale audit (AA/UA/CG/auto). Mapping is scaffolded here even if not yet
  // implemented, so the results schema is stable.
  std::string model_scale = "auto";
  std::string dt_unit = "simulation_time"; // physical unit is user-defined; this denotes simulation time.

  std::string config_path;
  std::string config_hash_fnv1a64_hex;
  std::string input_path;

  FileFingerprint config_fingerprint;
  FileFingerprint input_fingerprint;
  std::optional<MappingAudit> mapping;

  std::string output_dir;

  double dt = 0.0;
  int threads = 1;
  bool follow = false;

  // Topology audit (D)
  std::string topology_file;
  std::string topology_format;
  std::string topology_hash_fnv1a64_hex;
  FileFingerprint topology_fingerprint;
  std::vector<std::string> topology_loaded_sections;
  bool topology_derive_mol_from_bonds = false;
  std::optional<BondGraphAudit> bond_graph;

  // Polymer structural classifier (K7)
  std::optional<PolymerClassifierAudit> polymer_classifier;

  // Engineering loop controls (P0)
  std::size_t flush_every_frames = 0;
  double flush_every_seconds = 0.0;
  std::string checkpoint_out;
  std::string resume_from;
  std::uint64_t reader_offset = 0;

  std::size_t frames_processed = 0;
  double wall_seconds = 0.0;

  // coarse profiling
  double reader_seconds = 0.0;
  double group_setup_seconds = 0.0;

  // Group audit (C)
  std::vector<GroupAudit> atom_groups;
  std::vector<GroupAudit> topo_groups;
  std::vector<CombinedSelectionAudit> combined_selections;

  std::vector<MeasureDescriptor> measures;
  std::vector<std::pair<std::string, MeasureProfiling>> measure_profiling; // instance -> timings
};

inline std::string json_escape(const std::string& s) {
  std::ostringstream oss;
  for (unsigned char c : s) {
    switch (c) {
      case '\\': oss << "\\\\"; break;
      case '"':  oss << "\\\""; break;
      case '\b': oss << "\\b"; break;
      case '\f': oss << "\\f"; break;
      case '\n': oss << "\\n"; break;
      case '\r': oss << "\\r"; break;
      case '\t': oss << "\\t"; break;
      default:
        if (c < 0x20) {
          oss << "\\u" << std::hex << std::setw(4) << std::setfill('0') << (int)c << std::dec;
        } else {
          oss << c;
        }
    }
  }
  return oss.str();
}

inline void write_results_json(const fs::path& out_path, const ResultsIndex& idx) {
  auto q = [&](const std::string& s) {
    return std::string("\"") + json_escape(s) + "\"";
  };

  auto write_fingerprint = [&](std::ostream& ofs, const FileFingerprint& fp, int indent) {
    const std::string pad(indent, ' ');
    ofs << pad << "{\n";
    ofs << pad << "  \"path\": " << q(fp.path) << ",\n";
    ofs << pad << "  \"size_bytes\": " << fp.size_bytes << ",\n";
    ofs << pad << "  \"mtime_epoch_s\": " << fp.mtime_epoch_s << ",\n";
    ofs << pad << "  \"hash_kind\": " << q(fp.hash_kind) << ",\n";
    ofs << pad << "  \"hash_computed\": " << (fp.hash_computed ? "true" : "false") << ",\n";
    ofs << pad << "  \"hash_fnv1a64\": " << q(fp.hash_fnv1a64_hex) << "\n";
    ofs << pad << "}";
  };

  util::atomic_write_text(out_path, [&](std::ostream& ofs) {
    ofs << "{\n";
    ofs << "  \"schema_version\": " << q(idx.schema_version) << ",\n";
    ofs << "  \"pilots_version\": " << q(idx.pilots_version) << ",\n";
    ofs << "  \"run\": {\n";
    ofs << "    \"model_scale\": " << q(idx.model_scale) << ",\n";
    ofs << "    \"dt_unit\": " << q(idx.dt_unit) << ",\n";
    ofs << "    \"config_path\": " << q(idx.config_path) << ",\n";
    ofs << "    \"config_hash_fnv1a64\": " << q(idx.config_hash_fnv1a64_hex) << ",\n";
    ofs << "    \"input_path\": " << q(idx.input_path) << ",\n";
    ofs << "    \"config_fingerprint\": ";
    write_fingerprint(ofs, idx.config_fingerprint, 4);
    ofs << ",\n";
    ofs << "    \"input_fingerprint\": ";
    write_fingerprint(ofs, idx.input_fingerprint, 4);
    ofs << ",\n";

    if (idx.mapping) {
      const auto& m = *idx.mapping;
      ofs << "    \"mapping\": {\n";
      ofs << "      \"enabled\": " << (m.enabled ? "true" : "false") << ",\n";
      ofs << "      \"mode\": " << q(m.mode) << ",\n";
      ofs << "      \"position\": " << q(m.position) << ",\n";
      ofs << "      \"source_group\": " << q(m.source_group) << ",\n";
      ofs << "      \"file\": " << q(m.file) << ",\n";
      ofs << "      \"file_hash_fnv1a64\": " << q(m.file_hash_fnv1a64_hex) << ",\n";
      ofs << "      \"spec_hash_fnv1a64\": " << q(m.spec_hash_fnv1a64_hex) << "\n";
      ofs << "    },\n";
    }
    ofs << "    \"output_dir\": " << q(idx.output_dir) << ",\n";

    ofs << "    \"dt\": " << std::setprecision(17) << idx.dt << ",\n";
    ofs << "    \"threads\": " << idx.threads << ",\n";
    ofs << "    \"follow\": " << (idx.follow ? "true" : "false") << ",\n";

    ofs << "    \"flush_every_frames\": " << idx.flush_every_frames << ",\n";
    ofs << "    \"flush_every_seconds\": " << std::setprecision(17) << idx.flush_every_seconds << ",\n";
    ofs << "    \"checkpoint_out\": " << q(idx.checkpoint_out) << ",\n";
    ofs << "    \"resume_from\": " << q(idx.resume_from) << ",\n";
    ofs << "    \"reader_offset\": " << idx.reader_offset << ",\n";

    ofs << "    \"frames_processed\": " << idx.frames_processed << ",\n";
    ofs << "    \"wall_seconds\": " << std::setprecision(17) << idx.wall_seconds << ",\n";

    // Topology (D)
    ofs << "    \"topology\": {\n";
    ofs << "      \"file\": " << q(idx.topology_file) << ",\n";
    ofs << "      \"format\": " << q(idx.topology_format) << ",\n";
    ofs << "      \"hash_fnv1a64\": " << q(idx.topology_hash_fnv1a64_hex) << ",\n";
    ofs << "      \"fingerprint\": ";
    write_fingerprint(ofs, idx.topology_fingerprint, 6);
    ofs << ",\n";
    ofs << "      \"derive_mol_from_bonds\": " << (idx.topology_derive_mol_from_bonds ? "true" : "false") << ",\n";
    ofs << "      \"loaded_sections\": [";
    for (std::size_t i = 0; i < idx.topology_loaded_sections.size(); ++i) {
      ofs << q(idx.topology_loaded_sections[i]);
      if (i + 1 < idx.topology_loaded_sections.size()) ofs << ", ";
    }
    ofs << "]";
    if (idx.bond_graph) {
      const auto& bg = *idx.bond_graph;
      ofs << ",\n";
      ofs << "      \"bond_graph\": {\n";
      ofs << "        \"mode\": " << q(bg.mode) << ",\n";
      ofs << "        \"bond_types\": [";
      for (std::size_t i = 0; i < bg.bond_types.size(); ++i) {
        ofs << bg.bond_types[i];
        if (i + 1 < bg.bond_types.size()) ofs << ", ";
      }
      ofs << "],\n";
      ofs << "        \"atom_group\": " << q(bg.atom_group) << ",\n";
      ofs << "        \"edges\": " << bg.edges << "\n";
      ofs << "      }\n";
    } else {
      ofs << "\n";
    }
    ofs << "    }\n";
    ofs << "  },\n";

    ofs << "  \"profiling\": {\n";
    ofs << "    \"reader_seconds\": " << std::setprecision(17) << idx.reader_seconds << ",\n";
    ofs << "    \"group_setup_seconds\": " << std::setprecision(17) << idx.group_setup_seconds << ",\n";
    ofs << "    \"measures\": {\n";
    for (std::size_t i = 0; i < idx.measure_profiling.size(); ++i) {
      const auto& [name, mp] = idx.measure_profiling[i];
      ofs << "      " << q(name) << ": {"
          << "\"on_start_s\": " << std::setprecision(17) << mp.on_start_s << ", "
          << "\"on_frame_s\": " << std::setprecision(17) << mp.on_frame_s << ", "
          << "\"finalize_s\": " << std::setprecision(17) << mp.finalize_s << ", "
          << "\"frames\": " << mp.frames
          << "}";
      if (i + 1 < idx.measure_profiling.size()) ofs << ",";
      ofs << "\n";
    }
    ofs << "    }\n";
    ofs << "  },\n";

    // Polymer classifier (K7)
    if (idx.polymer_classifier) {
      const auto& pc = *idx.polymer_classifier;
      ofs << "  \"polymer_classifier\": {\n";
      ofs << "    \"graph_source\": " << q(pc.graph_source) << ",\n";
      ofs << "    \"tags\": [";
      for (std::size_t i = 0; i < pc.tags.size(); ++i) {
        ofs << q(pc.tags[i]);
        if (i + 1 < pc.tags.size()) ofs << ", ";
      }
      ofs << "],\n";
      ofs << "    \"hint\": " << q(pc.hint) << ",\n";
      ofs << "    \"summary\": {\n";
      ofs << "      \"n_components\": " << pc.n_components << ",\n";
      ofs << "      \"largest_component_size\": " << pc.largest_component_size << "\n";
      ofs << "    },\n";
      ofs << "    \"components\": [\n";
      for (std::size_t ci = 0; ci < pc.components.size(); ++ci) {
        const auto& c = pc.components[ci];
        ofs << "      {\n";
        ofs << "        \"id\": " << c.id << ",\n";
        ofs << "        \"size\": " << c.size << ",\n";
        ofs << "        \"n_ends\": " << c.n_ends << ",\n";
        ofs << "        \"n_branch_points\": " << c.n_branch_points << ",\n";
        ofs << "        \"cycle_rank\": " << c.cycle_rank << ",\n";
        ofs << "        \"tag\": " << q(c.tag) << "\n";
        ofs << "      }";
        if (ci + 1 < pc.components.size()) ofs << ",";
        ofs << "\n";
      }
      ofs << "    ]\n";
      ofs << "  },\n";
    }

    // Group audits (C)
    ofs << "  \"groups\": {\n";
    auto write_group_list = [&](const char* key, const std::vector<GroupAudit>& gs) {
      ofs << "    \"" << key << "\": [\n";
      for (std::size_t gi = 0; gi < gs.size(); ++gi) {
        const auto& g = gs[gi];
        ofs << "      {\n";
        ofs << "        \"name\": " << q(g.name) << ",\n";
        ofs << "        \"expr\": " << q(g.expr) << ",\n";
        ofs << "        \"is_dynamic\": " << (g.is_dynamic ? "true" : "false") << ",\n";
        ofs << "        \"samples\": " << g.samples << ",\n";
        ofs << "        \"size_min\": " << g.size_min << ",\n";
        ofs << "        \"size_max\": " << g.size_max << ",\n";
        ofs << "        \"size_mean\": " << std::setprecision(17) << g.size_mean << ",\n";
        ofs << "        \"size_var\": " << std::setprecision(17) << g.size_var << "\n";
        ofs << "      }";
        if (gi + 1 < gs.size()) ofs << ",";
        ofs << "\n";
      }
      ofs << "    ]";
    };
    write_group_list("atom_groups", idx.atom_groups);
    ofs << ",\n";
    write_group_list("topo_groups", idx.topo_groups);
    ofs << ",\n";
    ofs << "    \"combined_selections\": [\n";
    for (std::size_t si = 0; si < idx.combined_selections.size(); ++si) {
      const auto& s = idx.combined_selections[si];
      ofs << "      {\n";
      ofs << "        \"key\": " << q(s.key) << ",\n";
      ofs << "        \"is_dynamic\": " << (s.is_dynamic ? "true" : "false") << ",\n";
      ofs << "        \"samples\": " << s.samples << ",\n";
      ofs << "        \"size_min\": " << s.size_min << ",\n";
      ofs << "        \"size_max\": " << s.size_max << ",\n";
      ofs << "        \"size_mean\": " << std::setprecision(17) << s.size_mean << ",\n";
      ofs << "        \"size_var\": " << std::setprecision(17) << s.size_var << "\n";
      ofs << "      }";
      if (si + 1 < idx.combined_selections.size()) ofs << ",";
      ofs << "\n";
    }
    ofs << "    ]\n";
    ofs << "\n";
    ofs << "  },\n";

    ofs << "  \"measures\": [\n";
    for (std::size_t mi = 0; mi < idx.measures.size(); ++mi) {
      const auto& md = idx.measures[mi];
      ofs << "    {\n";
      ofs << "      \"instance\": " << q(md.instance) << ",\n";
      ofs << "      \"type\": " << q(md.type) << ",\n";
      ofs << "      \"selection\": " << q(md.selection) << ",\n";
      ofs << "      \"n_selected\": " << md.n_selected << ",\n";
      ofs << "      \"params\": {\n";
      std::size_t pk = 0;
      for (const auto& kv : md.params) {
        ofs << "        " << q(kv.first) << ": " << q(kv.second);
        if (++pk < md.params.size()) ofs << ",";
        ofs << "\n";
      }
      ofs << "      },\n";
      ofs << "      \"outputs\": [\n";
      for (std::size_t oi = 0; oi < md.outputs.size(); ++oi) {
        const auto& od = md.outputs[oi];
        ofs << "        {\n";
        ofs << "          \"path\": " << q(od.path) << ",\n";
        ofs << "          \"format\": " << q(od.format) << ",\n";
        ofs << "          \"x_axis\": " << q(od.x_axis) << ",\n";
        ofs << "          \"x_unit\": " << q(od.x_unit) << ",\n";
        ofs << "          \"columns\": [";
        for (std::size_t ci = 0; ci < od.columns.size(); ++ci) {
          ofs << q(od.columns[ci]);
          if (ci + 1 < od.columns.size()) ofs << ", ";
        }
        ofs << "],\n";

        const std::size_t ncols = od.columns.size();
        ofs << "          \"column_units\": [";
        for (std::size_t ci = 0; ci < ncols; ++ci) {
          if (ci < od.column_units.size()) {
            ofs << q(od.column_units[ci]);
          } else {
            ofs << q("");
          }
          if (ci + 1 < ncols) ofs << ", ";
        }
        ofs << "],\n";

        ofs << "          \"column_descriptions\": [";
        for (std::size_t ci = 0; ci < ncols; ++ci) {
          if (ci < od.column_descriptions.size()) {
            ofs << q(od.column_descriptions[ci]);
          } else {
            ofs << q("");
          }
          if (ci + 1 < ncols) ofs << ", ";
        }
        ofs << "]\n";
        ofs << "        }";
        if (oi + 1 < md.outputs.size()) ofs << ",";
        ofs << "\n";
      }
      ofs << "      ]\n";
      ofs << "    }";
      if (mi + 1 < idx.measures.size()) ofs << ",";
      ofs << "\n";
    }
    ofs << "  ]\n";
    ofs << "}\n";
  });
}

} // namespace pilots::output
