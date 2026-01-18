// TEMPLATE MEASURE (copy-paste starter)
//
// This file is intentionally NOT registered by default.
//
// Goal:
//   - Make it easy to add a new measure by copying this file to
//     src/measures/<MyMeasure>.cpp and editing in ONE place.
//   - With the current CMake setup (src/measures/*.cpp glob), dropping the new
//     file into src/measures/ is sufficient to compile it into the binary.
//
// How to create a new measure from this template:
//   1) Copy this file, e.g.
//        cp src/measures/template_measure.cpp src/measures/isf.cpp
//   2) Replace:
//        - kType string
//        - class name
//        - capabilities (requires_fields/topology/selection_policy/etc.)
//        - create() parsing (read your config keys)
//        - on_frame() physics
//   3) Enable registration by uncommenting the MeasureRegistrar line at the
//      bottom of the file.
//   4) Add a config section:
//        [measure.my_isf]
//        type = isf
//        enabled = true
//        ...
//
// Notes:
//   - Keep measure logic scientifically strict: if a required field/topology
//     section is missing, fail-fast.
//   - Prefer env.selection_provider to unify AtomGroup/TopoGroup selection.

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "pilots/measures/IMeasure.hpp"
#include "pilots/measures/MeasureRegistry.hpp"
#include "pilots/sdk/SDK.hpp"

namespace pilots {

namespace {

constexpr const char* kType = "template"; // CHANGE ME

class TemplateMeasure final : public IMeasure {
public:
  TemplateMeasure(std::string instance,
                  sdk::SelectionHandle selection,
                  sdk::DatasetBuilder dataset,
                  sdk::TextWriter writer)
      : instance_(std::move(instance)),
        sel_(std::move(selection)),
        dataset_(std::move(dataset)),
        writer_(std::move(writer)) {}

  std::string type() const override { return kType; }
  std::string instance_name() const override { return instance_; }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    // Example: write an initial header. For append mode, this only writes if the
    // file is missing/empty (avoids repeated headers during online flushes).
    writer_.write_with_header_if_empty(
        [&](std::ostream& os) { os << dataset_.header_line(); },
        [&](std::ostream& os) { (void)os; });
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    last_frame_ = frame_index;
    const SelectionView sel = sel_.view(frame, frame_index);
    // Example payload: just record selection size over time.
    // Replace with your physics.
    buffer_.push_back({frame_index, sel.size()});
  }

  void flush_partial() override {
    if (buffer_.empty()) return;
    // Default is safe/atomic rewrite. For high-frequency output you may prefer
    // append + checkpointed last_written_index, but that is a separate policy.
    writer_.write_with_header_if_empty(
        [&](std::ostream& os) { os << dataset_.header_line(); },
        [&](std::ostream& os) {
          for (const auto& row : buffer_) {
            os << row.frame_index << "\t" << row.n_selected << "\n";
          }
        });
  }

  void finalize() override {
    flush_partial();
  }

  output::MeasureDescriptor describe() const override {
    // If the selection is dynamic, n_selected is ambiguous; we keep the latest
    // observed size as a lightweight indicator.
    sdk::MeasureDescriptorBuilder mb(instance_, kType);
    mb.selection(sel_.key(), buffer_.empty() ? 0 : buffer_.back().n_selected);
    mb.output(dataset_);
    mb.param("note", "template_measure (copy and modify)");
    mb.param("writer_mode", sdk::text_write_mode_name(writer_.mode()));
    return mb.build();
  }

  void save_state(std::ostream& os) const override {
    sdk::StateWriter w(os);
    sdk::StateHeader h;
    h.magic = "PILOTSTEMPLATE";
    h.version = 1;
    h.measure_type = kType;
    h.instance = instance_;
    h.fingerprint = sel_.key();
    w.begin(h);
    w.write_u64(static_cast<std::uint64_t>(last_frame_));
    w.end();
  }

  void load_state(std::istream& is) override {
    sdk::StateReader r(is);
    const auto h = r.begin("PILOTSTEMPLATE", kType, instance_, sel_.key(), 1);
    (void)h;
    last_frame_ = static_cast<std::size_t>(r.read_u64());
    r.end();
  }

private:
  struct Row {
    std::size_t frame_index = 0;
    std::size_t n_selected = 0;
  };

  std::string instance_;
  sdk::SelectionHandle sel_;
  sdk::DatasetBuilder dataset_;
  sdk::TextWriter writer_;

  std::size_t last_frame_ = 0;
  std::vector<Row> buffer_;
};

static MeasureCapabilities caps_fn(const IniConfig& cfg,
                                   const std::string& section,
                                   const std::string& instance,
                                   const MeasureBuildEnv& env) {
  (void)cfg;
  (void)section;
  (void)instance;
  (void)env;
  MeasureCapabilities c;
  c.selection_policy = SelectionPolicy::AllowDynamic;
  c.requires_identity_consistent = false;
  // Example: declare required fields/topology here.
  // c.requires_dfields = {"xu", "yu", "zu"};
  // c.requires_topology_sections = {"bonds"};
  c.scale.aa = true;
  c.scale.ua = true;
  c.scale.cg = true;
  // Example: mapping placeholder (future).
  c.mapping.requires_mapping = false;
  // Provide group refs if your config references named groups beyond 'group'.
  // c.group_refs.push_back("some_group_key");
  return c;
}

static std::unique_ptr<IMeasure> create_fn(const IniConfig& cfg,
                                           const std::string& section,
                                           const std::string& instance,
                                           const MeasureBuildEnv& env,
                                           const SystemContext& sysctx) {
  (void)sysctx;
  // Unified selection: group/topo_group/combine.
  std::string group_ref = "";
  std::string topo_ref = "";
  std::string combine = "A & T";
  if (cfg.has_key(section, "group")) group_ref = cfg.get_string(section, "group");
  if (cfg.has_key(section, "topo_group")) topo_ref = cfg.get_string(section, "topo_group");
  if (cfg.has_key(section, "combine")) combine = cfg.get_string(section, "combine");

  if (!env.selection_provider) {
    throw std::runtime_error(std::string("TemplateMeasure: selection_provider is not available; set up groups/topology first (instance '") + instance + "')");
  }

  // Output path: measure.<instance>.output_prefix if present, otherwise <output_dir>/<instance>.dat
  std::string out_name = instance + ".dat";
  if (cfg.has_key(section, "output_prefix")) out_name = cfg.get_string(section, "output_prefix");
  const std::string out_path = (env.output_dir_general / out_name).string();

  // Example dataset contract.
  sdk::DatasetBuilder ds(out_path);
  ds.format("text")
    .x_axis("frame", "frames")
    .column("frame_index", "frames", "frame index (0-based within this run)")
    .column("n_selected", "count", "combined selection size");

  // Writer policy: default atomic rewrite.
  sdk::TextWriteMode wmode = sdk::TextWriteMode::AtomicRewrite;
  if (cfg.has_key(section, "writer_mode")) {
    const std::string wm = cfg.get_string(section, "writer_mode");
    if (wm == "append") wmode = sdk::TextWriteMode::Append;
  }
  sdk::TextWriter writer(out_path, wmode);

  sdk::SelectionHandle sel(instance, env.selection_provider, group_ref, topo_ref, combine);

  return std::make_unique<TemplateMeasure>(instance, std::move(sel), std::move(ds), std::move(writer));
}

} // namespace

// IMPORTANT:
// This template is NOT registered by default to avoid shipping a placeholder
// measure type in the production registry.
//
// To enable, uncomment the line below and change kType above.
//
// static MeasureRegistrar reg(kType, &caps_fn, &create_fn);

} // namespace pilots
