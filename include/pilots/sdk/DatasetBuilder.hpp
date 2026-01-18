#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "pilots/output/ResultsIndex.hpp"

namespace pilots::sdk {

struct ColumnSpec {
  std::string name;
  std::string unit;        // optional
  std::string description; // optional
};

// DatasetBuilder: declares a dataset's output contract.
//
// This object does NOT write data. It only standardizes:
//   - columns + (optional) units/descriptions
//   - x_axis semantics (frame/timestep/timebin/seconds...)
//   - results.json descriptor generation
//
// Measures can use this to keep output + audit consistent.
class DatasetBuilder {
public:
  DatasetBuilder() = default;

  explicit DatasetBuilder(std::string path)
      : path_(std::move(path)) {}

  DatasetBuilder& path(std::string p) {
    path_ = std::move(p);
    return *this;
  }

  const std::string& path() const { return path_; }

  DatasetBuilder& format(std::string fmt) {
    format_ = std::move(fmt);
    return *this;
  }

  DatasetBuilder& x_axis(std::string axis, std::string unit) {
    x_axis_ = std::move(axis);
    x_unit_ = std::move(unit);
    return *this;
  }

  DatasetBuilder& column(std::string name, std::string unit = {}, std::string description = {}) {
    if (name.empty()) throw std::runtime_error("DatasetBuilder: column name is empty");
    cols_.push_back(ColumnSpec{std::move(name), std::move(unit), std::move(description)});
    return *this;
  }

  const std::vector<ColumnSpec>& columns() const { return cols_; }

  // Generate a one-line header suitable for text outputs.
  // Example: "# lag\ttime(s)\tF...".
  std::string header_line(std::string_view prefix = "# ") const {
    std::string out(prefix);
    for (std::size_t i = 0; i < cols_.size(); ++i) {
      const auto& c = cols_[i];
      out += c.name;
      if (!c.unit.empty()) {
        out += "(";
        out += c.unit;
        out += ")";
      }
      if (i + 1 < cols_.size()) out += "\t";
    }
    out += "\n";
    return out;
  }

  pilots::output::OutputFileDescriptor build_descriptor() const {
    if (path_.empty()) throw std::runtime_error("DatasetBuilder: path is empty");
    if (cols_.empty()) throw std::runtime_error("DatasetBuilder: no columns declared");

    pilots::output::OutputFileDescriptor od;
    od.path = path_;
    od.format = format_;
    od.x_axis = x_axis_;
    od.x_unit = x_unit_;

    od.columns.reserve(cols_.size());
    od.column_units.reserve(cols_.size());
    od.column_descriptions.reserve(cols_.size());

    for (const auto& c : cols_) {
      od.columns.push_back(c.name);
      od.column_units.push_back(c.unit);
      od.column_descriptions.push_back(c.description);
    }

    return od;
  }

private:
  std::string path_;
  std::string format_ = "text";
  std::string x_axis_;
  std::string x_unit_;
  std::vector<ColumnSpec> cols_;
};

// Helper builder for output::MeasureDescriptor.
class MeasureDescriptorBuilder {
public:
  MeasureDescriptorBuilder(std::string instance, std::string type)
      : d_{} {
    d_.instance = std::move(instance);
    d_.type = std::move(type);
  }

  MeasureDescriptorBuilder& selection(std::string label, std::size_t n_selected) {
    d_.selection = std::move(label);
    d_.n_selected = n_selected;
    return *this;
  }

  MeasureDescriptorBuilder& output(const DatasetBuilder& ds) {
    d_.outputs.push_back(ds.build_descriptor());
    return *this;
  }

  MeasureDescriptorBuilder& param(std::string k, std::string v) {
    d_.params[std::move(k)] = std::move(v);
    return *this;
  }

  pilots::output::MeasureDescriptor build() const { return d_; }

private:
  pilots::output::MeasureDescriptor d_;
};

} // namespace pilots::sdk
