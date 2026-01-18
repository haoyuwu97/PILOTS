#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "pilots/core/Box.hpp"
#include "pilots/core/FieldSchema.hpp"

namespace pilots {

struct Frame {
  std::int64_t timestep = 0;
  std::size_t natoms = 0;
  Box box;

  // Optional core metadata presence (based on per-frame ATOMS header)
  bool has_type = false;
  bool has_mol = false;

  // Canonical per-atom metadata and coordinates (canonical order = first frame file order)
  std::vector<std::int64_t> id;
  // Optional arrays: allocated only when the corresponding field is present.
  std::vector<int> type;
  std::vector<std::int64_t> mol;

  // Unwrapped coordinates (SoA)
  std::vector<double> xu;
  std::vector<double> yu;
  std::vector<double> zu;

  // Extra per-atom scalar fields stored by field_id (schema owned by reader).
  const FieldSchema* extra_schema = nullptr;
  std::vector<std::vector<double>> extra_values;   // [field_id][i]
  std::vector<unsigned char> extra_present;        // [field_id] 0/1 for current frame

  void resize(std::size_t n) {
    natoms = n;
    id.resize(n);
    xu.resize(n);
    yu.resize(n);
    zu.resize(n);
    // Optional core arrays are allocated on-demand (after ATOMS header parsing).
    if (has_type) {
      type.resize(n);
    } else {
      type.clear();
    }
    if (has_mol) {
      mol.resize(n);
    } else {
      mol.clear();
    }
    // extra_values handled by reader via ensure_extra_storage
  }

  void ensure_type_storage() {
    if (!has_type) return;
    if (type.size() != natoms) type.resize(natoms);
  }

  void ensure_mol_storage() {
    if (!has_mol) return;
    if (mol.size() != natoms) mol.resize(natoms);
  }

  void bind_extra_schema(const FieldSchema* schema) {
    extra_schema = schema;
    if (!extra_schema) return;
    if (extra_values.size() < extra_schema->size()) {
      extra_values.resize(extra_schema->size());
    }
    if (extra_present.size() < extra_schema->size()) {
      extra_present.resize(extra_schema->size(), 0);
    }
  }

  void clear_extra_presence() {
    if (extra_present.empty()) return;
    std::fill(extra_present.begin(), extra_present.end(), static_cast<unsigned char>(0));
  }

  void mark_extra_present(std::size_t fid) {
    if (!extra_schema) throw std::runtime_error("Frame: extra_schema is not bound");
    bind_extra_schema(extra_schema);
    if (fid >= extra_present.size()) {
      throw std::runtime_error("Frame: extra field id out of range");
    }
    extra_present[fid] = 1;
  }

  std::vector<double>& ensure_extra_storage(std::size_t fid) {
    if (!extra_schema) throw std::runtime_error("Frame: extra_schema is not bound");
    bind_extra_schema(extra_schema);
    if (fid >= extra_values.size()) {
      throw std::runtime_error("Frame: extra field id out of range");
    }
    auto& v = extra_values[fid];
    if (v.size() != natoms) v.resize(natoms);
    return v;
  }

  std::size_t require_extra_id(const std::string& name) const {
    if (!extra_schema) {
      throw std::runtime_error("Frame: no extra schema bound; cannot resolve field '" + name + "'");
    }
    return extra_schema->require(name);
  }

  std::span<const double> require_extra_by_id(std::size_t fid) const {
    if (!extra_schema) throw std::runtime_error("Frame: extra_schema is not bound");
    if (fid >= extra_values.size()) throw std::runtime_error("Frame: extra field id out of range");
    if (fid >= extra_present.size() || extra_present[fid] == 0) {
      const std::string fname = extra_schema->info(fid).name;
      throw std::runtime_error("Frame: required field '" + fname + "' is not present in this frame's ATOMS header");
    }
    const auto& v = extra_values[fid];
    if (v.size() != natoms) {
      const std::string fname = extra_schema->info(fid).name;
      throw std::runtime_error("Frame: field '" + fname + "' has wrong length");
    }
    return {v.data(), v.size()};
  }

  std::span<const double> require_dfield(const std::string& name) const {
    if (name == "xu") {
      if (xu.size() != natoms) throw std::runtime_error("Frame: missing or invalid field 'xu'");
      return {xu.data(), xu.size()};
    }
    if (name == "yu") {
      if (yu.size() != natoms) throw std::runtime_error("Frame: missing or invalid field 'yu'");
      return {yu.data(), yu.size()};
    }
    if (name == "zu") {
      if (zu.size() != natoms) throw std::runtime_error("Frame: missing or invalid field 'zu'");
      return {zu.data(), zu.size()};
    }
    const std::size_t fid = require_extra_id(name);
    return require_extra_by_id(fid);
  }

  std::span<const double> require_dfield_by_id(std::size_t fid) const {
    return require_extra_by_id(fid);
  }

  std::span<const std::int64_t> require_i64field(const std::string& name) const {
    if (name == "id") {
      if (id.size() != natoms) throw std::runtime_error("Frame: missing or invalid field 'id'");
      return {id.data(), id.size()};
    }
    if (name == "mol") {
      if (!has_mol) throw std::runtime_error("Frame: required field 'mol' is not present in dump");
      if (mol.size() != natoms) throw std::runtime_error("Frame: missing or invalid field 'mol'");
      return {mol.data(), mol.size()};
    }
    throw std::runtime_error("Frame: unknown i64 field '" + name + "'");
  }

  std::span<const int> require_intfield(const std::string& name) const {
    if (name == "type") {
      if (!has_type) throw std::runtime_error("Frame: required field 'type' is not present in dump");
      if (type.size() != natoms) throw std::runtime_error("Frame: missing or invalid field 'type'");
      return {type.data(), type.size()};
    }
    throw std::runtime_error("Frame: unknown int field '" + name + "'");
  }
};

} // namespace pilots
