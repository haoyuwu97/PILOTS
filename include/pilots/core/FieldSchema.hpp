#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pilots {

// Minimal schema for per-atom numeric fields.
//
// Current design (baseline):
// - Only scalar double fields are supported for extra per-atom columns.
// - dtype/arity are kept for forward compatibility.
//
// The main goal is to resolve field names to integer ids once, so hot loops avoid
// string hash lookups.
class FieldSchema {
public:
  enum class DType : int { F64 = 1 };

  struct FieldInfo {
    std::string name;
    DType dtype = DType::F64;
    int arity = 1; // 1=scalar, 3=vector, 6=tensor... (future)
  };

  // Return existing field id, or create a new field entry.
  std::size_t ensure(const std::string& name, DType dtype = DType::F64, int arity = 1) {
    auto it = name2id_.find(name);
    if (it != name2id_.end()) {
      // Basic consistency checks for forward compatibility.
      const auto& fi = fields_[it->second];
      if (fi.dtype != dtype || fi.arity != arity) {
        throw std::runtime_error("FieldSchema: field '" + name + "' registered with conflicting dtype/arity");
      }
      return it->second;
    }
    const std::size_t id = fields_.size();
    fields_.push_back(FieldInfo{name, dtype, arity});
    name2id_.emplace(fields_.back().name, id);
    return id;
  }

  // Return field id if present, else throw.
  std::size_t require(const std::string& name) const {
    auto it = name2id_.find(name);
    if (it == name2id_.end()) {
      throw std::runtime_error("FieldSchema: required field '" + name + "' not found in schema");
    }
    return it->second;
  }

  bool has(const std::string& name) const {
    return name2id_.find(name) != name2id_.end();
  }

  const FieldInfo& info(std::size_t id) const {
    if (id >= fields_.size()) throw std::runtime_error("FieldSchema: invalid field id");
    return fields_[id];
  }

  std::size_t size() const { return fields_.size(); }

private:
  std::vector<FieldInfo> fields_;
  std::unordered_map<std::string, std::size_t> name2id_;
};

} // namespace pilots
