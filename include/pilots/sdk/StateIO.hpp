#pragma once

#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "pilots/util/BinaryIO.hpp"

namespace pilots::sdk {

// Measure state IO helpers.
//
// Goals:
//   - centralize magic/version/endianness checks
//   - provide safe POD/vector/string serialization
//   - give measures a place to store "state_fingerprint" metadata (so resume errors are fail-fast)
//
// This is intentionally light-weight: algorithms live in K; this is only for robust assembly.

inline constexpr std::uint32_t kStateEndianMarker = 0x01020304u;

struct StateHeader {
  std::string magic;
  std::uint32_t version = 1;
  std::string measure_type;
  std::string instance;
  std::string fingerprint; // user-defined canonical string; should match the running config
};

class StateWriter {
public:
  explicit StateWriter(std::ostream& os)
      : os_(os), bw_(os_) {}

  void begin(const StateHeader& h) {
    if (h.magic.empty()) throw std::runtime_error("StateWriter: magic is empty");
    pilots::util::write_magic(os_, h.magic);
    bw_.write_u32(h.version);
    bw_.write_u32(kStateEndianMarker);
    bw_.write_string(h.measure_type);
    bw_.write_string(h.instance);
    bw_.write_string(h.fingerprint);
  }

  void end(std::string_view end_magic = "PILOTSEND") {
    pilots::util::write_magic(os_, end_magic);
  }

  // Forwarders
  void write_u8(std::uint8_t v) { bw_.write_u8(v); }
  void write_u32(std::uint32_t v) { bw_.write_u32(v); }
  void write_u64(std::uint64_t v) { bw_.write_u64(v); }
  void write_i32(std::int32_t v) { bw_.write_i32(v); }
  void write_i64(std::int64_t v) { bw_.write_i64(v); }
  void write_f64(double v) { bw_.write_f64(v); }
  void write_string(std::string_view s) { bw_.write_string(s); }

  template <typename T>
  void write_vec_pod(const std::vector<T>& v) { bw_.write_vec_pod(v); }

private:
  std::ostream& os_;
  pilots::util::BinaryWriter bw_;
};

class StateReader {
public:
  explicit StateReader(std::istream& is)
      : is_(is), br_(is_) {}

  StateHeader begin(std::string_view expected_magic,
                    const std::string& expected_measure_type,
                    const std::string& expected_instance,
                    const std::string& expected_fingerprint,
                    std::uint32_t expected_version = 1) {
    pilots::util::require_magic(is_, expected_magic);
    StateHeader h;
    h.magic = std::string(expected_magic);
    h.version = br_.read_u32();
    if (h.version != expected_version) {
      throw std::runtime_error("StateReader: unsupported state version (expected " + std::to_string(expected_version) + ", got " + std::to_string(h.version) + ")");
    }
    const std::uint32_t endian = br_.read_u32();
    if (endian != kStateEndianMarker) {
      throw std::runtime_error("StateReader: endianness mismatch (state file not readable on this architecture)");
    }
    h.measure_type = br_.read_string();
    h.instance = br_.read_string();
    h.fingerprint = br_.read_string();

    if (!expected_measure_type.empty() && h.measure_type != expected_measure_type) {
      throw std::runtime_error("StateReader: measure_type mismatch in checkpoint (expected '" + expected_measure_type + "', got '" + h.measure_type + "')");
    }
    if (!expected_instance.empty() && h.instance != expected_instance) {
      throw std::runtime_error("StateReader: instance mismatch in checkpoint (expected '" + expected_instance + "', got '" + h.instance + "')");
    }
    if (!expected_fingerprint.empty() && h.fingerprint != expected_fingerprint) {
      throw std::runtime_error("StateReader: fingerprint mismatch in checkpoint (configuration changed?)");
    }

    return h;
  }

  void end(std::string_view end_magic = "PILOTSEND") {
    pilots::util::require_magic(is_, end_magic);
  }

  // Forwarders
  std::uint8_t read_u8() { return br_.read_u8(); }
  std::uint32_t read_u32() { return br_.read_u32(); }
  std::uint64_t read_u64() { return br_.read_u64(); }
  std::int32_t read_i32() { return br_.read_i32(); }
  std::int64_t read_i64() { return br_.read_i64(); }
  double read_f64() { return br_.read_f64(); }
  std::string read_string() { return br_.read_string(); }

  template <typename T>
  void read_vec_pod(std::vector<T>& v) { br_.read_vec_pod(v); }

private:
  std::istream& is_;
  pilots::util::BinaryReader br_;
};

} // namespace pilots::sdk
