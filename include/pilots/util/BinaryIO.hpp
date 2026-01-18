#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace pilots::util {

class BinaryWriter {
public:
  explicit BinaryWriter(std::ostream& os) : os_(os) {
    if (!os_) throw std::runtime_error("BinaryWriter: stream is not writable");
  }

  void write_bytes(const void* data, std::size_t n) {
    os_.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(n));
    if (!os_) throw std::runtime_error("BinaryWriter: write failed");
  }

  template <typename T>
  void write_pod(const T& v) {
    static_assert(std::is_trivially_copyable_v<T>, "write_pod requires trivially copyable type");
    write_bytes(&v, sizeof(T));
  }

  void write_u8(std::uint8_t v) { write_pod(v); }
  void write_u32(std::uint32_t v) { write_pod(v); }
  void write_u64(std::uint64_t v) { write_pod(v); }
  void write_i32(std::int32_t v) { write_pod(v); }
  void write_i64(std::int64_t v) { write_pod(v); }
  void write_f64(double v) { write_pod(v); }

  void write_string(std::string_view s) {
    if (s.size() > std::numeric_limits<std::uint32_t>::max()) {
      throw std::runtime_error("BinaryWriter: string too large");
    }
    write_u32(static_cast<std::uint32_t>(s.size()));
    if (!s.empty()) write_bytes(s.data(), s.size());
  }

  template <typename T>
  void write_vec_pod(const std::vector<T>& v) {
    static_assert(std::is_trivially_copyable_v<T>, "write_vec_pod requires trivially copyable type");
    if (v.size() > std::numeric_limits<std::uint64_t>::max()) {
      throw std::runtime_error("BinaryWriter: vector too large");
    }
    write_u64(static_cast<std::uint64_t>(v.size()));
    if (!v.empty()) write_bytes(v.data(), sizeof(T) * v.size());
  }

private:
  std::ostream& os_;
};

class BinaryReader {
public:
  explicit BinaryReader(std::istream& is) : is_(is) {
    if (!is_) throw std::runtime_error("BinaryReader: stream is not readable");
  }

  void read_bytes(void* data, std::size_t n) {
    is_.read(reinterpret_cast<char*>(data), static_cast<std::streamsize>(n));
    if (!is_) throw std::runtime_error("BinaryReader: read failed (truncated/corrupt file?)");
  }

  template <typename T>
  void read_pod(T& v) {
    static_assert(std::is_trivially_copyable_v<T>, "read_pod requires trivially copyable type");
    read_bytes(&v, sizeof(T));
  }

  std::uint8_t read_u8() { std::uint8_t v{}; read_pod(v); return v; }
  std::uint32_t read_u32() { std::uint32_t v{}; read_pod(v); return v; }
  std::uint64_t read_u64() { std::uint64_t v{}; read_pod(v); return v; }
  std::int32_t read_i32() { std::int32_t v{}; read_pod(v); return v; }
  std::int64_t read_i64() { std::int64_t v{}; read_pod(v); return v; }
  double read_f64() { double v{}; read_pod(v); return v; }

  std::string read_string() {
    const std::uint32_t n = read_u32();
    std::string s;
    s.resize(n);
    if (n > 0) read_bytes(s.data(), n);
    return s;
  }

  template <typename T>
  void read_vec_pod(std::vector<T>& v) {
    static_assert(std::is_trivially_copyable_v<T>, "read_vec_pod requires trivially copyable type");
    const std::uint64_t n = read_u64();
    if (n > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
      throw std::runtime_error("BinaryReader: vector size overflow");
    }
    v.resize(static_cast<std::size_t>(n));
    if (n > 0) read_bytes(v.data(), sizeof(T) * static_cast<std::size_t>(n));
  }

private:
  std::istream& is_;
};

inline void require_magic(std::istream& is, std::string_view magic) {
  std::string got;
  got.resize(magic.size());
  is.read(got.data(), static_cast<std::streamsize>(magic.size()));
  if (!is) throw std::runtime_error("BinaryReader: missing magic (truncated/corrupt file?)");
  if (std::string_view(got) != magic) {
    throw std::runtime_error("BinaryReader: magic mismatch (not a PILOTS state file, or wrong version)");
  }
}

inline void write_magic(std::ostream& os, std::string_view magic) {
  os.write(magic.data(), static_cast<std::streamsize>(magic.size()));
  if (!os) throw std::runtime_error("BinaryWriter: failed to write magic");
}

} // namespace pilots::util
