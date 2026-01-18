#pragma once

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

namespace pilots {

inline std::uint64_t fnv1a64_bytes(const void* data, std::size_t n) {
  const std::uint8_t* p = static_cast<const std::uint8_t*>(data);
  std::uint64_t h = 1469598103934665603ull;
  for (std::size_t i = 0; i < n; ++i) {
    h ^= static_cast<std::uint64_t>(p[i]);
    h *= 1099511628211ull;
  }
  return h;
}

inline std::uint64_t fnv1a64_update(std::uint64_t h, const void* data, std::size_t n) {
  const std::uint8_t* p = static_cast<const std::uint8_t*>(data);
  for (std::size_t i = 0; i < n; ++i) {
    h ^= static_cast<std::uint64_t>(p[i]);
    h *= 1099511628211ull;
  }
  return h;
}

inline std::uint64_t fnv1a64_file(const std::string& path) {
  std::ifstream ifs(path, std::ios::binary);
  if (!ifs) {
    throw std::runtime_error("failed to open file for hashing: " + path);
  }
  std::uint64_t h = 1469598103934665603ull;
  char buf[1 << 15];
  while (ifs) {
    ifs.read(buf, static_cast<std::streamsize>(sizeof(buf)));
    const std::streamsize got = ifs.gcount();
    if (got <= 0) break;
    h = fnv1a64_update(h, buf, static_cast<std::size_t>(got));
  }
  return h;
}

} // namespace pilots
