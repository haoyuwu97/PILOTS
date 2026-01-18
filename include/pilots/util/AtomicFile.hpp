#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace pilots::util {
namespace fs = std::filesystem;

inline fs::path make_tmp_path(const fs::path& out_path) {
  fs::path tmp = out_path;
  tmp += ".tmp";
  return tmp;
}

// Write to a temporary file in the same directory and then rename over the target.
// This makes outputs robust against partial writes (e.g. kill -9 / node preemption).
inline void atomic_rename_over(const fs::path& tmp_path, const fs::path& out_path) {
  std::error_code ec;
  fs::rename(tmp_path, out_path, ec);
  if (!ec) return;

  // Some platforms/filesystems don't overwrite existing paths on rename.
  fs::remove(out_path, ec);
  ec.clear();
  fs::rename(tmp_path, out_path, ec);
  if (ec) {
    throw std::runtime_error("atomic rename failed: '" + tmp_path.string() + "' -> '" + out_path.string() + "' (" + ec.message() + ")");
  }
}

// Convenience: atomic write for a text file.
template <typename WriteFn>
inline void atomic_write_text(const fs::path& out_path, WriteFn&& fn) {
  const fs::path tmp = make_tmp_path(out_path);
  {
    std::ofstream ofs(tmp);
    if (!ofs) throw std::runtime_error("failed to open temp file for atomic write: " + tmp.string());
    fn(ofs);
    ofs.flush();
    if (!ofs) throw std::runtime_error("failed while writing temp file: " + tmp.string());
  }
  atomic_rename_over(tmp, out_path);
}

// Convenience: atomic write for a binary file.
template <typename WriteFn>
inline void atomic_write_binary(const fs::path& out_path, WriteFn&& fn) {
  const fs::path tmp = make_tmp_path(out_path);
  {
    std::ofstream ofs(tmp, std::ios::binary);
    if (!ofs) throw std::runtime_error("failed to open temp file for atomic write: " + tmp.string());
    fn(ofs);
    ofs.flush();
    if (!ofs) throw std::runtime_error("failed while writing temp file: " + tmp.string());
  }
  atomic_rename_over(tmp, out_path);
}

} // namespace pilots::util
