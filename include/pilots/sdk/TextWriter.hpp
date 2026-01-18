#pragma once

#include <filesystem>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>

#include "pilots/util/AtomicFile.hpp"

namespace pilots::sdk {

namespace fs = std::filesystem;

enum class TextWriteMode {
  AtomicRewrite,
  Append,
};

inline std::string text_write_mode_name(TextWriteMode m) {
  switch (m) {
    case TextWriteMode::AtomicRewrite: return "atomic_rewrite";
    case TextWriteMode::Append: return "append";
  }
  return "atomic_rewrite";
}

class TextWriter {
public:
  TextWriter() = default;

  explicit TextWriter(fs::path path,
                      TextWriteMode mode = TextWriteMode::AtomicRewrite,
                      bool dry_run = false)
      : path_(std::move(path)), mode_(mode), dry_run_(dry_run) {}

  const fs::path& path() const { return path_; }
  TextWriteMode mode() const { return mode_; }
  bool dry_run() const { return dry_run_; }

  void write(const std::function<void(std::ostream&)>& fn) const {
    if (dry_run_) return;
    if (mode_ == TextWriteMode::AtomicRewrite) {
      pilots::util::atomic_write_text(path_, fn);
      return;
    }

    std::ofstream ofs(path_, std::ios::app);
    if (!ofs) throw std::runtime_error("TextWriter: failed to open for append: " + path_.string());
    fn(ofs);
    ofs.flush();
    if (!ofs) throw std::runtime_error("TextWriter: write failed: " + path_.string());
  }

  // For append mode, optionally write a header if the file is missing/empty.
  void write_with_header_if_empty(const std::function<void(std::ostream&)>& header_fn,
                                 const std::function<void(std::ostream&)>& body_fn) const {
    if (dry_run_) return;
    if (mode_ == TextWriteMode::AtomicRewrite) {
      pilots::util::atomic_write_text(path_, [&](std::ostream& os) {
        header_fn(os);
        body_fn(os);
      });
      return;
    }

    bool need_header = true;
    std::error_code ec;
    const auto sz = fs::file_size(path_, ec);
    if (!ec && sz > 0) need_header = false;

    std::ofstream ofs(path_, std::ios::app);
    if (!ofs) throw std::runtime_error("TextWriter: failed to open for append: " + path_.string());
    if (need_header) header_fn(ofs);
    body_fn(ofs);
    ofs.flush();
    if (!ofs) throw std::runtime_error("TextWriter: write failed: " + path_.string());
  }

private:
  fs::path path_;
  TextWriteMode mode_ = TextWriteMode::AtomicRewrite;
  bool dry_run_ = false;
};

} // namespace pilots::sdk
