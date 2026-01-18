#pragma once

#include <algorithm>
#include <cstddef>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pilots {

// Minimal INI parser:
// - Sections: [section.name]
// - Key: key = value
// - Comments: lines starting with '#' or ';'
// - Values: raw strings; surrounding quotes (single/double) are stripped.

class IniConfig {
public:
  explicit IniConfig(const std::filesystem::path& file) : file_(file) {
    parse_();
  }

  const std::filesystem::path& file_path() const { return file_; }
  std::filesystem::path base_dir() const { return file_.parent_path(); }

  bool has_section(const std::string& section) const {
    return data_.find(section) != data_.end();
  }

  // Deterministic list of section names (sorted).
  std::vector<std::string> section_names() const {
    std::vector<std::string> out;
    out.reserve(data_.size());
    for (const auto& kv : data_) out.push_back(kv.first);
    std::sort(out.begin(), out.end());
    return out;
  }

  bool has_key(const std::string& section, const std::string& key) const {
    auto it = data_.find(section);
    if (it == data_.end()) return false;
    return it->second.find(key) != it->second.end();
  }

  std::string get_string(const std::string& section, const std::string& key,
                         const std::optional<std::string>& def = std::nullopt) const {
    if (auto v = get_raw_(section, key)) {
      return *v;
    }
    if (def) return *def;
    throw std::runtime_error(err_prefix_() + "missing required key '" + key + "' in section [" + section + "]");
  }

  std::int64_t get_int64(const std::string& section, const std::string& key,
                         const std::optional<std::int64_t>& def = std::nullopt) const {
    std::string s = get_string(section, key, def ? std::optional<std::string>(std::to_string(*def)) : std::nullopt);
    try {
      std::size_t pos = 0;
      long long v = std::stoll(s, &pos);
      if (pos != s.size()) throw std::invalid_argument("trailing chars");
      return static_cast<std::int64_t>(v);
    } catch (...) {
      throw std::runtime_error(err_prefix_() + "failed to parse int64 for " + section + "." + key + " from value: '" + s + "'");
    }
  }

  std::size_t get_size(const std::string& section, const std::string& key,
                       const std::optional<std::size_t>& def = std::nullopt) const {
    std::string s = get_string(section, key, def ? std::optional<std::string>(std::to_string(*def)) : std::nullopt);
    try {
      std::size_t pos = 0;
      unsigned long long v = std::stoull(s, &pos);
      if (pos != s.size()) throw std::invalid_argument("trailing chars");
      return static_cast<std::size_t>(v);
    } catch (...) {
      throw std::runtime_error(err_prefix_() + "failed to parse size for " + section + "." + key + " from value: '" + s + "'");
    }
  }

  double get_double(const std::string& section, const std::string& key,
                    const std::optional<double>& def = std::nullopt) const {
    std::string s = get_string(section, key, def ? std::optional<std::string>(to_string_prec_(*def)) : std::nullopt);
    try {
      std::size_t pos = 0;
      double v = std::stod(s, &pos);
      if (pos != s.size()) throw std::invalid_argument("trailing chars");
      return v;
    } catch (...) {
      throw std::runtime_error(err_prefix_() + "failed to parse double for " + section + "." + key + " from value: '" + s + "'");
    }
  }

  bool get_bool(const std::string& section, const std::string& key,
                const std::optional<bool>& def = std::nullopt) const {
    auto s = get_string(section, key, def ? std::optional<std::string>(*def ? "true" : "false") : std::nullopt);
    for (auto& c : s) c = static_cast<char>(::tolower(c));
    if (s == "1" || s == "true" || s == "yes" || s == "on") return true;
    if (s == "0" || s == "false" || s == "no" || s == "off") return false;
    throw std::runtime_error(err_prefix_() + "failed to parse bool for " + section + "." + key + " from value: '" + s + "'");
  }

  // Parse comma-separated list. Whitespace around items is trimmed.
  std::vector<std::string> get_list(const std::string& section, const std::string& key,
                                    const std::optional<std::string>& def = std::nullopt) const {
    std::string s = get_string(section, key, def);
    std::vector<std::string> out;
    split_csv_(s, out);
    return out;
  }

  // Return the raw section map (useful for iterating over [groups])
  const std::unordered_map<std::string, std::string>& section(const std::string& section) const {
    auto it = data_.find(section);
    if (it == data_.end()) {
      throw std::runtime_error(err_prefix_() + "missing section [" + section + "]");
    }
    return it->second;
  }

private:
  std::filesystem::path file_;
  std::unordered_map<std::string, std::unordered_map<std::string, std::string>> data_;

  static std::string trim_(std::string s) {
    auto is_ws = [](unsigned char ch) { return ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n'; };
    std::size_t b = 0;
    while (b < s.size() && is_ws(static_cast<unsigned char>(s[b]))) ++b;
    std::size_t e = s.size();
    while (e > b && is_ws(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
  }

  static std::string strip_quotes_(std::string s) {
    if (s.size() >= 2) {
      const char a = s.front();
      const char b = s.back();
      if ((a == '"' && b == '"') || (a == '\'' && b == '\'')) {
        return s.substr(1, s.size() - 2);
      }
    }
    return s;
  }

  static void split_csv_(const std::string& s, std::vector<std::string>& out) {
    out.clear();
    std::string cur;
    for (char ch : s) {
      if (ch == ',') {
        auto t = trim_(cur);
        if (!t.empty()) out.push_back(t);
        cur.clear();
      } else {
        cur.push_back(ch);
      }
    }
    auto t = trim_(cur);
    if (!t.empty()) out.push_back(t);
  }

  static std::string to_string_prec_(double x) {
    // Keep stable precision for config-derived defaults.
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", x);
    return std::string(buf);
  }

  std::string err_prefix_() const {
    return std::string("IniConfig[") + file_.string() + "]: ";
  }

  std::optional<std::string> get_raw_(const std::string& section, const std::string& key) const {
    auto it = data_.find(section);
    if (it == data_.end()) return std::nullopt;
    auto it2 = it->second.find(key);
    if (it2 == it->second.end()) return std::nullopt;
    return it2->second;
  }

  void parse_() {
    std::ifstream ifs(file_);
    if (!ifs) {
      throw std::runtime_error(err_prefix_() + "failed to open config");
    }

    std::string section = "";
    std::string line;
    std::size_t lineno = 0;

    while (std::getline(ifs, line)) {
      ++lineno;
      std::string s = trim_(line);
      if (s.empty()) continue;
      if (s[0] == '#' || s[0] == ';') continue;

      if (s.front() == '[' && s.back() == ']') {
        section = trim_(s.substr(1, s.size() - 2));
        if (section.empty()) {
          throw std::runtime_error(err_prefix_() + "empty section header at line " + std::to_string(lineno));
        }
        (void)data_[section];
        continue;
      }

      auto eq = s.find('=');
      if (eq == std::string::npos) {
        throw std::runtime_error(err_prefix_() + "expected key=value at line " + std::to_string(lineno) + ": " + s);
      }

      std::string key = trim_(s.substr(0, eq));
      std::string val = trim_(s.substr(eq + 1));
      if (key.empty()) {
        throw std::runtime_error(err_prefix_() + "empty key at line " + std::to_string(lineno));
      }
      val = strip_quotes_(val);

      if (section.empty()) {
        throw std::runtime_error(err_prefix_() + "key outside any section at line " + std::to_string(lineno) + ": " + key);
      }

      data_[section][key] = val;
    }
  }
};

} // namespace pilots
