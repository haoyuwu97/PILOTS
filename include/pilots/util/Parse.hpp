#pragma once

#include <charconv>
#include <cstddef>
#include <cstdint>
#include <system_error>
#include <string_view>
#include <vector>

namespace pilots {

inline void split_ws(std::string_view s, std::vector<std::string_view>& out) {
  out.clear();
  std::size_t i = 0;
  const std::size_t n = s.size();
  while (i < n) {
    while (i < n && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r' || s[i] == '\n')) ++i;
    if (i >= n) break;
    std::size_t j = i;
    while (j < n && !(s[j] == ' ' || s[j] == '\t' || s[j] == '\r' || s[j] == '\n')) ++j;
    out.emplace_back(s.substr(i, j - i));
    i = j;
  }
}

inline bool next_token(const char*& p, const char* end, std::string_view& tok) {
  while (p < end && (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')) ++p;
  if (p >= end) {
    tok = std::string_view{};
    return false;
  }
  const char* start = p;
  while (p < end && !(*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')) ++p;
  tok = std::string_view(start, static_cast<std::size_t>(p - start));
  return true;
}

template <typename IntT>
inline bool parse_int(std::string_view tok, IntT& value) {
  const char* b = tok.data();
  const char* e = tok.data() + tok.size();
  auto res = std::from_chars(b, e, value);
  return res.ec == std::errc{};
}

inline bool parse_double(std::string_view tok, double& value) {
  const char* b = tok.data();
  const char* e = tok.data() + tok.size();
  auto res = std::from_chars(b, e, value);
  return res.ec == std::errc{};
}

} // namespace pilots
