#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pilots {

// Boolean expression over two atom sets A and T.
//
// Supported grammar:
//   expr   ::= term ( '|' term )*
//   term   ::= factor ( '&' factor )*
//   factor ::= '!' factor | '(' expr ')' | 'A' | 'T'
//
// Semantics:
// - 'A' and 'T' are sorted unique vectors of atom indices.
// - Complement is taken w.r.t. the universe [0, natoms).

namespace detail_combine {

inline std::string strip_ws(std::string_view s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
    out.push_back(c);
  }
  return out;
}

inline std::vector<std::size_t> intersect_sorted(const std::vector<std::size_t>& a,
                                                 const std::vector<std::size_t>& b) {
  std::vector<std::size_t> out;
  out.reserve((a.size() < b.size()) ? a.size() : b.size());
  std::size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] == b[j]) {
      out.push_back(a[i]);
      ++i; ++j;
    } else if (a[i] < b[j]) {
      ++i;
    } else {
      ++j;
    }
  }
  return out;
}

inline std::vector<std::size_t> unite_sorted(const std::vector<std::size_t>& a,
                                             const std::vector<std::size_t>& b) {
  std::vector<std::size_t> out;
  out.reserve(a.size() + b.size());
  std::size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] == b[j]) {
      out.push_back(a[i]);
      ++i; ++j;
    } else if (a[i] < b[j]) {
      out.push_back(a[i]);
      ++i;
    } else {
      out.push_back(b[j]);
      ++j;
    }
  }
  while (i < a.size()) out.push_back(a[i++]);
  while (j < b.size()) out.push_back(b[j++]);
  return out;
}

inline std::vector<std::size_t> complement_sorted(const std::vector<std::size_t>& a,
                                                  std::size_t natoms) {
  std::vector<std::size_t> out;
  out.reserve(natoms - ((a.size() < natoms) ? a.size() : natoms));
  std::size_t j = 0;
  for (std::size_t i = 0; i < natoms; ++i) {
    if (j < a.size() && a[j] == i) {
      ++j;
      continue;
    }
    out.push_back(i);
  }
  return out;
}

struct Node {
  enum class Kind { VarA, VarT, And, Or, Not } kind = Kind::VarA;
  std::unique_ptr<Node> lhs;
  std::unique_ptr<Node> rhs;
};

class Parser {
public:
  explicit Parser(std::string s) : s_(std::move(s)) {}

  std::unique_ptr<Node> parse() {
    pos_ = 0;
    auto n = parse_expr_();
    if (pos_ != s_.size()) {
      throw std::runtime_error("CombineExpr: unexpected trailing token at pos " + std::to_string(pos_));
    }
    return n;
  }

private:
  std::string s_;
  std::size_t pos_ = 0;

  char peek_() const { return (pos_ < s_.size()) ? s_[pos_] : '\0'; }
  char get_() { return (pos_ < s_.size()) ? s_[pos_++] : '\0'; }

  std::unique_ptr<Node> parse_expr_() {
    auto left = parse_term_();
    while (peek_() == '|') {
      get_();
      auto right = parse_term_();
      auto n = std::make_unique<Node>();
      n->kind = Node::Kind::Or;
      n->lhs = std::move(left);
      n->rhs = std::move(right);
      left = std::move(n);
    }
    return left;
  }

  std::unique_ptr<Node> parse_term_() {
    auto left = parse_factor_();
    while (peek_() == '&') {
      get_();
      auto right = parse_factor_();
      auto n = std::make_unique<Node>();
      n->kind = Node::Kind::And;
      n->lhs = std::move(left);
      n->rhs = std::move(right);
      left = std::move(n);
    }
    return left;
  }

  std::unique_ptr<Node> parse_factor_() {
    const char c = peek_();
    if (c == '!') {
      get_();
      auto inner = parse_factor_();
      auto n = std::make_unique<Node>();
      n->kind = Node::Kind::Not;
      n->lhs = std::move(inner);
      return n;
    }
    if (c == '(') {
      get_();
      auto inner = parse_expr_();
      if (peek_() != ')') {
        throw std::runtime_error("CombineExpr: expected ')' at pos " + std::to_string(pos_));
      }
      get_();
      return inner;
    }
    if (c == 'A' || c == 'a') {
      get_();
      auto n = std::make_unique<Node>();
      n->kind = Node::Kind::VarA;
      return n;
    }
    if (c == 'T' || c == 't') {
      get_();
      auto n = std::make_unique<Node>();
      n->kind = Node::Kind::VarT;
      return n;
    }
    throw std::runtime_error("CombineExpr: expected one of A,T,!,() at pos " + std::to_string(pos_));
  }
};

inline std::vector<std::size_t> eval(const Node& n,
                                    const std::vector<std::size_t>& A,
                                    const std::vector<std::size_t>& T,
                                    std::size_t natoms) {
  switch (n.kind) {
    case Node::Kind::VarA: return A;
    case Node::Kind::VarT: return T;
    case Node::Kind::And: {
      if (!n.lhs || !n.rhs) throw std::runtime_error("CombineExpr: malformed And");
      auto a = eval(*n.lhs, A, T, natoms);
      auto b = eval(*n.rhs, A, T, natoms);
      return intersect_sorted(a, b);
    }
    case Node::Kind::Or: {
      if (!n.lhs || !n.rhs) throw std::runtime_error("CombineExpr: malformed Or");
      auto a = eval(*n.lhs, A, T, natoms);
      auto b = eval(*n.rhs, A, T, natoms);
      return unite_sorted(a, b);
    }
    case Node::Kind::Not: {
      if (!n.lhs) throw std::runtime_error("CombineExpr: malformed Not");
      auto a = eval(*n.lhs, A, T, natoms);
      return complement_sorted(a, natoms);
    }
  }
  return A;
}

} // namespace detail_combine

inline std::vector<std::size_t> combine_sets(const std::vector<std::size_t>& A,
                                             const std::vector<std::size_t>& T,
                                             std::size_t natoms,
                                             const std::string& expr_raw) {
  const std::string expr = detail_combine::strip_ws(expr_raw);
  const std::string use = expr.empty() ? std::string("A&T") : expr;
  detail_combine::Parser p(use);
  auto ast = p.parse();
  return detail_combine::eval(*ast, A, T, natoms);
}

} // namespace pilots
