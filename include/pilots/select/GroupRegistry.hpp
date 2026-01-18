#pragma once

#include <algorithm>
#include <cstddef>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pilots/core/Frame.hpp"
#include "pilots/select/Selection.hpp"
#include "pilots/select/SelectionView.hpp"

namespace pilots {

// Platform-level group registry:
// - Supports base selector groups (type/mol/id/region) via Selection.hpp
// - Supports composite groups referencing other groups with &, |, and parentheses
// - Supports explicit selection_mode: static|dynamic (default: static)
// - Caches static groups once (computed on the first frame)
// - Caches dynamic groups per frame (by frame_index) to avoid duplicate computation
class GroupRegistry {
public:
  struct SizeStats {
    std::size_t samples = 0;
    std::size_t min = 0;
    std::size_t max = 0;
    double mean = 0.0;
    double m2 = 0.0; // sum of squared deviations (Welford)

    void add(std::size_t x) {
      if (samples == 0) {
        samples = 1;
        min = max = x;
        mean = static_cast<double>(x);
        m2 = 0.0;
        return;
      }
      ++samples;
      if (x < min) min = x;
      if (x > max) max = x;
      const double dx = static_cast<double>(x) - mean;
      mean += dx / static_cast<double>(samples);
      const double dx2 = static_cast<double>(x) - mean;
      m2 += dx * dx2;
    }

    double variance() const {
      return (samples > 1) ? (m2 / static_cast<double>(samples - 1)) : 0.0;
    }
  };

  struct GroupDef {
    enum class Kind { Selector, Expr } kind = Kind::Selector;
    bool dynamic = false;
    std::string expr; // selector string or boolean expression
  };

  explicit GroupRegistry(std::unordered_map<std::string, std::string> raw_defs)
  : raw_defs_(std::move(raw_defs)) {
    parse_all_();
  }

  bool has(const std::string& name) const {
    return defs_.find(name) != defs_.end();
  }

  bool is_dynamic(const std::string& name) const {
    auto it = defs_.find(name);
    if (it == defs_.end()) return false;
    return it->second.dynamic;
  }

  // Prepare and cache all static groups using the provided frame (typically the first frame).
  void prepare_static(const Frame& frame0) {
    prepared_ = true;
    // Always provide built-in all group.
    static_cache_["all"] = make_all_selection(frame0, "all");
    // Optional override of 'all' via a user selector (must be static selector).
    if (auto it = defs_.find("all"); it != defs_.end()) {
      const GroupDef& def = it->second;
      if (def.dynamic) {
        throw std::runtime_error("GroupRegistry: group 'all' cannot be dynamic");
      }
      if (def.kind != GroupDef::Kind::Selector) {
        throw std::runtime_error("GroupRegistry: group 'all' override must be a selector expression");
      }
      static_cache_["all"] = build_selection_from_selector(frame0, "all", def.expr);
    }

    // Evaluate all other static defs.
    std::vector<std::string> stack;
    for (const auto& [name, def] : defs_) {
      if (name == "all") continue;
      if (!def.dynamic) {
        (void)eval_static_(name, frame0, stack);
      }
    }
  }

  // Get a selection by name for a given frame.
  // - For static groups, returns the cached selection.
  // - For dynamic groups, computes/caches per frame_index.
  const Selection& get(const std::string& name, const Frame& frame, std::size_t frame_index) {
    if (!prepared_) {
      throw std::runtime_error("GroupRegistry: prepare_static() must be called before get()");
    }

    // Built-in all
    if (name == "all") {
      auto it = static_cache_.find("all");
      if (it == static_cache_.end()) {
        throw std::runtime_error("GroupRegistry: internal error: missing built-in 'all'");
      }
      return it->second;
    }

    auto it = defs_.find(name);
    if (it == defs_.end()) {
      throw std::runtime_error("GroupRegistry: unknown group name: '" + name + "'");
    }
    const GroupDef& def = it->second;
    if (!def.dynamic) {
      auto sit = static_cache_.find(name);
      if (sit == static_cache_.end()) {
        throw std::runtime_error("GroupRegistry: static group not prepared: '" + name + "'");
      }
      return sit->second;
    }

    // Dynamic group: per-frame cache
    auto& entry = dynamic_cache_[name];
    if (entry.last_frame_index == frame_index) {
      return entry.sel;
    }
    entry.sel = eval_dynamic_(name, frame);
    entry.last_frame_index = frame_index;
    entry.size_stats.add(entry.sel.idx.size());
    return entry.sel;
  }

  SelectionView get_view(const std::string& name, const Frame& frame, std::size_t frame_index) {
    const Selection& s = get(name, frame, frame_index);
    return view_of(s);
  }

  struct GroupAuditInfo {
    std::string name;
    std::string expr;
    bool is_dynamic = false;
    std::size_t static_size = 0;
    SizeStats size_stats; // samples=0 => never evaluated
  };

  std::vector<GroupAuditInfo> audit() const {
    std::vector<GroupAuditInfo> out;
    out.reserve(defs_.size() + 1);

    // built-in 'all'
    {
      GroupAuditInfo g;
      g.name = "all";
      g.expr = (defs_.find("all") != defs_.end()) ? defs_.at("all").expr : "type:*";
      g.is_dynamic = false;
      if (auto it = static_cache_.find("all"); it != static_cache_.end()) {
        g.static_size = it->second.idx.size();
      }
      out.push_back(std::move(g));
    }

    for (const auto& kv : defs_) {
      const std::string& name = kv.first;
      if (name == "all") continue;
      const GroupDef& def = kv.second;
      GroupAuditInfo g;
      g.name = name;
      g.expr = def.expr;
      g.is_dynamic = def.dynamic;
      if (!def.dynamic) {
        if (auto it = static_cache_.find(name); it != static_cache_.end()) {
          g.static_size = it->second.idx.size();
        }
      } else {
        if (auto it = dynamic_cache_.find(name); it != dynamic_cache_.end()) {
          g.size_stats = it->second.size_stats;
        }
      }
      out.push_back(std::move(g));
    }

    std::sort(out.begin(), out.end(), [](const GroupAuditInfo& a, const GroupAuditInfo& b) {
      return a.name < b.name;
    });
    return out;
  }

private:
  // --- Expression AST ---
  struct Node {
    enum class Kind { Name, And, Or } kind = Kind::Name;
    std::string name; // for Kind::Name
    std::unique_ptr<Node> lhs;
    std::unique_ptr<Node> rhs;
  };

  struct DynEntry {
    Selection sel;
    std::size_t last_frame_index = static_cast<std::size_t>(-1);
    SizeStats size_stats;
  };

  std::unordered_map<std::string, std::string> raw_defs_;
  std::unordered_map<std::string, GroupDef> defs_;
  std::unordered_map<std::string, std::unique_ptr<Node>> ast_;

  bool prepared_ = false;
  std::unordered_map<std::string, Selection> static_cache_;
  std::unordered_map<std::string, DynEntry> dynamic_cache_;

  // --- Parsing helpers ---
  static std::string trim_(const std::string& s) {
    std::size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    std::size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e-1]))) --e;
    return s.substr(b, e - b);
  }

  static bool starts_with_(const std::string& s, const std::string& pref) {
    return s.size() >= pref.size() && s.compare(0, pref.size(), pref) == 0;
  }

  static std::vector<std::string> split_semicolon_(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : s) {
      if (ch == ';') {
        auto t = trim_(cur);
        if (!t.empty()) out.push_back(t);
        cur.clear();
      } else {
        cur.push_back(ch);
      }
    }
    auto t = trim_(cur);
    if (!t.empty()) out.push_back(t);
    return out;
  }

  static std::string join_semicolon_(const std::vector<std::string>& parts, std::size_t from) {
    std::string out;
    for (std::size_t i = from; i < parts.size(); ++i) {
      if (!out.empty()) out += "; ";
      out += parts[i];
    }
    return out;
  }

  // Tokenizer for boolean group expressions: identifiers, &, |, (, )
  struct Tok {
    enum class K { Ident, And, Or, LPar, RPar, End } k = K::End;
    std::string s;
  };

  struct Lexer {
    const std::string& src;
    std::size_t i = 0;
    explicit Lexer(const std::string& s) : src(s) {}

    Tok next() {
      while (i < src.size() && std::isspace(static_cast<unsigned char>(src[i]))) ++i;
      if (i >= src.size()) return Tok{Tok::K::End, ""};
      const char ch = src[i];
      if (ch == '&') { ++i; return Tok{Tok::K::And, "&"}; }
      if (ch == '|') { ++i; return Tok{Tok::K::Or, "|"}; }
      if (ch == '(') { ++i; return Tok{Tok::K::LPar, "("}; }
      if (ch == ')') { ++i; return Tok{Tok::K::RPar, ")"}; }

      // identifier: [A-Za-z0-9_\-]
      std::size_t j = i;
      while (j < src.size()) {
        char c = src[j];
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-') {
          ++j;
        } else {
          break;
        }
      }
      if (j == i) {
        throw std::runtime_error("GroupRegistry: invalid character in group expression: '" + std::string(1, src[i]) + "'");
      }
      Tok t;
      t.k = Tok::K::Ident;
      t.s = src.substr(i, j - i);
      i = j;
      return t;
    }
  };

  struct Parser {
    Lexer lex;
    Tok cur;
    explicit Parser(const std::string& s) : lex(s) { cur = lex.next(); }

    void eat(Tok::K k) {
      if (cur.k != k) {
        throw std::runtime_error("GroupRegistry: parse error in expression");
      }
      cur = lex.next();
    }

    std::unique_ptr<Node> parse_expr() {
      auto n = parse_term();
      while (cur.k == Tok::K::Or) {
        eat(Tok::K::Or);
        auto rhs = parse_term();
        auto parent = std::make_unique<Node>();
        parent->kind = Node::Kind::Or;
        parent->lhs = std::move(n);
        parent->rhs = std::move(rhs);
        n = std::move(parent);
      }
      return n;
    }

    std::unique_ptr<Node> parse_term() {
      auto n = parse_factor();
      while (cur.k == Tok::K::And) {
        eat(Tok::K::And);
        auto rhs = parse_factor();
        auto parent = std::make_unique<Node>();
        parent->kind = Node::Kind::And;
        parent->lhs = std::move(n);
        parent->rhs = std::move(rhs);
        n = std::move(parent);
      }
      return n;
    }

    std::unique_ptr<Node> parse_factor() {
      if (cur.k == Tok::K::Ident) {
        auto n = std::make_unique<Node>();
        n->kind = Node::Kind::Name;
        n->name = cur.s;
        eat(Tok::K::Ident);
        return n;
      }
      if (cur.k == Tok::K::LPar) {
        eat(Tok::K::LPar);
        auto n = parse_expr();
        eat(Tok::K::RPar);
        return n;
      }
      throw std::runtime_error("GroupRegistry: parse error: expected name or '(' ");
    }
  };

  void parse_all_() {
    defs_.clear();
    ast_.clear();

    for (const auto& [name_raw, expr_raw] : raw_defs_) {
      const std::string name = trim_(name_raw);
      if (name.empty()) continue;
      if (name == "all") {
        // allow overriding all via selector
      }

      GroupDef def;
      std::string expr = trim_(expr_raw);
      if (expr.empty()) {
        throw std::runtime_error("GroupRegistry: empty definition for group '" + name + "'");
      }

      // Extract optional selection_mode directive from leading ';'-separated clauses.
      const auto parts = split_semicolon_(expr);
      std::size_t start = 0;
      if (!parts.empty() && starts_with_(parts[0], "selection_mode:")) {
        const std::string v = trim_(parts[0].substr(std::string("selection_mode:").size()));
        if (v == "dynamic") def.dynamic = true;
        else if (v == "static") def.dynamic = false;
        else throw std::runtime_error("GroupRegistry: invalid selection_mode for group '" + name + "': " + v);
        start = 1;
      }
      expr = join_semicolon_(parts, start);
      if (expr.empty()) {
        throw std::runtime_error("GroupRegistry: empty expression after selection_mode for group '" + name + "'");
      }

      // Decide whether this is a base selector or a boolean expression.
      // Heuristic: selector expressions contain ':' (type:..., region:...).
      if (expr.find(':') != std::string::npos) {
        def.kind = GroupDef::Kind::Selector;
        def.expr = expr;
      } else {
        def.kind = GroupDef::Kind::Expr;
        def.expr = expr;
        Parser p(expr);
        auto root = p.parse_expr();
        if (p.cur.k != Tok::K::End) {
          throw std::runtime_error("GroupRegistry: trailing tokens in group expression for '" + name + "'");
        }
        ast_[name] = std::move(root);
      }

      defs_[name] = std::move(def);
    }
  }

  static std::vector<std::size_t> intersect_sorted_(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
    std::vector<std::size_t> out;
    out.reserve(std::min(a.size(), b.size()));
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

  static std::vector<std::size_t> unite_sorted_(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
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

  const Selection& eval_static_(const std::string& name, const Frame& frame0, std::vector<std::string>& stack) {
    if (auto it = static_cache_.find(name); it != static_cache_.end()) {
      return it->second;
    }
    auto dit = defs_.find(name);
    if (dit == defs_.end()) {
      throw std::runtime_error("GroupRegistry: unknown referenced group '" + name + "'");
    }
    const GroupDef& def = dit->second;
    if (def.dynamic) {
      throw std::runtime_error("GroupRegistry: static group depends on dynamic group '" + name + "'");
    }
    // cycle detection
    for (const auto& s : stack) {
      if (s == name) {
        throw std::runtime_error("GroupRegistry: cycle detected in group definitions at '" + name + "'");
      }
    }
    stack.push_back(name);

    Selection sel;
    if (def.kind == GroupDef::Kind::Selector) {
      sel = build_selection_from_selector(frame0, name, def.expr);
    } else {
      auto ait = ast_.find(name);
      if (ait == ast_.end() || !ait->second) {
        throw std::runtime_error("GroupRegistry: missing AST for group '" + name + "'");
      }
      sel = eval_expr_static_(*ait->second, frame0, stack);
      sel.name = name;
    }

    auto [it, ok] = static_cache_.emplace(name, std::move(sel));
    if (!ok) {
      throw std::runtime_error("GroupRegistry: duplicate static cache insert for '" + name + "'");
    }
    stack.pop_back();
    return it->second;
  }

  Selection eval_expr_static_(const Node& node, const Frame& frame0, std::vector<std::string>& stack) {
    if (node.kind == Node::Kind::Name) {
      const Selection& ref = eval_static_(node.name, frame0, stack);
      return ref; // copy
    }
    if (!node.lhs || !node.rhs) {
      throw std::runtime_error("GroupRegistry: malformed AST");
    }
    Selection a = eval_expr_static_(*node.lhs, frame0, stack);
    Selection b = eval_expr_static_(*node.rhs, frame0, stack);
    Selection out;
    if (node.kind == Node::Kind::And) {
      out.idx = intersect_sorted_(a.idx, b.idx);
    } else if (node.kind == Node::Kind::Or) {
      out.idx = unite_sorted_(a.idx, b.idx);
    }
    return out;
  }

  Selection eval_dynamic_(const std::string& name, const Frame& frame) {
    std::vector<std::string> stack;
    return eval_dynamic_with_stack_(name, frame, stack);
  }

  Selection eval_dynamic_with_stack_(const std::string& name, const Frame& frame, std::vector<std::string>& stack) {
    // cycle detection (dynamic expressions can reference other dynamic groups)
    for (const auto& s : stack) {
      if (s == name) {
        throw std::runtime_error("GroupRegistry: cycle detected in dynamic group evaluation at '" + name + "'");
      }
    }
    stack.push_back(name);

    auto dit = defs_.find(name);
    if (dit == defs_.end()) throw std::runtime_error("GroupRegistry: unknown group '" + name + "'");
    const GroupDef& def = dit->second;
    if (!def.dynamic) {
      // static groups should have been cached
      auto it = static_cache_.find(name);
      if (it == static_cache_.end()) throw std::runtime_error("GroupRegistry: static group not prepared: '" + name + "'");
      stack.pop_back();
      return it->second;
    }

    Selection sel;
    if (def.kind == GroupDef::Kind::Selector) {
      sel = build_selection_from_selector(frame, name, def.expr);
    } else {
      auto ait = ast_.find(name);
      if (ait == ast_.end() || !ait->second) {
        throw std::runtime_error("GroupRegistry: missing AST for dynamic group '" + name + "'");
      }
      sel = eval_expr_dynamic_(*ait->second, frame, stack);
      sel.name = name;
    }
    stack.pop_back();
    return sel;
  }

  Selection eval_expr_dynamic_(const Node& node, const Frame& frame, std::vector<std::string>& stack) {
    if (node.kind == Node::Kind::Name) {
      if (node.name == "all") {
        auto sit = static_cache_.find("all");
        if (sit == static_cache_.end()) throw std::runtime_error("GroupRegistry: missing built-in 'all'");
        return sit->second;
      }
      auto it = defs_.find(node.name);
      if (it == defs_.end()) {
        throw std::runtime_error("GroupRegistry: unknown referenced group '" + node.name + "'");
      }
      if (it->second.dynamic) {
        // recurse dynamically
        return eval_dynamic_with_stack_(node.name, frame, stack);
      }
      // static referenced
      auto sit = static_cache_.find(node.name);
      if (sit == static_cache_.end()) {
        throw std::runtime_error("GroupRegistry: referenced static group not prepared: '" + node.name + "'");
      }
      return sit->second;
    }
    if (!node.lhs || !node.rhs) throw std::runtime_error("GroupRegistry: malformed AST");
    Selection a = eval_expr_dynamic_(*node.lhs, frame, stack);
    Selection b = eval_expr_dynamic_(*node.rhs, frame, stack);
    Selection out;
    if (node.kind == Node::Kind::And) out.idx = intersect_sorted_(a.idx, b.idx);
    else if (node.kind == Node::Kind::Or) out.idx = unite_sorted_(a.idx, b.idx);
    return out;
  }
};

} // namespace pilots
