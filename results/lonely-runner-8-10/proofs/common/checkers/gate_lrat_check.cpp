// gate_lrat_check: independent checker for RUP-only text LRAT proofs (any number of vars).
// Gate-SAT agent, 2026-10-21.  Written from the LRAT format definition, not from Astra's code.
//
// For each addition "id l1..ln 0 h1..hm 0": assume every li false, then walk the hints in
// order; each hint clause must be satisfied (skipped, no inference), unit under the current
// assignment (its last literal is set true), or falsified (conflict, step accepted).  A step
// with no conflict is rejected.  Negative hints (RAT) are rejected.  Deletions "id d ids 0"
// remove clauses (sound: only shrinks the database).  Clause ids must strictly increase.
// Exit 0 and prints "<steps> <max_live_clauses>" only if the empty clause was derived.
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

static void die(const char* m, long long line) {
  std::fprintf(stderr, "REJECT line %lld: %s\n", line, m);
  std::exit(1);
}

struct Reader {
  FILE* f; char buf[1 << 16]; size_t n = 0, i = 0; bool eof = false;
  explicit Reader(FILE* f_) : f(f_) {}
  int get() {
    if (i == n) { if (eof) return -1; n = std::fread(buf, 1, sizeof buf, f); i = 0; if (!n) { eof = true; return -1; } }
    return (unsigned char)buf[i++];
  }
};

int main(int argc, char** argv) {
  if (argc != 3) { std::fprintf(stderr, "usage: gate_lrat_check f.cnf f.lrat\n"); return 2; }
  FILE* fc = std::fopen(argv[1], "rb"); FILE* fp = std::fopen(argv[2], "rb");
  if (!fc || !fp) { std::fprintf(stderr, "cannot open\n"); return 2; }
  std::unordered_map<long long, std::vector<int>> db;
  long long nv = -1, nc = -1;
  // ---- CNF (comment lines start with 'c' in column 0)
  {
    std::vector<int> cur; long long cnt = 0; char line[1 << 16]; bool hdr = false;
    std::string acc;
    while (std::fgets(line, sizeof line, fc)) {
      acc += line;
      if (acc.back() != '\n' && !std::feof(fc)) continue;
      const char* q = acc.c_str();
      while (*q && std::isspace((unsigned char)*q)) ++q;
      if (*q == 'c' || *q == 0) { acc.clear(); continue; }
      if (*q == 'p') {
        if (std::sscanf(q, "p cnf %lld %lld", &nv, &nc) != 2 || hdr) die("bad header", 0);
        hdr = true; acc.clear(); continue;
      }
      if (!hdr) die("clause before header", 0);
      char* e;
      while (true) {
        while (*q && std::isspace((unsigned char)*q)) ++q;
        if (!*q) break;
        long long l = std::strtoll(q, &e, 10);
        if (e == q) die("bad cnf token", 0);
        q = e;
        if (l == 0) { db[++cnt] = cur; cur.clear(); }
        else { if (std::llabs(l) > nv) die("cnf literal out of range", 0); cur.push_back((int)l); }
      }
      acc.clear();
    }
    if (!hdr || !cur.empty() || cnt != nc) die("cnf clause count mismatch", 0);
  }
  // ---- proof
  std::vector<signed char> val(nv + 1, 0);  // +1 true, -1 false, by variable
  std::vector<int> trail;
  auto litval = [&](int l) -> int { int v = val[std::abs(l)]; return l > 0 ? v : -v; };
  long long last = nc, steps = 0, line = 0; size_t live = db.size(), maxlive = live;
  bool empty = false;
  Reader r(fp);
  std::vector<long long> nums; std::string tok;
  int ch = 0;
  while (true) {
    // read one line into tokens
    nums.clear(); tok.clear(); bool isdel = false, any = false;
    while ((ch = r.get()) != -1 && ch != '\n') {
      if (ch == ' ' || ch == '\t' || ch == '\r') {
        if (!tok.empty()) { if (tok == "d") isdel = true; else nums.push_back(std::atoll(tok.c_str())); tok.clear(); }
      } else { tok.push_back((char)ch); any = true; }
    }
    if (!tok.empty()) { if (tok == "d") isdel = true; else nums.push_back(std::atoll(tok.c_str())); }
    ++line;
    if (any && !nums.empty()) {
      if (isdel) {
        for (size_t j = 1; j < nums.size() && nums[j]; ++j) { if (db.erase(nums[j])) --live; }
      } else {
        long long id = nums[0];
        if (id <= last) die("non-increasing clause id", line);
        size_t j = 1; std::vector<int> cl;
        for (; j < nums.size() && nums[j] != 0; ++j) {
          if (std::llabs(nums[j]) > nv) die("literal out of range", line);
          cl.push_back((int)nums[j]);
        }
        if (j >= nums.size()) die("missing clause terminator", line);
        ++j;
        trail.clear(); bool conflict = false;
        for (int l : cl) {
          int v = litval(l);
          if (v == -1) continue;
          if (v == 1) { conflict = true; break; }  // tautology (l and -l in clause)
          val[std::abs(l)] = l > 0 ? -1 : 1; trail.push_back(std::abs(l));
        }
        bool ended = false, hint = false;
        for (; j < nums.size(); ++j) {
          long long h = nums[j];
          if (h == 0) { ended = true; break; }
          if (h < 0) die("RAT hint not supported", line);
          hint = true;
          if (conflict) continue;
          auto it = db.find(h);
          if (it == db.end()) die("unknown hint id", line);
          int unit = 0, nun = 0; bool sat = false;
          for (int l : it->second) {
            int v = litval(l);
            if (v == 1) { sat = true; break; }
            if (v == 0) { if (unit != l) { ++nun; unit = l; } }
          }
          if (sat) continue;
          if (nun == 0) { conflict = true; continue; }
          if (nun > 1) die("hint is not unit", line);
          val[std::abs(unit)] = unit > 0 ? 1 : -1; trail.push_back(std::abs(unit));
        }
        for (int v : trail) val[v] = 0;
        if (!ended) die("missing hint terminator", line);
        if (!conflict) die("no RUP conflict", line);
        if (!hint && !cl.empty()) die("no hints", line);
        db[id] = cl; ++live; if (live > maxlive) maxlive = live;
        last = id; ++steps;
        if (cl.empty()) { empty = true; break; }
      }
    }
    if (ch == -1) break;
  }
  if (!empty) { std::fprintf(stderr, "REJECT: no empty clause derived\n"); return 1; }
  std::printf("%lld %zu\n", steps, maxlive);
  return 0;
}
