// level1: C++ port of proofs/common/level1.py (same search, same leaf order; see the proof
// there).  Prints one leaf per line "I c1 c2 .. | E e1 e2 ..", then "# nodes dead leaves".
// usage: level1 K P
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
typedef unsigned long long u64;
static int K, P, H, W;
static std::vector<std::vector<u64>> cov;   // cov[c] bitset over T1
static std::vector<u64> FULL;
static long long nodes = 0, dead = 0, nleaves = 0;
static std::vector<int> I;
static std::vector<char> inE;
static int pc(const u64* a, const u64* b) { int s = 0; for (int w = 0; w < W; ++w) s += __builtin_popcountll(a[w] & b[w]); return s; }

static void dfs(std::vector<u64>& U) {
  ++nodes;
  bool empty = true; int nu = 0;
  for (int w = 0; w < W; ++w) { if (U[w]) empty = false; nu += __builtin_popcountll(U[w]); }
  if (empty) {
    ++nleaves;
    printf("I");
    for (int c : I) printf(" %d", c);
    printf(" | E");
    for (int c = 2; c <= H; ++c) if (inE[c]) printf(" %d", c);
    printf("\n");
    return;
  }
  if ((int)I.size() == K - 1) { ++dead; return; }
  std::vector<char> inI(H + 1, 0);
  for (int c : I) inI[c] = 1;
  int r = K - 1 - (int)I.size(), mx = 0;
  for (int c = 2; c <= H; ++c) if (!inE[c] && !inI[c]) mx = std::max(mx, pc(cov[c].data(), U.data()));
  if ((long long)r * mx < nu) { ++dead; return; }
  // pick the uncovered time (lowest bit order) with the fewest candidates; first minimum wins
  std::vector<int> best; bool have = false;
  for (int w = 0; w < W && !(have && best.empty()); ++w) {
    u64 m = U[w];
    while (m) {
      u64 low = m & (~m + 1); m ^= low;
      std::vector<int> cand;
      for (int c = 2; c <= H; ++c) if (!inE[c] && !inI[c] && (cov[c][w] & low)) cand.push_back(c);
      if (!have || cand.size() < best.size()) { best = cand; have = true; if (best.empty()) break; }
    }
  }
  if (best.empty()) { ++dead; return; }
  std::vector<int> added;
  for (int c : best) {
    std::vector<u64> U2(W);
    for (int w = 0; w < W; ++w) U2[w] = U[w] & ~cov[c][w];
    I.push_back(c);
    dfs(U2);
    I.pop_back();
    inE[c] = 1; added.push_back(c);
  }
  for (int c : added) inE[c] = 0;
}

int main(int argc, char** argv) {
  if (argc != 3) { fprintf(stderr, "usage: level1 K P\n"); return 2; }
  K = atoi(argv[1]); P = atoi(argv[2]); H = (P - 1) / 2;
  std::vector<int> T1;
  for (int s = 1; s <= H; ++s) if ((long long)(K + 1) * s >= P) T1.push_back(s);
  W = ((int)T1.size() + 63) / 64;
  cov.assign(H + 1, std::vector<u64>(W, 0));
  for (int c = 2; c <= H; ++c)
    for (size_t i = 0; i < T1.size(); ++i) {
      long long r = (long long)T1[i] * c % P; if (P - r < r) r = P - r;
      if ((long long)(K + 1) * r < P) cov[c][i / 64] |= 1ULL << (i % 64);
    }
  std::vector<u64> U(W, 0);
  for (size_t i = 0; i < T1.size(); ++i) U[i / 64] |= 1ULL << (i % 64);
  inE.assign(H + 1, 0);
  dfs(U);
  printf("# %lld %lld %lld\n", nodes, dead, nleaves);
  return 0;
}
