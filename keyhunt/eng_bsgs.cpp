/*
 * Experimental bounded/segmented BSGS benchmark.
 *
 * This intentionally lives outside keyhunt.cpp: the production search modes
 * remain unchanged while the memory/time trade-off can be measured on a
 * finite, user-supplied interval and a known test public key.
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <malloc.h>
#include <fcntl.h>
#include <fstream>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <sys/mman.h>
#include <sys/stat.h>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Point.h"
#include "secp256k1/Int.h"

namespace {

struct BabyRow {
  std::array<unsigned char, 33> key;
  uint64_t j;
};

Secp256K1 curve;
Point target;
Point search_target;
Int range_start;
uint64_t range_size = 0;
uint64_t segment_size = 100000;
uint64_t requested_baby_rows = 0;
bool glv_orbit_mode = false;
unsigned worker_count = 1;
std::atomic<uint64_t> segments_done{0};
std::atomic<uint64_t> baby_rows{0};
std::atomic<uint64_t> raw_baby_rows{0};
std::atomic<uint64_t> giant_steps{0};
std::atomic<uint64_t> giant_steps_scheduled{0};
std::atomic<uint64_t> giant_steps_planned{0};
std::atomic<uint64_t> invariant_hits{0};
// Build/search accounting split - comparable across table modes and
// kill/resume cycles. Names state what each counter actually counts:
//   merged_rows       - input rows folded into an existing row while
//                       BUILDING the table (orbit collapse: equal canonical
//                       scalar or equal K3 key). Expected by construction,
//                       not a search event: 2 per head on a 3-row expanded
//                       file, 0 on a canonical file and on the continuous
//                       window. Session-pure, never restored by --resume.
//   insert_collisions - rows whose K3 key equals an existing row's key but
//                       whose stored scalar differs - a genuine cross-orbit
//                       collision. Must stay 0 on full C6 orbits (equal K3
//                       key implies exactly one orbit); nonzero means a
//                       broken invariant: the row is dropped and stderr
//                       warns. Detectable only in --canonical-k3 mode,
//                       where the stored scalar is the orbit minimum.
//   search_hits       - giant-scan bucket matches, one per matched giant
//                       position regardless of bucket row count (restored
//                       by --resume, cumulative across sessions).
//   invariant_hits stays the legacy row-level counter.
std::atomic<uint64_t> merged_rows{0};
std::atomic<uint64_t> insert_collisions{0};
std::atomic<uint64_t> search_hits{0};
// Giant steps restored from the checkpoint by --resume; subtracted from
// the cumulative counter so steps/s rates stay per-session honest.
std::atomic<uint64_t> restored_giant_steps{0};
// Phase-granular RSS timeline (kB). External progress samples every 60 s
// miss malloc transients (glibc mmap blocks are freed and returned to the
// OS before the first sample), so the build phases sample VmRSS at their
// boundaries and VmHWM (high-water mark) records the true peak so far.
// current_rss_kb() is defined next to the progress ticker; VmHWM here.
static long peak_rss_kb() {
  FILE *f = std::fopen("/proc/self/status", "r");
  if (!f) return -1;
  char line[256];
  long hwm = -1;
  while (std::fgets(line, sizeof(line), f)) {
    if (std::sscanf(line, "VmHWM: %ld kB", &hwm) == 1) break;
  }
  std::fclose(f);
  return hwm;
}
static long rss_after_load_kb = 0, hwm_after_load_kb = 0;
static long rss_after_prepass_kb = 0, hwm_after_prepass_kb = 0;
static long rss_after_build_kb = 0, hwm_after_build_kb = 0;
std::atomic<uint64_t> exact_candidate_checks{0};
std::atomic<uint64_t> orbit_factor_checks{0};
std::atomic<uint64_t> signature_routes{0};
std::atomic<uint64_t> signature_anomalies{0};
std::atomic<uint64_t> orbit_duplicate_rows{0};
std::atomic<uint64_t> signature_pair_collapses{0};
std::atomic<uint64_t> signature_pair_mismatches{0};
std::atomic<uint64_t> partner_recoveries{0};
// Honest EC accounting: every scalar multiplication (ComputePublicKey) and
// every real point addition (add_safe with both operands non-zero) counts
// here, so builds and search paths compare on the same unit.
std::atomic<uint64_t> ec_scalar_mults{0};
std::atomic<uint64_t> ec_point_adds{0};
// EC work of the exact candidate checks (the reach-on-a-bucket-match path).
// Counted additionally to the global counters so the printout can split
// build / giant-search / hit EC-ops.
std::atomic<uint64_t> hit_ec_scalar_mults{0};
// Microseconds of the program run at the moment of the successful recovery
// (steady_clock); grouped with the counters because the printout of the
// split accounting references it.
std::atomic<uint64_t> found_run_us{0};
std::atomic<uint64_t> lambda_calculations{0};
std::atomic<uint64_t> priority_slices{0};
std::atomic<uint64_t> fallback_slices{0};
std::atomic<bool> found{false};
std::mutex found_mutex;
std::vector<BabyRow> shared_table;
uint64_t shared_baby_rows = 0;
Int found_scalar;
Int found_shift;
Int found_shifted_scalar;
uint64_t found_segment = 0;
uint64_t found_local_offset = 0;
Int glv_lambda;
Int glv_lambda2;
Int active_shift;
uint64_t active_range_size = 0;
bool enforce_cover_range = false;
Int cover_start;
Int cover_end;
// --checkpoint-file PATH / --checkpoint-steps N / --resume: atomic progress
// checkpoint for long --cover-start/--cover-size runs. Granularity is the
// external cover chunk (one full giant pass over range_size), written after
// every N completed chunks (default 1: every chunk) plus at the end of the
// run. --resume continues from source_start and restores the search-phase
// counters, so the printed totals stay cumulative across kills/sessions.
static std::string checkpoint_path;
static uint64_t checkpoint_steps = 1;
static bool checkpoint_resume = false;
bool lambda_priority = false;
uint64_t lambda_sector_count = 16;
bool phase_signature = false;
// --partner-table: exact-verify out-of-interval GLV orbit branches. A K3 hit
// never fires spuriously (roots of X^3-c are exactly one orbit's x-triple),
// so any branch candidate that passes the EC check IS the target scalar even
// when it lies outside the searched interval (an orbit partner of a row).
bool partner_table_mode = false;
// --canonical-k3: store one canonical row per C6 orbit. The row scalar is the
// orbit minimum c = min_r min(lambda^r*k mod n, n - lambda^r*k mod n); K3=x^3
// and the (d01,d12) signature depend on x-coordinates only, so they are
// unchanged on the canonical representative (up to the phi-rotation the
// router already absorbs). Rows of one orbit collapse at build time:
// one row and one signature per orbit, zero runtime dedup logic.
bool canonical_k3_mode = false;
// --table-scalars FILE: build the baby table from an explicit scalar list
// instead of the continuous window [0, m). This is the research mode where
// K3 buckets genuinely contain orbit partners - the full-range representative
// table anticipated by the routing-economy measurements.
std::string table_scalars_path;
Int beta;
Int beta_g;
Int beta2;

// Parallel to shared_table (same order): (d01, d12) of the row point, used
// by --phase-signature to route a K3 bucket match to one GLV branch pair.
std::vector<std::array<unsigned char, 64>> shared_signatures;
// Parallel to shared_table (GLV mode only): the row scalar as a full 32-byte
// value mod n. Plain mode stores j < m; --canonical-k3 stores the orbit
// minimum, which can lie anywhere in [0, n/2].
std::vector<std::array<unsigned char, 32>> shared_scalars;

Point add_safe(Point &a, Point &b) {
  if (a.isZero()) return Point(b);
  if (b.isZero()) return Point(a);
  ec_point_adds.fetch_add(1, std::memory_order_relaxed);
  if (a.x.IsEqual(&b.x)) {
    if (a.y.IsEqual(&b.y)) return curve.DoubleDirect(a);
    Int y_sum(&a.y);
    y_sum.Add(&b.y);
    y_sum.Mod(&curve.P);
    if (y_sum.IsZero()) {
      Point infinity;
      infinity.Clear();
      return infinity;
    }
  }
  return curve.AddDirect(a, b);
}

// Counted scalar multiplication wrapper: every build row, anchor/step/shift
// point and exact candidate check goes through here, so the printed
// "EC scalar multiplications" covers the whole program run.
Point ec_mult(Int &k) {
  ec_scalar_mults.fetch_add(1, std::memory_order_relaxed);
  return curve.ComputePublicKey(&k);
}

// Exact candidate check wrapper: the same EC work counted additionally into
// the hit-phase counter (giant = search − exact checks in the printout).
Point hit_ec_mult(Int &k) {
  hit_ec_scalar_mults.fetch_add(1, std::memory_order_relaxed);
  return ec_mult(k);
}

std::chrono::steady_clock::time_point program_started;

uint64_t run_us() {
  return std::chrono::duration_cast<std::chrono::microseconds>(
      std::chrono::steady_clock::now() - program_started).count();
}

std::array<unsigned char, 33> point_key(Point &p) {
  std::array<unsigned char, 33> raw{};
  if (p.isZero()) {
    return raw;
  }
  raw[0] = p.y.IsEven() ? 0x02 : 0x03;
  alignas(8) unsigned char x_bytes[32];
  p.x.Get32Bytes(x_bytes);
  std::memcpy(raw.data() + 1, x_bytes, 32);
  return raw;
}

std::array<unsigned char, 33> invariant_key_from_x(Int &x) {
  std::array<unsigned char, 33> raw{};
  Int cube;
  cube.ModCube(&x);
  alignas(8) unsigned char cube_bytes[32];
  cube.Get32Bytes(cube_bytes);
  std::memcpy(raw.data() + 1, cube_bytes, 32);
  return raw;
}

std::array<unsigned char, 33> invariant_key(Point &p) {
  if (p.isZero()) return std::array<unsigned char, 33>{};
  return invariant_key_from_x(p.x);
}

Int lambda_value_for_x(Int x) {
  if (x.IsZero()) {
    Int zero;
    zero.SetInt32(0);
    return zero;
  }
  Int beta_x;
  beta_x.ModMulK1(&beta, &x);
  Int left_g;
  left_g.ModSub(&curve.G.x, &beta_g);
  Int left_x;
  left_x.ModSub(&x, &beta_x);
  Int numerator;
  numerator.ModMulK1(&left_g, &left_x);

  Int right_g;
  right_g.ModSub(&curve.G.x, &beta_x);
  Int right_x;
  right_x.ModSub(&x, &beta_g);
  Int denominator;
  denominator.ModMulK1(&right_g, &right_x);
  if (denominator.IsZero()) {
    Int zero;
    zero.SetInt32(0);
    return zero;
  }
  denominator.ModInv();
  numerator.ModMulK1(&denominator);
  return numerator;
}

std::array<Int, 3> lambda_triangle_values_from_x(Int &x0) {
  Int x1;
  x1.ModMulK1(&beta, &x0);
  Int x2;
  x2.ModMulK1(&beta, &x1);
  return {lambda_value_for_x(x0), lambda_value_for_x(x1),
          lambda_value_for_x(x2)};
}

std::array<Int, 3> lambda_triangle_values(Point &point) {
  if (point.isZero()) {
    Int zero;
    zero.SetInt32(0);
    return lambda_triangle_values_from_x(zero);
  }
  return lambda_triangle_values_from_x(point.x);
}

// Phase signature (d01, d12) of the ordered Lambda triple of the point with
// x-coordinate x:
//   d01 = Lambda(phi(p)) - Lambda(p), d12 = Lambda(phi^2(p)) - Lambda(phi(p)).
// It is invariant under p -> -p and cyclically permuted by phi, so a bucket
// match on K3=x^3 can be routed to one (+/- lambda^r) branch pair.
std::array<unsigned char, 64> phase_signature_from_x(Int &x0) {
  std::array<unsigned char, 64> raw{};
  std::array<Int, 3> values = lambda_triangle_values_from_x(x0);
  Int d01;
  d01.ModSub(&values[1], &values[0]);
  Int d12;
  d12.ModSub(&values[2], &values[1]);
  alignas(8) unsigned char tmp[32];
  d01.Get32Bytes(tmp);
  std::memcpy(raw.data(), tmp, 32);
  d12.Get32Bytes(tmp);
  std::memcpy(raw.data() + 32, tmp, 32);
  return raw;
}

std::array<unsigned char, 64> phase_signature_of(Point &p) {
  if (p.isZero()) return std::array<unsigned char, 64>{};
  return phase_signature_from_x(p.x);
}

// Canonical representative of the C6 orbit of the scalar s (mod n):
//   c = min over r in {0,1,2} of min(lambda^r*s mod n, n - lambda^r*s mod n).
// Sets r_out to the rotation with c == +/-lambda^{r_out}*s (mod n); the sign
// itself is irrelevant because K3 = x^3 and the (d01,d12) signature depend on
// x-coordinates only, and negation does not change x. The map s -> c is
// constant on orbits and distinct across orbits, so "equal canonical scalar"
// is exactly "same orbit".
Int canonical_orbit_scalar(Int &s, int &r_out) {
  Int reps[3];
  reps[0].Set(&s);
  reps[0].Mod(&curve.order);
  reps[1].Set(&glv_lambda);
  reps[1].ModMulK1order(&s);
  reps[2].Set(&glv_lambda2);
  reps[2].ModMulK1order(&s);
  Int best;
  bool have = false;
  r_out = 0;
  for (int r = 0; r < 3; ++r) {
    Int v(reps[r]);
    Int opposite(&curve.order);
    opposite.Sub(&v);  // n - v, lies in (0, n] for v in [0, n)
    bool take_v = v.IsLower(&opposite) || opposite.IsZero();
    Int &candidate = take_v ? v : opposite;
    if (!have || candidate.IsLower(&best)) {
      have = true;
      best.Set(&candidate);
      r_out = r;
    }
  }
  return best;
}

uint64_t ceil_sqrt_u64(uint64_t n) {
  if (n <= 1) return 1;
  uint64_t r = static_cast<uint64_t>(std::sqrt(static_cast<long double>(n)));
  if (r == 0) r = 1;
  while (r != UINT64_MAX && r < (n - 1) / r + 1) ++r;
  while (r > 1 && (r - 1) >= (n - 1) / (r - 1) + 1) --r;
  return r;
}

bool parse_u64(const char *s, uint64_t &out) {
  if (!s || !*s) return false;
  // Reject leading '+'/'-' and blanks: strtoull silently wraps "-1" into
  // ULLONG_MAX, which used to be accepted as a valid width.
  for (const char *p = s; *p; ++p) {
    if (*p == '+' || *p == '-' || *p == ' ' || *p == '\t') return false;
  }
  char *end = nullptr;
  unsigned long long value = std::strtoull(s, &end, 0);
  if (*end != '\0') return false;
  out = static_cast<uint64_t>(value);
  return true;
}

bool parse_int(const std::string &text, Int &out) {
  if (text.size() > 2 && text[0] == '0' &&
      (text[1] == 'x' || text[1] == 'X')) {
    out.SetBase16(text.c_str() + 2);
    return true;
  }
  bool decimal = !text.empty();
  for (char ch : text) decimal = decimal && ch >= '0' && ch <= '9';
  if (!decimal) return false;
  out.SetBase10(text.c_str());
  return true;
}

bool parse_signed_int(const std::string &text, Int &out) {
  if (text.empty()) return false;
  if (text[0] != '-') return parse_int(text, out);
  if (text.size() == 1) return false;
  Int magnitude;
  if (!parse_int(text.substr(1), magnitude)) return false;
  out.Set(&curve.order);
  out.Sub(&magnitude);
  out.Mod(&curve.order);
  return true;
}

// factor * value (mod group order) for full-size scalars.
void set_scaled_mod(Int &out, Int &factor, Int &value) {
  out.Set(&factor);
  out.ModMulK1order(&value);
}

void subtract_order_mod(Int &out, Int &value, Int &shift) {
  out.Set(&value);
  if (out.IsLower(&shift)) out.Add(&curve.order);
  out.Sub(&shift);
  out.Mod(&curve.order);
}

bool candidate_in_cover(const Int &candidate) {
  if (!enforce_cover_range) return true;
  Int value(candidate);
  return !value.IsLower(&cover_start) && value.IsLower(&cover_end);
}

bool check_orbit_candidate(Int &base_scalar, Int &row_scalar,
                           Int &segment_end, Int &candidate_out,
                           uint64_t &offset_out) {
  Int factors[6];
  factors[0].Set(&row_scalar);
  factors[0].Mod(&curve.order);
  set_scaled_mod(factors[1], glv_lambda, row_scalar);
  set_scaled_mod(factors[2], glv_lambda2, row_scalar);
  for (int i = 0; i < 3; ++i) {
    factors[i + 3].Set(&curve.order);
    factors[i + 3].Sub(&factors[i]);
  }
  for (Int &factor : factors) {
    orbit_factor_checks.fetch_add(1, std::memory_order_relaxed);
    Int candidate;
    candidate.Set(&base_scalar);
    candidate.Add(&factor);
    candidate.Mod(&curve.order);
    bool in_interval = !candidate.IsLower(&range_start) &&
                       candidate.IsLower(&segment_end);
    if (!in_interval && !partner_table_mode) continue;
    exact_candidate_checks.fetch_add(1, std::memory_order_relaxed);
    Point check = hit_ec_mult(candidate);
    if (check.equals(search_target)) {
      Int shifted_candidate(&candidate);
      Int recovered;
      subtract_order_mod(recovered, candidate, active_shift);
      Point original_check = hit_ec_mult(recovered);
      if (!original_check.equals(target) ||
          (!partner_table_mode && !candidate_in_cover(recovered))) {
        continue;
      }
      if (!in_interval) {
        partner_recoveries.fetch_add(1, std::memory_order_relaxed);
      }
      candidate_out.Set(&recovered);
      Int delta(&shifted_candidate);
      delta.Sub(&range_start);
      // Partner candidates sit far outside the 64-bit window offset range.
      offset_out = in_interval ? delta.GetInt64() : 0;
      return true;
    }
  }
  return false;
}

// Route a K3 bucket row by its (d01,d12) signature: identify the rotation r
// with giant = phi^r(+-row point) and set base_factor = lambda^r * row_scalar.
// Returns r, or -1 when the row signature is degenerate; the caller falls back
// to the full six-branch check, so coverage never depends on the router.
int route_signature_pair(size_t row_index,
                         std::array<unsigned char, 64> &giant_sig,
                         Int &base_factor) {
  const std::array<unsigned char, 64> &row_sig =
      shared_signatures[row_index];
  alignas(8) unsigned char r01_bytes[32];
  alignas(8) unsigned char r12_bytes[32];
  std::memcpy(r01_bytes, row_sig.data(), 32);
  std::memcpy(r12_bytes, row_sig.data() + 32, 32);
  Int r01, r12, rsum_neg;
  r01.Set32Bytes(r01_bytes);
  r12.Set32Bytes(r12_bytes);
  // rsum_neg = -(r01 + r12) mod P, kept canonical in [0, P): Int::ModNeg
  // maps 0 to P, which would break the byte comparison in the degenerate case.
  rsum_neg.Set(&r01);
  rsum_neg.Add(&r12);
  rsum_neg.Mod(&curve.P);
  if (!rsum_neg.IsZero()) {
    rsum_neg.Neg();
    rsum_neg.Add(&curve.P);
  }
  alignas(8) unsigned char g01_bytes[32];
  alignas(8) unsigned char g12_bytes[32];
  std::memcpy(g01_bytes, giant_sig.data(), 32);
  std::memcpy(g12_bytes, giant_sig.data() + 32, 32);
  Int g01, g12;
  g01.Set32Bytes(g01_bytes);
  g12.Set32Bytes(g12_bytes);

  int route = -1;
  if (g01.IsEqual(&r01) && g12.IsEqual(&r12)) {
    route = 0;
  } else if (g01.IsEqual(&r12) && g12.IsEqual(&rsum_neg)) {
    route = 1;
  } else if (g01.IsEqual(&rsum_neg) && g12.IsEqual(&r01)) {
    route = 2;
  }
  if (route < 0) return -1;

  Int row_scalar;
  row_scalar.Set32Bytes(shared_scalars[row_index].data());
  if (route == 0) base_factor.Set(&row_scalar);
  else if (route == 1) set_scaled_mod(base_factor, glv_lambda, row_scalar);
  else set_scaled_mod(base_factor, glv_lambda2, row_scalar);
  return route;
}

// Exact-check the routed candidate pair {+base_factor, -base_factor}.
bool check_factor_pair(Int &base_factor, Int &base_scalar,
                       Int &segment_end, Int &candidate_out,
                       uint64_t &offset_out) {
  Int neg_factor;
  neg_factor.Set(&curve.order);
  neg_factor.Sub(&base_factor);
  Int factors[2] = {base_factor, neg_factor};
  for (Int &factor : factors) {
    orbit_factor_checks.fetch_add(1, std::memory_order_relaxed);
    Int candidate;
    candidate.Set(&base_scalar);
    candidate.Add(&factor);
    candidate.Mod(&curve.order);
    bool in_interval = !candidate.IsLower(&range_start) &&
                       candidate.IsLower(&segment_end);
    if (!in_interval && !partner_table_mode) continue;
    exact_candidate_checks.fetch_add(1, std::memory_order_relaxed);
    Point check = hit_ec_mult(candidate);
    if (check.equals(search_target)) {
      Int shifted_candidate(&candidate);
      Int recovered;
      subtract_order_mod(recovered, candidate, active_shift);
      Point original_check = hit_ec_mult(recovered);
      if (!original_check.equals(target) ||
          (!partner_table_mode && !candidate_in_cover(recovered))) {
        continue;
      }
      if (!in_interval) {
        partner_recoveries.fetch_add(1, std::memory_order_relaxed);
      }
      candidate_out.Set(&recovered);
      Int delta(&shifted_candidate);
      delta.Sub(&range_start);
      // Partner candidates sit far outside the 64-bit window offset range.
      offset_out = in_interval ? delta.GetInt64() : 0;
      return true;
    }
  }
  return false;
}

// Sort the table together with its parallel arrays by key, apply the
// canonical one-row-per-orbit merge, and refresh the unique-row counters.
void finalize_table(size_t raw_count) {
  const bool sig_mode = glv_orbit_mode && phase_signature;
  std::vector<size_t> order(raw_count);
  for (size_t i = 0; i < raw_count; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [](size_t a, size_t b) {
    return shared_table[a].key < shared_table[b].key;
  });
  std::vector<BabyRow> sorted_rows(raw_count);
  std::vector<std::array<unsigned char, 64>> sorted_sigs;
  std::vector<std::array<unsigned char, 32>> sorted_scalars;
  if (sig_mode) sorted_sigs.resize(raw_count);
  if (glv_orbit_mode) sorted_scalars.resize(raw_count);
  for (size_t w = 0; w < raw_count; ++w) {
    size_t idx = order[w];
    sorted_rows[w] = shared_table[idx];
    if (sig_mode) sorted_sigs[w] = shared_signatures[idx];
    if (glv_orbit_mode) sorted_scalars[w] = shared_scalars[idx];
  }
  shared_table.swap(sorted_rows);
  if (sig_mode) shared_signatures.swap(sorted_sigs);
  if (glv_orbit_mode) shared_scalars.swap(sorted_scalars);

  size_t final_count = raw_count;
  if (glv_orbit_mode && canonical_k3_mode) {
    // Equal K3 keys are exactly one orbit (roots of X^3-c are
    // x, beta*x, beta^2*x), and the canonical scalar is constant on the
    // orbit, so adjacent equal-key rows are duplicates by construction.
    size_t w = 0;
    for (size_t i = 0; i < raw_count; ++i) {
      if (w > 0 && shared_table[i].key == shared_table[w - 1].key) {
        // Equal K3 key: same orbit when the stored scalars agree (they are
        // orbit minima here), a genuine cross-orbit collision otherwise.
        if (std::memcmp(shared_scalars[i].data(),
                        shared_scalars[w - 1].data(), 32) != 0) {
          insert_collisions.fetch_add(1, std::memory_order_relaxed);
          std::fprintf(stderr,
                       "[W] insert collision: equal K3 key with a different "
                       "stored scalar (row dropped)\n");
        }
        merged_rows.fetch_add(1, std::memory_order_relaxed);
        continue;
      }
      shared_table[w] = shared_table[i];
      if (sig_mode) shared_signatures[w] = shared_signatures[i];
      shared_scalars[w] = shared_scalars[i];
      ++w;
    }
    final_count = w;
    shared_table.resize(final_count);
    if (sig_mode) shared_signatures.resize(final_count);
    shared_scalars.resize(final_count);
  }

  uint64_t unique = final_count == 0 ? 0 : 1;
  for (size_t i = 1; i < final_count; ++i) {
    if (!(shared_table[i].key == shared_table[i - 1].key)) ++unique;
  }
  raw_baby_rows.store(static_cast<uint64_t>(raw_count),
                      std::memory_order_relaxed);
  baby_rows.store(unique, std::memory_order_relaxed);
  // Generic build-phase merges (duplicate keys): counted only outside the
  // canonical path, where finalize_table already merged those rows above
  // (raw_count - unique == compaction drops there; adding both would
  // double-count one collapse).
  if (!(glv_orbit_mode && canonical_k3_mode)) {
    merged_rows.fetch_add(raw_count > unique ? raw_count - unique : 0,
                          std::memory_order_relaxed);
  }
  shared_baby_rows = final_count;
}

// ---------------- CBT1: canonical baby table, binary format v1 -------------
// payload_kind 1 (giant_grid): rows of {K3 key: 33 B} + {canonical scalar:
// 32 B}, 65 B/row, sorted ascending by key. The file is the serialized
// search state of the in-RAM canonical table: the giant scan binary-searches
// the mmap'd keys directly, so a cold start needs no build EC work and no
// in-RAM table (RSS stays at the scan working set even at H=4M, where the
// text loader would hold ~576 MB of Int vectors).
// Header: 256 bytes, big-endian, fixed offsets:
//    0 magic "CBT1" |  4 canonical_version |  8 payload_kind
//   12 header_size  | 16 row_size         | 20 curve_id (1 = secp256k1)
//   24 orbit_size(6)| 28 flags (bit0: keys sorted)
//   32 row_count    | 40 stride           | 48 m | 56 H | 64 window
//   72 rows_per_head_source | 76 reserved
//   80 n (32 B) | 112 p (32 B) | 144 lambda (32 B)
//  176 created_unix | 184..239 reserved
//  240 payload_fnv1a64 | 248 header_fnv1a64 (bytes 0..247)
// Reject policy: any mismatch (version, kind, sizes, curve params, checksums,
// sortedness, row_count vs file size) exits with code 2 - never "read how
// much there is". canonical_version future-proofs the format: a file written
// by an older generator cannot silently load as a newer layout.
static constexpr char CBT1_MAGIC[4] = {'C', 'B', 'T', '1'};
static constexpr uint32_t CBT1_VERSION = 1;
static constexpr uint32_t CBT1_PAYLOAD_GIANT_GRID = 1;
static constexpr uint32_t CBT1_CURVE_SECP256K1 = 1;
static constexpr uint32_t CBT1_HEADER_SIZE = 256;
static constexpr uint32_t CBT1_ROW_SIZE = 65;
static std::string table_binary_path;       // --table-binary (mmap load)
static std::string dump_table_binary_path;  // --dump-table-binary (write)
static const unsigned char *file_rows = nullptr;  // payload base (row 0)
static uint64_t file_row_count = 0;

struct MappedTable {
  const unsigned char *base = nullptr;
  size_t size = 0;
  ~MappedTable() {
    if (base) ::munmap(const_cast<unsigned char *>(base), size);
  }
};
static MappedTable table_map;

uint64_t fnv1a_bytes(const unsigned char *data, size_t n) {
  uint64_t hash = 0xcbf29ce484222325ULL;
  for (size_t i = 0; i < n; ++i) {
    hash ^= data[i];
    hash *= 0x100000001b3ULL;
  }
  return hash;
}

void put_u32be(unsigned char *d, uint32_t v) {
  d[0] = (unsigned char)(v >> 24); d[1] = (unsigned char)(v >> 16);
  d[2] = (unsigned char)(v >> 8);  d[3] = (unsigned char)v;
}

void put_u64be(unsigned char *d, uint64_t v) {
  for (int i = 0; i < 8; ++i) d[i] = (unsigned char)(v >> (56 - 8 * i));
}

uint32_t get_u32be(const unsigned char *d) {
  return (uint32_t(d[0]) << 24) | (uint32_t(d[1]) << 16) |
         (uint32_t(d[2]) << 8) | uint32_t(d[3]);
}

uint64_t get_u64be(const unsigned char *d) {
  uint64_t v = 0;
  for (int i = 0; i < 8; ++i) v = (v << 8) | d[i];
  return v;
}

int load_table_binary(const std::string &path) {
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0) {
    std::fprintf(stderr, "[E] table binary: cannot open %s\n", path.c_str());
    return 2;
  }
  struct stat st {};
  if (::fstat(fd, &st) != 0 ||
      st.st_size < static_cast<off_t>(CBT1_HEADER_SIZE)) {
    std::fprintf(stderr, "[E] table binary: %s too small for a CBT1 header\n",
                 path.c_str());
    ::close(fd);
    return 2;
  }
  const size_t size = static_cast<size_t>(st.st_size);
  void *base = ::mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
  ::close(fd);
  if (base == MAP_FAILED) {
    std::fprintf(stderr, "[E] table binary: mmap failed for %s\n",
                 path.c_str());
    return 2;
  }
  const unsigned char *h = static_cast<const unsigned char *>(base);
  auto reject = [&](const char *reason) {
    std::fprintf(stderr, "[E] table binary %s: %s (reject)\n",
                 path.c_str(), reason);
    ::munmap(base, size);
    return 2;
  };
  if (std::memcmp(h, CBT1_MAGIC, 4) != 0) return reject("bad magic");
  const uint32_t version = get_u32be(h + 4);
  if (version != CBT1_VERSION) {
    std::fprintf(stderr,
                 "[E] table binary %s: canonical_version %u != %u - the file "
                 "was written by a different format generation\n",
                 path.c_str(), version, CBT1_VERSION);
    ::munmap(base, size);
    return 2;
  }
  if (get_u32be(h + 8) != CBT1_PAYLOAD_GIANT_GRID)
    return reject("unsupported payload_kind");
  if (get_u32be(h + 12) != CBT1_HEADER_SIZE) return reject("header_size");
  if (get_u32be(h + 16) != CBT1_ROW_SIZE) return reject("row_size");
  if (get_u32be(h + 20) != CBT1_CURVE_SECP256K1) return reject("curve_id");
  const uint64_t rows = get_u64be(h + 32);
  const uint64_t stride = get_u64be(h + 40);
  if (rows == 0) return reject("row_count is zero");
  if (CBT1_HEADER_SIZE + rows * CBT1_ROW_SIZE != size ||
      (size - CBT1_HEADER_SIZE) / CBT1_ROW_SIZE != rows)
    return reject("row_count does not match file size");
  alignas(8) unsigned char expect[32];
  curve.order.Get32Bytes(expect);
  if (std::memcmp(h + 80, expect, 32) != 0)
    return reject("group order mismatch");
  curve.P.Get32Bytes(expect);
  if (std::memcmp(h + 112, expect, 32) != 0)
    return reject("field prime mismatch");
  glv_lambda.Get32Bytes(expect);
  if (std::memcmp(h + 144, expect, 32) != 0)
    return reject("GLV lambda mismatch");
  if (fnv1a_bytes(h, 248) != get_u64be(h + 248))
    return reject("header checksum mismatch");
  const unsigned char *payload = h + CBT1_HEADER_SIZE;
  if (fnv1a_bytes(payload, static_cast<size_t>(rows) * CBT1_ROW_SIZE) !=
      get_u64be(h + 240))
    return reject("payload checksum mismatch");
  for (uint64_t i = 1; i < rows; ++i) {
    if (std::memcmp(payload + i * CBT1_ROW_SIZE,
                    payload + (i - 1) * CBT1_ROW_SIZE, 33) < 0) {
      std::fprintf(stderr,
                   "[E] table binary %s: keys not sorted ascending at row "
                   "%llu (reject)\n",
                   path.c_str(), static_cast<unsigned long long>(i));
      ::munmap(base, size);
      return 2;
    }
  }
  table_map.base = static_cast<const unsigned char *>(base);
  table_map.size = size;
  file_rows = payload;
  file_row_count = rows;
  // Validation touched every payload page; drop them from the working set
  // so the scan re-faults pages on demand (read-only mapping faults back
  // from the file). Keeps RSS at the scan working set, not the file size.
  ::madvise(const_cast<unsigned char *>(table_map.base), table_map.size,
            MADV_DONTNEED);
  std::cout << "  table binary file: " << path << " (" << rows
            << " rows, payload giant_grid, 65 B/row, checksum + sortedness"
            << " validated)\n"
            << "  table binary H: " << get_u64be(h + 56)
            << "  rows/head(source): " << get_u32be(h + 72) << "\n";
  if (stride != rows) {
    std::cout << "  table binary stride field: " << stride
              << " (row_count " << rows << " wins as the giant stride)\n";
  }
  return 0;
}

int dump_table_binary(const std::string &path) {
  const size_t rows = shared_table.size();
  if (rows == 0 || rows != shared_scalars.size()) {
    std::fprintf(stderr,
                 "[E] dump table binary: table/scalars size mismatch "
                 "(%zu vs %zu)\n",
                 rows, shared_scalars.size());
    return 2;
  }
  std::vector<unsigned char> buf(CBT1_HEADER_SIZE +
                                 rows * CBT1_ROW_SIZE);
  unsigned char *h = buf.data();
  std::memcpy(h, CBT1_MAGIC, 4);
  put_u32be(h + 4, CBT1_VERSION);
  put_u32be(h + 8, CBT1_PAYLOAD_GIANT_GRID);
  put_u32be(h + 12, CBT1_HEADER_SIZE);
  put_u32be(h + 16, CBT1_ROW_SIZE);
  put_u32be(h + 20, CBT1_CURVE_SECP256K1);
  put_u32be(h + 24, 6);            // orbit_size: full C6 GLV orbit
  put_u32be(h + 28, 1);            // flags: keys sorted ascending
  put_u64be(h + 32, rows);
  put_u64be(h + 40, shared_baby_rows);  // giant stride == stored rows
  put_u64be(h + 48, rows);              // m
  put_u64be(h + 56, 0);                 // H unknown for a generic table
  put_u64be(h + 64, 0);                 // window unknown
  put_u32be(h + 72, 1);                 // canonical source rows/head
  curve.order.Get32Bytes(h + 80);
  curve.P.Get32Bytes(h + 112);
  glv_lambda.Get32Bytes(h + 144);
  put_u64be(h + 176, static_cast<uint64_t>(std::time(nullptr)));
  unsigned char *payload = h + CBT1_HEADER_SIZE;
  for (size_t i = 0; i < rows; ++i) {
    unsigned char *row = payload + i * CBT1_ROW_SIZE;
    std::memcpy(row, shared_table[i].key.data(), 33);
    std::memcpy(row + 33, shared_scalars[i].data(), 32);
    if (i > 0 && std::memcmp(row, row - CBT1_ROW_SIZE, 33) < 0) {
      std::fprintf(stderr,
                   "[E] dump table binary: table not sorted ascending at "
                   "row %zu\n",
                   i);
      return 2;
    }
  }
  put_u64be(h + 240, fnv1a_bytes(payload, rows * CBT1_ROW_SIZE));
  put_u64be(h + 248, fnv1a_bytes(h, 248));
  // Atomic write (tmp + fsync + rename), same discipline as checkpoints.
  const std::string tmp = path + ".tmp";
  int fd = ::open(tmp.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (fd < 0) {
    std::fprintf(stderr, "[E] dump table binary: cannot open %s\n",
                 tmp.c_str());
    return 2;
  }
  size_t off = 0;
  while (off < buf.size()) {
    ssize_t n = ::write(fd, buf.data() + off, buf.size() - off);
    if (n <= 0) {
      std::fprintf(stderr, "[E] dump table binary: write failed\n");
      ::close(fd);
      return 2;
    }
    off += static_cast<size_t>(n);
  }
  ::fsync(fd);
  ::close(fd);
  if (::rename(tmp.c_str(), path.c_str()) != 0) {
    std::fprintf(stderr, "[E] dump table binary: rename failed\n");
    return 2;
  }
  std::cout << "  dumped canonical table: " << path << " (" << rows
            << " rows, " << buf.size() << " bytes, payload giant_grid)\n";
  return 0;
}

// Per-row GLV payload: K3 key, 32-byte row scalar and (optionally) the phase
// signature. In canonical mode the stored scalar is the orbit minimum and the
// key/signature are taken at the phi-rotated x-coordinate beta^r * x — the
// x-coordinate of the canonical point (negation does not change x), so no
// extra scalar multiplication is needed.
struct GlvRowData {
  std::array<unsigned char, 33> key;
  std::array<unsigned char, 32> scalar;
  std::array<unsigned char, 64> sig;  // meaningful when phase_signature
};

GlvRowData build_glv_row(Int &row_scalar, Int &row_x, bool want_sig) {
  GlvRowData out;
  Int stored(&row_scalar);
  Int sig_x(&row_x);
  if (canonical_k3_mode) {
    int rotation = 0;
    Int canon = canonical_orbit_scalar(row_scalar, rotation);
    stored.Set(&canon);
    for (int i = 0; i < rotation; ++i) {
      sig_x.ModMulK1(&beta, &sig_x);
    }
  }
  out.key = invariant_key_from_x(sig_x);
  alignas(8) unsigned char scalar_bytes[32];
  stored.Get32Bytes(scalar_bytes);
  std::memcpy(out.scalar.data(), scalar_bytes, 32);
  if (want_sig) out.sig = phase_signature_from_x(sig_x);
  return out;
}

void build_shared_table(uint64_t m) {
  shared_table.resize(static_cast<size_t>(m));
  const bool sig_mode = glv_orbit_mode && phase_signature;
  if (glv_orbit_mode) shared_scalars.resize(static_cast<size_t>(m));
  if (sig_mode) shared_signatures.resize(static_cast<size_t>(m));
  std::vector<std::thread> builders;
  builders.reserve(worker_count);
  for (unsigned worker = 0; worker < worker_count; ++worker) {
    uint64_t begin = m * worker / worker_count;
    uint64_t end = m * (worker + 1) / worker_count;
    if (begin == end) continue;
    builders.emplace_back([begin, end, sig_mode]() {
      Point current;
      if (begin == 0) {
        current.Clear();
      } else {
        Int scalar;
        scalar.SetInt64(begin);
        current = ec_mult(scalar);
      }
      for (uint64_t j = begin; j < end; ++j) {
        if (glv_orbit_mode) {
          Int row_scalar;
          row_scalar.SetInt64(j);
          GlvRowData row = build_glv_row(row_scalar, current.x, sig_mode);
          shared_table[static_cast<size_t>(j)] = {row.key, j};
          shared_scalars[static_cast<size_t>(j)] = row.scalar;
          if (sig_mode) {
            shared_signatures[static_cast<size_t>(j)] = row.sig;
          }
        } else {
          shared_table[static_cast<size_t>(j)] = {point_key(current), j};
        }
        current = add_safe(current, curve.G);
      }
    });
  }
  for (auto &builder : builders) builder.join();
  finalize_table(static_cast<size_t>(m));
}

// --table-scalars: build the table from an explicit scalar list (one row per
// scalar). The point of each scalar is computed directly; the giant plan then
// strides by the stored row count, so targets of the form
// start + i*m + s (s in the orbit set of the table) remain recoverable.
// In --canonical-k3 mode main() has already canonicalized and deduped the
// list (direct canonical pass before any EC work), so the builder receives
// exactly one scalar per orbit.
void build_table_from_scalars(std::vector<Int> &scalars) {
  size_t m = scalars.size();
  const bool sig_mode = glv_orbit_mode && phase_signature;
  shared_table.resize(m);
  if (glv_orbit_mode) shared_scalars.resize(m);
  if (sig_mode) shared_signatures.resize(m);
  std::vector<std::thread> builders;
  builders.reserve(worker_count);
  for (unsigned worker = 0; worker < worker_count; ++worker) {
    size_t begin = m * worker / worker_count;
    size_t end = m * (worker + 1) / worker_count;
    if (begin == end) continue;
    builders.emplace_back([begin, end, sig_mode, &scalars]() {
      for (size_t i = begin; i < end; ++i) {
        Point row_point = ec_mult(scalars[i]);
        GlvRowData row = build_glv_row(scalars[i], row_point.x, sig_mode);
        shared_table[i] = {row.key, 0};
        shared_scalars[i] = row.scalar;
        if (sig_mode) shared_signatures[i] = row.sig;
      }
    });
  }
  for (auto &builder : builders) builder.join();
  finalize_table(m);
}

void search_slice(uint64_t local_start, uint64_t local_size,
                  uint64_t segment) {
  if (local_size == 0 || local_start >= active_range_size) return;
  local_size = std::min(local_size, active_range_size - local_start);
  uint64_t m = shared_baby_rows;
  uint64_t giant_count = (local_size - 1) / m + 1;
  giant_steps_planned.fetch_add(giant_count, std::memory_order_relaxed);

  Int anchor_scalar(&range_start);
  anchor_scalar.Add(local_start);
  Int segment_end(&anchor_scalar);
  segment_end.Add(local_size);
  Point anchor = ec_mult(anchor_scalar);
  Point neg_anchor = curve.Negation(anchor);
  Point giant_zero = add_safe(search_target, neg_anchor);

  Int m_int(m);
  Point step = ec_mult(m_int);
  Point neg_step = curve.Negation(step);

  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (unsigned worker = 0; worker < worker_count; ++worker) {
    uint64_t begin = giant_count * worker / worker_count;
    uint64_t end = giant_count * (worker + 1) / worker_count;
    if (begin == end) continue;
    workers.emplace_back([&, begin, end, segment, local_start]() {
      Point giant = giant_zero;
      Int base_scalar(&anchor_scalar);
      if (begin != 0) {
        Int delta;
        delta.SetInt64(begin);
        delta.Mult(m);
        delta.Mod(&curve.order);
        base_scalar.Add(&delta);
        base_scalar.Mod(&curve.order);
        Point shift = ec_mult(delta);
        Point neg_shift = curve.Negation(shift);
        giant = add_safe(giant, neg_shift);
      }
      std::array<unsigned char, 64> giant_sig{};
      bool giant_sig_ready = false;

      // Table backend: in-RAM rows or CBT1 mmap rows (65 B stride = 33-byte
      // key + 32-byte canonical scalar). Both are sorted ascending by key;
      // memcmp order == unsigned lexicographic order of
      // std::array<unsigned char, 33>, so one index-based binary search
      // serves both backends.
      const bool file_mode = file_rows != nullptr;
      const size_t row_count = file_mode
          ? static_cast<size_t>(file_row_count)
          : shared_table.size();
      auto row_key = [&](size_t idx) -> const unsigned char * {
        return file_mode ? file_rows + idx * 65
                         : shared_table[idx].key.data();
      };

      for (uint64_t i = begin; i < end &&
                              !found.load(std::memory_order_relaxed); ++i) {
        giant_steps_scheduled.fetch_add(1, std::memory_order_relaxed);
        giant_steps.fetch_add(1, std::memory_order_relaxed);
        auto query_key = glv_orbit_mode ? invariant_key(giant) : point_key(giant);
        auto key_matches = [&](size_t idx) -> bool {
          return idx < row_count &&
                 std::memcmp(row_key(idx), query_key.data(), 33) == 0;
        };
        size_t lo = 0, hi = row_count;
        while (lo < hi) {
          size_t mid = lo + (hi - lo) / 2;
          if (std::memcmp(row_key(mid), query_key.data(), 33) < 0) lo = mid + 1;
          else hi = mid;
        }
        size_t row_index = lo;
        if (key_matches(row_index)) {
          // One search hit per matched giant position (bucket), independent
          // of how many rows the bucket holds - the cross-mode comparable
          // hit count.
          search_hits.fetch_add(1, std::memory_order_relaxed);
        }
        bool bucket_checked = false;
        // Signature-path bucket state: the first routed pair and whether the
        // full six-branch set has been checked (anomaly fallback).
        bool six_branch_done = false;
        bool have_pair = false;
        Int pair_b;
        Int pair_neg_b;
        while (key_matches(row_index) &&
               !found.load(std::memory_order_relaxed)) {
          if (glv_orbit_mode) {
            invariant_hits.fetch_add(1, std::memory_order_relaxed);
            if (!phase_signature && bucket_checked) {
              // Every row of a K3=x^3 bucket lies in one C6 orbit (the roots
              // of X^3-c over F_p are exactly x, beta*x, beta^2*x), so the
              // six GLV branches are the same scalar set for every row of
              // the bucket. The first row's check already covered them.
              // In non-GLV mode rows give distinct candidates (base + j),
              // so no dedup is applied there.
              orbit_duplicate_rows.fetch_add(1, std::memory_order_relaxed);
              ++row_index;
              continue;
            }
          }
          Int candidate;
          uint64_t local_offset = 0;
          bool match = false;
          if (glv_orbit_mode) {
            if (phase_signature && six_branch_done) {
              // The six-branch set of this orbit was already checked; any
              // routed pair is a subset of it.
              orbit_duplicate_rows.fetch_add(1, std::memory_order_relaxed);
              ++row_index;
              continue;
            }
            Int row_scalar;
            if (file_mode) {
              row_scalar.Set32Bytes(
                  const_cast<unsigned char *>(file_rows + row_index * 65 +
                                              33));
            } else {
              row_scalar.Set32Bytes(shared_scalars[row_index].data());
            }
            if (phase_signature) {
              if (!giant_sig_ready) {
                giant_sig = phase_signature_of(giant);
                giant_sig_ready = true;
              }
              Int base_factor;
              int route = route_signature_pair(row_index, giant_sig,
                                               base_factor);
              if (route < 0) {
                signature_anomalies.fetch_add(1, std::memory_order_relaxed);
                match = check_orbit_candidate(base_scalar, row_scalar,
                                              segment_end, candidate,
                                              local_offset);
                six_branch_done = true;
              } else {
                signature_routes.fetch_add(1, std::memory_order_relaxed);
                Int neg_factor;
                neg_factor.Set(&curve.order);
                neg_factor.Sub(&base_factor);
                if (have_pair) {
                  bool same_pair =
                      (base_factor.IsEqual(&pair_b) &&
                       neg_factor.IsEqual(&pair_neg_b)) ||
                      (base_factor.IsEqual(&pair_neg_b) &&
                       neg_factor.IsEqual(&pair_b));
                  if (same_pair) {
                    // All rows of a K3 bucket route to the same candidate
                    // pair {+o, -o} (validated in
                    // phase_signature_economy_experiment.py: 100% collapse),
                    // so the first row's exact check already covered them.
                    signature_pair_collapses.fetch_add(
                        1, std::memory_order_relaxed);
                    ++row_index;
                    continue;
                  }
                  // Structural-invariant violation: verify defensively.
                  signature_pair_mismatches.fetch_add(
                      1, std::memory_order_relaxed);
                } else {
                  have_pair = true;
                  pair_b.Set(&base_factor);
                  pair_neg_b.Set(&neg_factor);
                }
                match = check_factor_pair(base_factor, base_scalar,
                                          segment_end, candidate,
                                          local_offset);
              }
            } else {
              match = check_orbit_candidate(base_scalar, row_scalar,
                                            segment_end, candidate,
                                            local_offset);
            }
          } else {
            uint64_t j = shared_table[row_index].j;
            local_offset = local_start + i * m + j;
            if (local_offset < local_start + local_size) {
              candidate.Set(&base_scalar);
              candidate.Add(j);
              exact_candidate_checks.fetch_add(1, std::memory_order_relaxed);
              Point check = hit_ec_mult(candidate);
              match = check.equals(search_target);
              if (match) {
                Int recovered;
                subtract_order_mod(recovered, candidate, active_shift);
                Point original_check = hit_ec_mult(recovered);
                match = original_check.equals(target) && candidate_in_cover(recovered);
                if (match) {
                  candidate.Set(&recovered);
                }
              }
            }
          }
          if (match) {
            std::lock_guard<std::mutex> lock(found_mutex);
            if (!found.load()) {
              found_scalar.Set(&candidate);
              Int shifted_candidate(&candidate);
              shifted_candidate.Add(&active_shift);
              shifted_candidate.Mod(&curve.order);
              found_shifted_scalar.Set(&shifted_candidate);
              found_shift.Set(&active_shift);
              found_segment = segment;
              found_local_offset = local_offset;
              found_run_us.store(run_us(), std::memory_order_relaxed);
              found.store(true);
            }
            break;
          }
          if (glv_orbit_mode) {
            bucket_checked = true;
          }
          ++row_index;
        }
        giant = add_safe(giant, neg_step);
        base_scalar.Add(&m_int);
        base_scalar.Mod(&curve.order);
      }
    });
  }
  for (auto &worker : workers) worker.join();
  segments_done.fetch_add(1, std::memory_order_relaxed);
}

void search_segment(uint64_t segment) {
  uint64_t local_start = segment * segment_size;
  if (local_start >= active_range_size) return;
  search_slice(local_start, std::min(segment_size, active_range_size - local_start),
               segment);
}

std::vector<uint64_t> lambda_priority_order(Point &point) {
  std::vector<uint64_t> result;
  if (!lambda_priority || lambda_sector_count == 0) return result;
  lambda_calculations.fetch_add(3, std::memory_order_relaxed);
  std::array<Int, 3> values = lambda_triangle_values(point);
  for (Int &value : values) {
    value.Mult(lambda_sector_count);
    value.Div(&curve.P);
    uint64_t sector = value.GetInt64();
    if (sector >= lambda_sector_count) sector = lambda_sector_count - 1;
    if (std::find(result.begin(), result.end(), sector) == result.end()) {
      result.push_back(sector);
    }
  }
  return result;
}

void search_current_target() {
  uint64_t sector_count = lambda_priority
      ? std::min(lambda_sector_count, active_range_size == 0 ? uint64_t(1) : active_range_size)
      : 0;
  if (!lambda_priority) {
    uint64_t segment_count = active_range_size / segment_size +
                             (active_range_size % segment_size != 0);
    for (uint64_t segment = 0; segment < segment_count &&
                                 !found.load(std::memory_order_relaxed); ++segment) {
      search_segment(segment);
    }
    return;
  }

  std::vector<uint64_t> priority = lambda_priority_order(search_target);
  std::vector<bool> used(static_cast<size_t>(sector_count), false);
  auto search_sector = [&](uint64_t sector, bool is_priority) {
    if (sector >= sector_count || used[static_cast<size_t>(sector)]) return;
    used[static_cast<size_t>(sector)] = true;
    uint64_t begin = (active_range_size * sector) / sector_count;
    uint64_t end = (active_range_size * (sector + 1)) / sector_count;
    if (end <= begin) return;
    if (is_priority) priority_slices.fetch_add(1, std::memory_order_relaxed);
    else fallback_slices.fetch_add(1, std::memory_order_relaxed);
    for (uint64_t pos = begin; pos < end &&
                               !found.load(std::memory_order_relaxed); ) {
      uint64_t amount = std::min(segment_size, end - pos);
      search_slice(pos, amount, pos / segment_size);
      pos += amount;
    }
  };
  for (uint64_t sector : priority) {
    if (found.load(std::memory_order_relaxed)) break;
    search_sector(sector, true);
  }
  for (uint64_t sector = 0; sector < sector_count &&
                             !found.load(std::memory_order_relaxed); ++sector) {
    search_sector(sector, false);
  }
}

// ---------------------------------------------------------------------------
// --cover checkpoint: atomic progress file for long automatic-coverage runs.
//
// Granularity is the external cover chunk (one run_target call over range_size
// at --cover granularity). The file is written to <path>.tmp, fsynced and
// renamed - a killed process leaves the last completed chunk recorded and the
// file never appears half-written. On --resume the plan fingerprint is
// re-checked (target/start/size/cover/modes/table file) so a checkpoint of a
// different plan cannot silently mislead the continuation, and the
// search-phase counters are restored, so printed totals stay cumulative
// across kills/sessions. Build-phase counters (canonical merge, table rows)
// are not restored: the table itself is rebuilt by design.
// ---------------------------------------------------------------------------

constexpr unsigned kCheckpointVersion = 1;

std::string int_hex_64(Int &v) {
  // zero-padded 64-digit hex, human-comparable across sessions
  char *raw = v.GetBase16();
  std::string digits(raw);
  std::free(raw);
  return std::string(digits.size() < 64 ? 64 - digits.size() : 0, '0') + digits;
}

uint64_t fnv1a(const std::string &s) {
  uint64_t hash = 0xcbf29ce484222325ULL;
  for (unsigned char ch : s) {
    hash ^= ch;
    hash *= 0x100000001b3ULL;
  }
  return hash;
}

std::string cover_fingerprint(const std::string &target_hex,
                              const std::string &start_text,
                              uint64_t size, const std::string &cover_start_text,
                              const std::string &cover_size_text,
                              uint64_t segment, const std::string &table_file,
                              bool glv, bool sig, bool partner, bool canon,
                              bool priority) {
  std::string s;
  s.reserve(256);
  s += target_hex;
  s += "|";
  s += start_text;
  s += "|" + std::to_string(size);
  s += "|" + cover_start_text;
  s += "|" + cover_size_text;
  s += "|" + std::to_string(segment);
  s += "|" + table_file;
  s += "|" + std::to_string(glv);
  s += "|" + std::to_string(sig);
  s += "|" + std::to_string(partner);
  s += "|" + std::to_string(canon);
  s += "|" + std::to_string(priority);
  char buf[32];
  snprintf(buf, sizeof(buf), "%016llx",
           static_cast<unsigned long long>(fnv1a(s)));
  return buf;
}

// (name, getter) pairs describe every restored search-phase counter.
struct CheckpointCounter {
  const char *name;
  const std::atomic<uint64_t> &counter;
};

std::vector<CheckpointCounter> checkpoint_counter_list() {
  return {
      {"giant steps", giant_steps},
      {"giant steps scheduled", giant_steps_scheduled},
      {"giant steps planned", giant_steps_planned},
      {"invariant hits", invariant_hits},
      {"search hits", search_hits},
      {"orbit factor checks", orbit_factor_checks},
      {"exact candidate checks", exact_candidate_checks},
      {"EC scalar multiplications", ec_scalar_mults},
      {"EC point additions", ec_point_adds},
      {"hit EC scalar multiplications", hit_ec_scalar_mults},
      {"signature routes", signature_routes},
      {"signature anomalies", signature_anomalies},
      {"signature pair collapses", signature_pair_collapses},
      {"signature pair mismatches", signature_pair_mismatches},
      {"orbit duplicate rows", orbit_duplicate_rows},
      {"partner recoveries", partner_recoveries},
      {"lambda calculations", lambda_calculations},
      {"lambda priority slices", priority_slices},
      {"lambda fallback slices", fallback_slices},
      {"segments done", segments_done},
  };
}

void atomic_write_text(const std::string &path, const std::string &content) {
  std::string tmp = path + ".tmp";
  int fd = ::open(tmp.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (fd < 0) {
    std::fprintf(stderr, "[W] checkpoint: cannot open %s\n", tmp.c_str());
    return;
  }
  size_t off = 0;
  while (off < content.size()) {
    ssize_t n = ::write(fd, content.data() + off, content.size() - off);
    if (n <= 0) {
      ::close(fd);
      std::fprintf(stderr, "[W] checkpoint: write failed for %s\n",
                   tmp.c_str());
      return;
    }
    off += static_cast<size_t>(n);
  }
  if (::fsync(fd) != 0) {
    std::fprintf(stderr, "[W] checkpoint: fsync failed for %s\n", tmp.c_str());
  }
  ::close(fd);
  if (::rename(tmp.c_str(), path.c_str()) != 0) {
    std::fprintf(stderr, "[W] checkpoint: rename failed for %s\n", tmp.c_str());
  }
}

void save_cover_checkpoint(const std::string &fingerprint,
                           const Int &source_start, uint64_t chunks_done,
                           double cumulative_elapsed, bool found_state) {
  std::string text;
  text.reserve(2048);
  text += "cover checkpoint version: "
      + std::to_string(kCheckpointVersion) + "\n";
  text += "plan fingerprint: " + fingerprint + "\n";
  text += "next source chunk start: 0x" + int_hex_64(
      const_cast<Int &>(source_start)) + "\n";
  text += "chunks completed: " + std::to_string(chunks_done) + "\n";
  char buf[64];
  snprintf(buf, sizeof(buf), "cumulative elapsed seconds: %.3f\n",
           cumulative_elapsed);
  text += buf;
  text += "found: " + std::to_string(found_state ? 1 : 0) + "\n";
  if (found_state) {
    text += "recovered scalar: 0x" + int_hex_64(found_scalar) + "\n";
    text += "recovered private key: 0x" + int_hex_64(found_scalar) + "\n";
  }
  for (CheckpointCounter &cc : checkpoint_counter_list()) {
    text += std::string("counter ") + cc.name + ": "
        + std::to_string(cc.counter.load(std::memory_order_relaxed)) + "\n";
  }
  atomic_write_text(checkpoint_path, text);
}

bool parse_checkpoint_file(const std::string &path,
                           std::map<std::string, std::string> &fields) {
  // The main file is only created by rename after a full fsync, so it is
  // always whole; a leftover <path>.tmp means the CURRENT save did not
  // finish and is ignored (the previous completed save stays on disk).
  std::ifstream file(path);
  if (!file) return false;
  std::string line;
  while (std::getline(file, line)) {
    size_t colon = line.find(':');
    if (colon == std::string::npos) continue;
    std::string key = line.substr(0, colon);
    std::string value = line.substr(colon + 1);
    while (!value.empty() && (value.front() == ' ' || value.front() == '\t')) {
      value.erase(value.begin());
    }
    fields[key] = value;
  }
  return true;
}

// Reads and validates the checkpoint file. Counter restoration is kept in a
// separate step (restore_checkpoint_counters) because the engine takes the
// fresh build-phase EC snapshot after the table build and BEFORE resuming, so
// the printed build line stays pure and the search totals stay cumulative.
bool read_cover_checkpoint(const std::string &fingerprint,
                           const std::map<std::string, std::string> &_unused_fields,
                           Int &source_start, uint64_t &chunks_done,
                           double &cumulative_elapsed, bool &found_state,
                           std::map<std::string, std::string> &fields,
                           std::string &error) {
  (void)_unused_fields;
  if (!parse_checkpoint_file(checkpoint_path, fields)) {
    error = "no readable checkpoint file: " + checkpoint_path;
    return false;
  }
  auto field = [&fields](const char *k, std::string &out) -> bool {
    auto it = fields.find(k);
    if (it == fields.end()) return false;
    out = it->second;
    return true;
  };
  std::string version, saved_fingerprint, source_hex, chunks_text, elapsed_text,
      found_text;
  if (!field("cover checkpoint version", version) ||
      !field("plan fingerprint", saved_fingerprint) ||
      !field("next source chunk start", source_hex) ||
      !field("chunks completed", chunks_text) ||
      !field("cumulative elapsed seconds", elapsed_text) ||
      !field("found", found_text)) {
    error = "checkpoint file is missing required fields";
    return false;
  }
  if (std::strtoul(version.c_str(), nullptr, 10) != kCheckpointVersion) {
    error = "checkpoint version mismatch: " + version;
    return false;
  }
  if (saved_fingerprint != fingerprint) {
    error = "checkpoint belongs to a different plan (fingerprint " +
            saved_fingerprint + ", expected " + fingerprint + ")";
    return false;
  }
  if (!parse_int(source_hex, source_start)) {
    error = "checkpoint source_start is not a valid number: " + source_hex;
    return false;
  }
  chunks_done = std::strtoull(chunks_text.c_str(), nullptr, 10);
  cumulative_elapsed = std::strtod(elapsed_text.c_str(), nullptr);
  found_state = std::strtoul(found_text.c_str(), nullptr, 10) != 0;
  if (found_state) {
    auto it = fields.find("recovered scalar");
    if (it != fields.end() && !parse_int(it->second, found_scalar)) {
      error = "checkpoint recovered scalar is not a valid number: "
              + it->second;
      return false;
    }
  }
  return true;
}

// Restores the search-phase counters into the live atomics (fetch_add).
void restore_checkpoint_counters(
    const std::map<std::string, std::string> &fields) {
  for (CheckpointCounter &cc : checkpoint_counter_list()) {
    auto it = fields.find(std::string("counter ") + cc.name);
    if (it == fields.end()) continue;  // counters added since the save
    uint64_t value = std::strtoull(it->second.c_str(), nullptr, 10);
    if (std::string(cc.name) == "giant steps") {
      restored_giant_steps.store(value, std::memory_order_relaxed);
    }
    const_cast<std::atomic<uint64_t> &>(cc.counter).fetch_add(
        value, std::memory_order_relaxed);
  }
}

// Progress ticker for long runs: one stdout line every progress_interval
// seconds of search time - per-session giant steps/s (cumulative minus the
// --resume restored offset) and current RSS, for external stability logs.
static std::atomic<uint64_t> progress_interval{60};
static std::atomic<bool> progress_stop{false};

static long current_rss_kb() {
  FILE *f = std::fopen("/proc/self/status", "r");
  if (!f) return -1;
  char line[256];
  long rss = -1;
  while (std::fgets(line, sizeof(line), f)) {
    if (std::sscanf(line, "VmRSS: %ld kB", &rss) == 1) break;
  }
  std::fclose(f);
  return rss;
}

static void progress_ticker(
    std::chrono::steady_clock::time_point search_t0) {
  const uint64_t interval = progress_interval.load();
  if (interval == 0) return;
  while (!progress_stop.load(std::memory_order_relaxed)) {
    for (unsigned i = 0;
         i < static_cast<unsigned>(interval) * 10 &&
         !progress_stop.load(std::memory_order_relaxed);
         ++i) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    if (progress_stop.load(std::memory_order_relaxed)) break;
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - search_t0).count();
    const uint64_t session_steps =
        giant_steps.load() - restored_giant_steps.load();
    std::printf("progress: elapsed %.1fs | giant steps %llu | "
                "steps/s %.1f | rss kB %ld\n",
                elapsed, static_cast<unsigned long long>(session_steps),
                elapsed > 0 ? session_steps / elapsed : 0.0,
                current_rss_kb());
    std::fflush(stdout);
  }
}

void usage(const char *name) {
  std::fprintf(stderr,
      "Usage: %s --target PUBKEY --start HEX --size N [options]\n"
      "  -q, --target PUBKEY   SEC1 compressed/uncompressed public key\n"
      "  -a, --start VALUE     interval start (hex 0x... or decimal)\n"
      "  -w, --size N          interval width (finite bounded test)\n"
      "  -s, --segment N       local segment width (default 100000)\n"
      "      --baby-rows N     choose segment width as N^2\n"
      "      --shift T         search Q+T*G before the original target (repeatable)\n"
      "      --glv-orbit       use K3=x^3 quotient plus exact GLV branches\n"
      "      --phase-signature route GLV branches by the (d01,d12) Lambda signature (requires --glv-orbit)\n"
      "      --canonical-k3   store one canonical row per C6 orbit: scalar min(k, n-k) over lambda^r*k (requires --glv-orbit)\n"
      "      --table-scalars FILE  build baby rows from explicit hex scalars, one per line (requires --glv-orbit)\n"
      "      --table-binary FILE  mmap a CBT1 canonical table as the table backend (payload giant_grid)\n"
      "      --dump-table-binary FILE  write the built table as a CBT1 file (65 B/row, sorted, checksummed)\n"
      "      --partner-table  exact-verify out-of-interval orbit branches (requires --glv-orbit)\n"
      "      --lambda-priority order local slices by the Lambda triangle, then fallback\n"
      "      --lambda-sectors N number of local numeric slices (default 16)\n"
      "      --cover-start V   external coverage interval start (Int/hex)\n"
      "      --cover-size V    external coverage interval width (Int/hex)\n"
      "      --checkpoint-file PATH  atomic progress checkpoint for --cover runs (one per cover chunk)\n"
      "      --checkpoint-steps N    write the checkpoint every N chunks (default 1)\n"
      "      --resume          resume the --cover plan from --checkpoint-file (fingerprint is verified)\n"
      "      --progress-interval N  print a progress line (steps/s, RSS) every N seconds, 0 = off (default 60)\n"
      "      --lambda-dump     print the exact three Lambda values and exit\n"
      "  -t, --threads N       worker threads (default 1)\n"
      "  -h, --help            show this help\n", name);
}

}  // namespace

int main(int argc, char **argv) {
  const option options[] = {
      {"target", required_argument, nullptr, 'q'},
      {"start", required_argument, nullptr, 'a'},
      {"size", required_argument, nullptr, 'w'},
      {"segment", required_argument, nullptr, 's'},
      {"baby-rows", required_argument, nullptr, 'b'},
      {"shift", required_argument, nullptr, 'd'},
      {"glv-orbit", no_argument, nullptr, 'g'},
      {"lambda-priority", no_argument, nullptr, 'l'},
      {"lambda-sectors", required_argument, nullptr, 'L'},
      {"phase-signature", no_argument, nullptr, 'P'},
      {"partner-table", no_argument, nullptr, 'p'},
      {"canonical-k3", no_argument, nullptr, 'k'},
      {"table-scalars", required_argument, nullptr, 'S'},
      {"table-binary", required_argument, nullptr, 'B'},
      {"dump-table-binary", required_argument, nullptr, 'T'},
      {"cover-start", required_argument, nullptr, 'c'},
      {"cover-size", required_argument, nullptr, 'C'},
      {"checkpoint-file", required_argument, nullptr, 'x'},
      {"checkpoint-steps", required_argument, nullptr, 'y'},
      {"resume", no_argument, nullptr, 'z'},
      {"progress-interval", required_argument, nullptr, 'R'},
      {"lambda-dump", no_argument, nullptr, 'D'},
      {"threads", required_argument, nullptr, 't'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}};
  std::string target_hex;
  std::string start_text;
  std::vector<std::string> shift_texts;
  std::string cover_start_text;
  std::string cover_size_text;
  bool lambda_dump = false;
  int c;
  while ((c = getopt_long(argc, argv, "q:a:w:s:t:b:d:glL:PpkS:B:T:c:C:x:y:zR:Dh", options, nullptr)) != -1) {
    switch (c) {
      case 'q': target_hex = optarg; break;
      case 'a': start_text = optarg; break;
      case 'w': if (!parse_u64(optarg, range_size)) { usage(argv[0]); return 2; } break;
      case 's': if (!parse_u64(optarg, segment_size)) { usage(argv[0]); return 2; } break;
      case 'b': if (!parse_u64(optarg, requested_baby_rows) || requested_baby_rows == 0) { usage(argv[0]); return 2; } break;
      case 'd': shift_texts.emplace_back(optarg); break;
      case 'g': glv_orbit_mode = true; break;
      case 'l': lambda_priority = true; break;
      case 'L': if (!parse_u64(optarg, lambda_sector_count) || lambda_sector_count == 0) { usage(argv[0]); return 2; } break;
      case 'P': phase_signature = true; break;
      case 'p': partner_table_mode = true; break;
      case 'k': canonical_k3_mode = true; break;
      case 'S': table_scalars_path = optarg; break;
      case 'B': table_binary_path = optarg; break;
      case 'T': dump_table_binary_path = optarg; break;
      case 'c': cover_start_text = optarg; break;
      case 'C': cover_size_text = optarg; break;
      case 'x': checkpoint_path = optarg; break;
      case 'y': if (!parse_u64(optarg, checkpoint_steps) || checkpoint_steps == 0) { usage(argv[0]); return 2; } break;
      case 'z': checkpoint_resume = true; break;
      case 'R': {
        uint64_t n;
        if (!parse_u64(optarg, n)) {
          usage(argv[0]);
          return 2;
        }
        progress_interval.store(n, std::memory_order_relaxed);
        break;
      }
      case 'D': lambda_dump = true; break;
      case 't': { uint64_t n; if (!parse_u64(optarg, n) || n == 0 || n > 4096) { usage(argv[0]); return 2; } worker_count = static_cast<unsigned>(n); break; }
      case 'h': usage(argv[0]); return 0;
      default: usage(argv[0]); return 2;
    }
  }
  if (target_hex.empty() || start_text.empty() || range_size == 0 || segment_size == 0) {
    usage(argv[0]); return 2;
  }
  if (phase_signature && !glv_orbit_mode) {
    std::fprintf(stderr, "[E] --phase-signature requires --glv-orbit\n");
    return 2;
  }
  if (partner_table_mode && !glv_orbit_mode) {
    std::fprintf(stderr, "[E] --partner-table requires --glv-orbit\n");
    return 2;
  }
  if (canonical_k3_mode && !glv_orbit_mode) {
    std::fprintf(stderr, "[E] --canonical-k3 requires --glv-orbit\n");
    return 2;
  }
  if (!table_scalars_path.empty() && !glv_orbit_mode) {
    std::fprintf(stderr, "[E] --table-scalars requires --glv-orbit\n");
    return 2;
  }
  if (!table_scalars_path.empty() && requested_baby_rows != 0) {
    std::fprintf(stderr, "[E] --baby-rows cannot be combined with --table-scalars\n");
    return 2;
  }
  if (requested_baby_rows != 0) {
    if (requested_baby_rows > std::numeric_limits<uint64_t>::max() / requested_baby_rows) {
      std::fprintf(stderr, "[E] baby row count is too large\n");
      return 2;
    }
    segment_size = requested_baby_rows * requested_baby_rows;
  }

  bool automatic_coverage = !cover_start_text.empty() || !cover_size_text.empty();
  if (automatic_coverage && (cover_start_text.empty() || cover_size_text.empty())) {
    std::fprintf(stderr, "[E] --cover-start and --cover-size must be supplied together\n");
    return 2;
  }
  if (automatic_coverage && !shift_texts.empty()) {
    std::fprintf(stderr, "[E] manual --shift cannot be combined with automatic coverage\n");
    return 2;
  }
  if (checkpoint_resume && checkpoint_path.empty()) {
    std::fprintf(stderr, "[E] --resume requires --checkpoint-file\n");
    return 2;
  }
  if (!checkpoint_path.empty() && !automatic_coverage) {
    std::fprintf(stderr, "[E] --checkpoint-file/--resume require --cover-start/--cover-size\n");
    return 2;
  }

  curve.Init();
  beta.SetBase16("7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE");
  beta_g.ModMulK1(&beta, &curve.G.x);
  beta2.ModMulK1(&beta, &beta);
  glv_lambda.SetBase16("5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72");
  glv_lambda2.SetBase16("AC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE");
  // --table-binary: mmap a CBT1 file as the table backend. Mutually
  // exclusive with --table-scalars; requires --glv-orbit (the payload keys
  // are K3 invariant keys). Phase-signature routing is unavailable in file
  // mode (the payload carries keys and scalars, not row signatures), but the
  // hit path degrades to the direct orbit check with checks/hit still 1.0.
  if (!table_binary_path.empty()) {
    if (!table_scalars_path.empty()) {
      std::fprintf(stderr,
                   "[E] --table-binary cannot be combined with --table-scalars\n");
      return 2;
    }
    if (!dump_table_binary_path.empty()) {
      std::fprintf(stderr,
                   "[E] --table-binary cannot be combined with "
                   "--dump-table-binary\n");
      return 2;
    }
    if (!glv_orbit_mode) {
      std::fprintf(stderr,
                   "[E] --table-binary requires --glv-orbit (payload keys are "
                   "K3 invariant keys)\n");
      return 2;
    }
    int rc = load_table_binary(table_binary_path);
    if (rc != 0) return rc;
    if (phase_signature) {
      std::fprintf(stderr,
                   "[W] --table-binary has no row signatures: disabling "
                   "phase-signature routing (hit path falls back to the "
                   "direct orbit check)\n");
      phase_signature = false;
    }
  }
  if (!parse_int(start_text, range_start)) {
    std::fprintf(stderr, "[E] start must be hexadecimal (0x...) or decimal\n");
    return 2;
  }
  Int range_end(&range_start);
  range_end.Add(range_size);
  if (!range_start.IsStrictPositive() || range_end.IsGreater(&curve.order)) {
    std::fprintf(stderr, "[E] invalid interval: start + size must be <= group order\n");
    return 2;
  }
  Int cover_size;
  double checkpoint_elapsed_before = 0.0;  // from --resume, for final print
  if (automatic_coverage) {
    if (!parse_int(cover_start_text, cover_start) ||
        !parse_int(cover_size_text, cover_size) ||
        !cover_start.IsStrictPositive() || !cover_size.IsStrictPositive()) {
      std::fprintf(stderr, "[E] cover interval must contain positive Int values\n");
      return 2;
    }
    cover_end.Set(&cover_start);
    cover_end.Add(&cover_size);
    if (cover_end.IsGreater(&curve.order)) {
      std::fprintf(stderr, "[E] cover-start + cover-size must be <= group order\n");
      return 2;
    }
  }
  std::vector<char> mutable_key(target_hex.begin(), target_hex.end());
  mutable_key.push_back('\0');
  bool compressed = false;
  if (!curve.ParsePublicKeyHex(mutable_key.data(), target, compressed)) {
    std::fprintf(stderr, "[E] invalid SEC1 public key\n");
    return 2;
  }
  if (lambda_dump) {
    std::array<Int, 3> values = lambda_triangle_values(target);
    for (size_t i = 0; i < values.size(); ++i) {
      char *value = values[i].GetBase16();
      std::cout << "lambda" << i << ": 0x" << value << "\n";
      std::free(value);
    }
    return 0;
  }

  std::vector<Int> requested_shifts;
  requested_shifts.reserve(shift_texts.size());
  for (const std::string &text : shift_texts) {
    Int shift;
    if (!parse_signed_int(text, shift)) {
      std::fprintf(stderr, "[E] shift must be a signed hexadecimal or decimal scalar\n");
      return 2;
    }
    requested_shifts.push_back(shift);
  }

  std::vector<Int> table_scalars;
  if (!table_scalars_path.empty()) {
    std::ifstream in(table_scalars_path);
    if (!in) {
      std::fprintf(stderr, "[E] cannot open table scalars file: %s\n",
                   table_scalars_path.c_str());
      return 2;
    }
    std::string line;
    uint64_t line_number = 0;
    while (std::getline(in, line)) {
      ++line_number;
      size_t comment = line.find('#');
      if (comment != std::string::npos) line.erase(comment);
      size_t first = line.find_first_not_of(" \t\r\n");
      if (first == std::string::npos) continue;
      size_t last = line.find_last_not_of(" \t\r\n");
      line = line.substr(first, last - first + 1);
      Int value;
      if (!parse_int(line, value)) {
        std::fprintf(stderr, "[E] invalid scalar on line %llu of %s\n",
                     static_cast<unsigned long long>(line_number),
                     table_scalars_path.c_str());
        return 2;
      }
      value.Mod(&curve.order);
      if (value.IsZero()) {
        std::fprintf(stderr, "[E] zero scalar on line %llu of %s\n",
                     static_cast<unsigned long long>(line_number),
                     table_scalars_path.c_str());
        return 2;
      }
      table_scalars.push_back(value);
    }
    if (table_scalars.empty()) {
      std::fprintf(stderr, "[E] table scalars file is empty: %s\n",
                   table_scalars_path.c_str());
      return 2;
    }
  }
  rss_after_load_kb = current_rss_kb();
  hwm_after_load_kb = peak_rss_kb();

  // --canonical-k3 + --table-scalars: direct canonical pass. Every input
  // scalar is canonicalized first (two order multiplications per scalar,
  // no EC work) and equal canonical scalars - exactly one orbit - are
  // deduped BEFORE any point is computed, so build cost and peak memory
  // follow the stored table, not the raw input list. One point per orbit
  // instead of one point per input row; merged_rows counts input
  // rows dropped before the first scalar multiplication.
  const uint64_t raw_table_scalars =
      table_scalars.empty() ? 0 : static_cast<uint64_t>(table_scalars.size());
  if (!table_scalars.empty() && glv_orbit_mode && canonical_k3_mode) {
    const size_t raw_input = table_scalars.size();
    // In-place canonical pass: canonicalize every input scalar (two order
    // multiplications per scalar, no EC work), then sort+unique in the same
    // buffer. The previous two-vector version (table_scalars + canon) held
    // 2 x raw_input Ints and left the merged table living in a
    // raw_input-capacity buffer after the swap - the whole +61 MB steady
    // RSS tax of merge mode at H=400k. In-place + shrink_to_fit keeps the
    // peak at one input buffer and the steady state at the stored table.
    for (Int &value : table_scalars) {
      int rotation_unused;
      value = canonical_orbit_scalar(value, rotation_unused);
    }
    std::sort(table_scalars.begin(), table_scalars.end(),
              [](const Int &a, const Int &b) {
                return const_cast<Int &>(a).IsLower(const_cast<Int *>(&b));
              });
    table_scalars.erase(
        std::unique(table_scalars.begin(), table_scalars.end(),
                    [](const Int &a, const Int &b) {
                      return const_cast<Int &>(a).IsEqual(const_cast<Int *>(&b));
                    }),
        table_scalars.end());
    // Pre-pass drops: input rows folded into their orbit representative
    // before any point is computed. Deterministic canonicalization makes
    // these merges by construction - they are never search events.
    merged_rows.fetch_add(raw_input - table_scalars.size(),
                          std::memory_order_relaxed);
    table_scalars.shrink_to_fit();
  }
  rss_after_prepass_kb = current_rss_kb();
  hwm_after_prepass_kb = peak_rss_kb();

  uint64_t segment_count = range_size / segment_size + (range_size % segment_size != 0);
  uint64_t m = table_scalars.empty()
      ? (file_rows != nullptr
             ? file_row_count
             : ceil_sqrt_u64(std::min(segment_size, range_size)))
      : static_cast<uint64_t>(table_scalars.size());
  shared_baby_rows = m;
  std::cout << "Experimental bounded segmented BSGS\n"
            << "  interval width: " << range_size << "\n"
            << "  segments: " << segment_count << "\n"
            << "  segment width: " << segment_size << "\n"
            << "  workers: " << worker_count << "\n"
            << "  GLV orbit quotient: " << (glv_orbit_mode ? "true" : "false") << "\n"
            << "  phase signature routing: " << (phase_signature ? "true" : "false") << "\n"
            << "  canonical K3 storage: " << (canonical_k3_mode ? "true" : "false") << "\n"
            << "  partner table: " << (partner_table_mode ? "true" : "false") << "\n"
            << "  shifted targets: " << requested_shifts.size()
            << " (original target is always the fallback)\n"
            << "  baby rows (shared): " << m << "\n";
  if (!table_scalars.empty()) {
    std::cout << "  table scalars file: " << table_scalars_path << " ("
              << table_scalars.size() << " rows)\n";
  }
  if (file_rows != nullptr) {
    std::cout << "  table backend: mmap (CBT1 payload, no in-RAM table, "
              << "no build EC work)\n";
  }
  std::cout << "  estimated shared table bytes: " << (m * sizeof(BabyRow)) << "\n"
            << "  table copies: 1 (read-only, shared by workers)\n";
  if (lambda_priority) {
    std::cout << "  Lambda priority: true\n"
              << "  Lambda sectors: " << lambda_sector_count << "\n";
  }
  if (automatic_coverage) {
    char *cover_start_hex = cover_start.GetBase16();
    char *cover_size_hex = cover_size.GetBase16();
    std::cout << "  external cover start: 0x" << cover_start_hex << "\n"
              << "  external cover size: 0x" << cover_size_hex << "\n";
    std::free(cover_start_hex);
    std::free(cover_size_hex);
  }
  if (requested_baby_rows != 0) {
    std::cout << "  requested baby rows: " << requested_baby_rows << "\n";
    if (m < requested_baby_rows) {
      std::cout << "  effective baby rows limited by interval: " << m << "\n"
                << "  requested segment width: "
                << (requested_baby_rows * requested_baby_rows) << "\n";
    }
  }

  std::cout << "  building baby table with " << worker_count
            << " builder threads..." << std::endl;
  auto started = std::chrono::steady_clock::now();
  auto table_started = started;
  program_started = started;
  if (!table_scalars.empty()) {
    build_table_from_scalars(table_scalars);
  } else if (file_rows != nullptr) {
    // CBT1 file mode: the mmap payload IS the table. No build, no in-RAM
    // rows; the giant scan binary-searches the file directly. The stride is
    // the stored row count by construction (stride == merged rows).
    m = file_row_count;
    shared_baby_rows = m;
  } else {
    build_shared_table(m);
  }
  double table_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - table_started).count();
  rss_after_build_kb = current_rss_kb();
  hwm_after_build_kb = peak_rss_kb();
  // Return freed build transients (finalize_table sorted copies, the
  // shrink_to_fit victim buffer) to the OS. glibc raises the mmap
  // threshold after a large mmap'd block is freed, so those copies land in
  // the arena and stay cached as resident pages (~30 MB at H=400k) - the
  // printed RSS timeline should reflect live data, not allocator caching.
  malloc_trim(0);
  rss_after_build_kb = current_rss_kb();
  if (canonical_k3_mode) {
    std::cout << "  table rows stored: " << shared_baby_rows
              << " (build merged rows " << merged_rows.load() << ")\n";
  }
  if (!table_scalars_path.empty()) {
    // Hard invariant behind "giant stride == post-merge rows": the stored
    // table plus the rows merged away while building (pre-pass canonical
    // drops + finalize_table compaction) must reconstruct the raw scalar
    // file exactly. A mismatch would mean the giant plan strides a
    // different grid than the table covers.
    const uint64_t stored = shared_baby_rows;
    const uint64_t merged = merged_rows.load();
    // Canonical mode removes rows from the table (pre-pass dedup +
    // finalize_table compaction), so raw == stored + merged. Non-canonical
    // mode keeps every row (duplicate keys share a bucket), so nothing is
    // removed and raw == stored; merged_rows there counts duplicate keys.
    const uint64_t expected_raw =
        canonical_k3_mode ? stored + merged : stored;
    if (expected_raw != raw_table_scalars) {
      std::fprintf(stderr,
                   "[E] canonical stride accounting broken: stored %llu + "
                   "merged %llu != raw %llu input rows\n",
                   static_cast<unsigned long long>(stored),
                   static_cast<unsigned long long>(merged),
                   static_cast<unsigned long long>(raw_table_scalars));
      std::exit(2);
    }
    std::cout << "  giant stride (post-merge rows): " << stored << "\n"
              << std::flush;
  } else if (file_rows != nullptr) {
    // CBT1 file mode: stride == stored row_count by construction (the
    // payload IS the post-merge table), asserted by the loader's
    // row_count-vs-file-size check.
    std::cout << "  giant stride (post-merge rows): " << shared_baby_rows
              << "\n"
              << std::flush;
  }
  if (!dump_table_binary_path.empty()) {
    int rc = dump_table_binary(dump_table_binary_path);
    if (rc != 0) return rc;
  }

  auto search_started = std::chrono::steady_clock::now();
  // Fresh-build EC snapshot: taken BEFORE the checkpoint counters are
  // restored, so "(build)" lines below stay this-session-pure and any
  // prior-session work appears only in the cumulative totals.
  const uint64_t build_mults_snapshot = ec_scalar_mults.load();
  const uint64_t build_adds_snapshot = ec_point_adds.load();
  std::thread ticker(progress_ticker, search_started);
  auto run_target = [&](Int &shift, uint64_t local_size, bool cover_check) {
    active_shift.Set(&shift);
    active_range_size = local_size;
    enforce_cover_range = cover_check;
    Int normalized_shift(&active_shift);
    normalized_shift.Mod(&curve.order);
    Point shift_point;
    if (normalized_shift.IsZero()) shift_point.Clear();
    else shift_point = ec_mult(normalized_shift);
    search_target = add_safe(target, shift_point);
    search_current_target();
  };

  if (automatic_coverage) {
    const std::string plan_fingerprint = cover_fingerprint(
        target_hex, start_text, range_size, cover_start_text,
        cover_size_text, segment_size, table_scalars_path, glv_orbit_mode,
        phase_signature, partner_table_mode, canonical_k3_mode,
        lambda_priority);
    Int source_start(&cover_start);
    uint64_t chunks_before = 0;
    double elapsed_before = 0.0;
    bool resumed_found = false;    if (checkpoint_resume) {
      std::map<std::string, std::string> saved_fields;
      std::string resume_error;
      if (!read_cover_checkpoint(plan_fingerprint, saved_fields, source_start,
                                 chunks_before, elapsed_before, resumed_found,
                                 saved_fields, resume_error)) {
        std::fprintf(stderr, "[E] --resume: %s\n", resume_error.c_str());
        progress_stop.store(true, std::memory_order_relaxed);
        ticker.join();
        return 2;
      }
      restore_checkpoint_counters(saved_fields);
      checkpoint_elapsed_before = elapsed_before;
      if (resumed_found) {
        std::cout << "  checkpoint resumed: the target was already "
                  << "recovered in an earlier session, nothing to search\n"
                  << "  recovered private key: 0x"
                  << int_hex_64(found_scalar) << "\n"
                  << "  recovery correct: true (from checkpoint)\n";
        progress_stop.store(true, std::memory_order_relaxed);
        ticker.join();
        return 0;
      }
      char *resumed_hex = source_start.GetBase16();
      std::cout << "  checkpoint resumed: chunks completed before this "
                   "session: " << chunks_before
                << ", next source chunk start: 0x" << resumed_hex << "\n"
                << "  checkpoint resume: search counters restored "
                << "(printed totals are cumulative across sessions)\n"
                << std::endl;
      std::free(resumed_hex);
    }
    const auto cover_session_started = std::chrono::steady_clock::now();
    uint64_t automatic_chunks = 0;
    while (source_start.IsLower(&cover_end)) {
      Int remaining(&cover_end);
      remaining.Sub(&source_start);
      Int local_capacity(range_size);
      uint64_t chunk_size;
      if (remaining.IsGreaterOrEqual(&local_capacity)) chunk_size = range_size;
      else chunk_size = remaining.GetInt64();
      if (chunk_size == 0) break;
      Int shift;
      subtract_order_mod(shift, range_start, source_start);
      if (automatic_chunks < 8 || automatic_chunks % 1000 == 0) {
        char *source_hex = source_start.GetBase16();
        std::cout << "  automatic source chunk start: 0x" << source_hex
                  << " size: " << chunk_size << "\n" << std::flush;
        std::free(source_hex);
      }
      run_target(shift, chunk_size, true);
      source_start.Add(chunk_size);
      ++automatic_chunks;
      const bool final_chunk = !source_start.IsLower(&cover_end);
      if (!checkpoint_path.empty() &&
          (final_chunk || automatic_chunks % checkpoint_steps == 0)) {
        double cumulative = elapsed_before
            + std::chrono::duration<double>(std::chrono::steady_clock::now()
                - cover_session_started).count();
        save_cover_checkpoint(plan_fingerprint, source_start,
                              chunks_before + automatic_chunks, cumulative,
                              found.load(std::memory_order_relaxed));
        if (automatic_chunks < 8 || automatic_chunks % 1000 == 0) {
          char *next_hex = source_start.GetBase16();
          std::cout << "  checkpoint saved: chunks "
                    << (chunks_before + automatic_chunks)
                    << ", next source chunk start: 0x" << next_hex << "\n"
                    << std::flush;
          std::free(next_hex);
        }
      }
      if (found.load(std::memory_order_relaxed)) {
        break;
      }
    }
    std::cout << "  automatic chunks attempted: " << automatic_chunks << "\n";
    if (!found.load(std::memory_order_relaxed)) {
      Int zero_shift;
      zero_shift.SetInt32(0);
      run_target(zero_shift, range_size, true);
      if (!checkpoint_path.empty() && found.load(std::memory_order_relaxed)) {
        // the fallback pass (zero shift) found the target: record it
        double cumulative = elapsed_before
            + std::chrono::duration<double>(std::chrono::steady_clock::now()
                - cover_session_started).count();
        save_cover_checkpoint(plan_fingerprint, source_start,
                              chunks_before + automatic_chunks, cumulative,
                              true);
      }
    }
  } else {
    for (Int requested : requested_shifts) {
      if (found.load(std::memory_order_relaxed)) break;
      Int shift(&requested);
      run_target(shift, range_size, false);
    }
    if (!found.load(std::memory_order_relaxed)) {
      Int zero_shift;
      zero_shift.SetInt32(0);
      run_target(zero_shift, range_size, false);
    }
  }
  progress_stop.store(true, std::memory_order_relaxed);
  ticker.join();
  double search_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - search_started).count();
  double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();

  const uint64_t total_mults = ec_scalar_mults.load();
  const uint64_t total_adds = ec_point_adds.load();
  const uint64_t hit_mults = hit_ec_scalar_mults.load();
  const uint64_t search_mults = total_mults - build_mults_snapshot;
  std::cout << "  segments completed: " << segments_done.load() << "\n"
            << "  raw baby rows built: " << raw_baby_rows.load() << "\n"
            << "  unique baby rows built: " << baby_rows.load() << "\n"
            << "  duplicate baby keys: "
            << (raw_baby_rows.load() - baby_rows.load()) << "\n"
            << "  giant steps: " << giant_steps.load() << "\n"
            << "  giant steps scheduled: " << giant_steps_scheduled.load() << "\n"
            << "  invariant hits: " << invariant_hits.load() << "\n"
            << "  search hits: " << search_hits.load() << "\n"
            << "  merged rows: " << merged_rows.load() << "\n"
            << "  insert collisions: " << insert_collisions.load() << "\n"
            << "  orbit factors scheduled: " << orbit_factor_checks.load() << "\n"
            << "  exact candidate checks: " << exact_candidate_checks.load() << "\n"
            << "  EC scalar multiplications: " << total_mults << "\n"
            << "  EC scalar multiplications (build): " << build_mults_snapshot << "\n"
            << "  EC scalar multiplications (giant search, cumulative): "
            << (search_mults - hit_mults) << "\n"
            << "  EC scalar multiplications (exact candidate checks): "
            << hit_mults << "\n"
            << "  EC point additions: " << total_adds << "\n"
            << "  EC point additions (build): " << build_adds_snapshot << "\n"
            << "  EC point additions (giant search, cumulative): "
            << (total_adds - build_adds_snapshot) << "\n"
            << "  total EC operations: "
            << (ec_scalar_mults.load() + ec_point_adds.load()) << "\n"
            << "  Lambda calculations: " << lambda_calculations.load() << "\n"
            << "  Lambda priority slices: " << priority_slices.load() << "\n"
            << "  Lambda fallback slices: " << fallback_slices.load() << "\n";
  if (phase_signature) {
    std::cout << "  signature routes: " << signature_routes.load() << "\n"
              << "  signature anomalies (6-branch fallback): "
              << signature_anomalies.load() << "\n"
              << "  signature pair collapses: "
              << signature_pair_collapses.load() << "\n"
              << "  signature pair mismatches: "
              << signature_pair_mismatches.load() << "\n";
  }
  if (glv_orbit_mode) {
    std::cout << "  orbit duplicate rows skipped: "
              << orbit_duplicate_rows.load() << "\n";
    if (partner_table_mode) {
      std::cout << "  partner recoveries: "
                << partner_recoveries.load() << "\n";
    }
  }
  std::cout << "  giant steps planned: " << giant_steps_planned.load() << "\n"
            << "  table build seconds: " << table_seconds << "\n"
            << "  RSS after table load: " << rss_after_load_kb
            << " kB (peak " << hwm_after_load_kb << " kB)\n";
  if (hwm_after_prepass_kb != 0) {
    std::cout << "  RSS after canonical pre-pass: " << rss_after_prepass_kb
              << " kB (peak " << hwm_after_prepass_kb << " kB)\n";
  }
  std::cout << "  RSS after table build: " << rss_after_build_kb
            << " kB (peak " << hwm_after_build_kb << " kB)\n"
            << "  giant search seconds: " << search_seconds << "\n"
            << "  giant steps/second: "
            << (search_seconds > 0
                    ? (giant_steps.load() - restored_giant_steps.load()) /
                          search_seconds
                    : 0)
            << "\n"
            << "  total seconds: " << seconds << "\n";
  if (automatic_coverage && checkpoint_resume) {
    std::cout << "  cumulative elapsed seconds (checkpoint): "
              << (checkpoint_elapsed_before + seconds) << "\n";
  }
  if (found.load()) {
    std::cout << "  time to solution seconds: "
              << (static_cast<double>(found_run_us.load()) / 1e6)
              << " (includes table build)\n"
              << "  exact candidate checks per hit: "
              << (search_hits.load() > 0
                      ? static_cast<double>(exact_candidate_checks.load()) /
                            search_hits.load()
                      : static_cast<double>(exact_candidate_checks.load()))
              << " (exact " << exact_candidate_checks.load()
              << " / search hits " << search_hits.load() << ")\n";
  }
  if (found.load()) {
    char *hex = found_scalar.GetBase16();
    std::string short_hex(hex);
    std::cout << "  recovered scalar: 0x" << short_hex << "\n"
            << "  recovered private key: 0x"
              << std::string(short_hex.size() < 64 ? 64 - short_hex.size() : 0, '0')
              << short_hex << "\n"
              << "  shifted target scalar: 0x";
    char *shifted_hex = found_shifted_scalar.GetBase16();
    std::cout << std::string(shifted_hex).append("\n")
              << "  applied shift: 0x";
    char *shift_hex = found_shift.GetBase16();
    std::cout << std::string(shift_hex).append("\n")
              << "  recovered segment: " << found_segment << "\n"
              << "  recovered local offset: " << found_local_offset << "\n"
              << "  recovery correct: true\n";
    std::free(shifted_hex);
    std::free(shift_hex);
    std::free(hex);
  } else {
    std::cout << "  recovery correct: false (no scalar in interval)\n";
  }
  return found.load() ? 0 : 1;
}
