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
#include <fstream>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
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
std::atomic<uint64_t> exact_candidate_checks{0};
std::atomic<uint64_t> orbit_factor_checks{0};
std::atomic<uint64_t> signature_routes{0};
std::atomic<uint64_t> signature_anomalies{0};
std::atomic<uint64_t> orbit_duplicate_rows{0};
std::atomic<uint64_t> signature_pair_collapses{0};
std::atomic<uint64_t> signature_pair_mismatches{0};
std::atomic<uint64_t> canonical_rows_merged{0};
std::atomic<uint64_t> partner_recoveries{0};
// Honest EC accounting: every scalar multiplication (ComputePublicKey) and
// every real point addition (add_safe with both operands non-zero) counts
// here, so builds and search paths compare on the same unit.
std::atomic<uint64_t> ec_scalar_mults{0};
std::atomic<uint64_t> ec_point_adds{0};
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
    Point check = ec_mult(candidate);
    if (check.equals(search_target)) {
      Int shifted_candidate(&candidate);
      Int recovered;
      subtract_order_mod(recovered, candidate, active_shift);
      Point original_check = ec_mult(recovered);
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
    Point check = ec_mult(candidate);
    if (check.equals(search_target)) {
      Int shifted_candidate(&candidate);
      Int recovered;
      subtract_order_mod(recovered, candidate, active_shift);
      Point original_check = ec_mult(recovered);
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
        canonical_rows_merged.fetch_add(1, std::memory_order_relaxed);
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
  shared_baby_rows = final_count;
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

      for (uint64_t i = begin; i < end &&
                              !found.load(std::memory_order_relaxed); ++i) {
        giant_steps_scheduled.fetch_add(1, std::memory_order_relaxed);
        giant_steps.fetch_add(1, std::memory_order_relaxed);
        auto query_key = glv_orbit_mode ? invariant_key(giant) : point_key(giant);
        auto it = std::lower_bound(shared_table.begin(), shared_table.end(), query_key,
            [](const BabyRow &row, const std::array<unsigned char, 33> &needle) {
              return row.key < needle;
            });
        bool bucket_checked = false;
        // Signature-path bucket state: the first routed pair and whether the
        // full six-branch set has been checked (anomaly fallback).
        bool six_branch_done = false;
        bool have_pair = false;
        Int pair_b;
        Int pair_neg_b;
        while (it != shared_table.end() && it->key == query_key &&
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
              ++it;
              continue;
            }
          }
          size_t row_index = static_cast<size_t>(it - shared_table.begin());
          Int candidate;
          uint64_t local_offset = 0;
          bool match = false;
          if (glv_orbit_mode) {
            if (phase_signature && six_branch_done) {
              // The six-branch set of this orbit was already checked; any
              // routed pair is a subset of it.
              orbit_duplicate_rows.fetch_add(1, std::memory_order_relaxed);
              ++it;
              continue;
            }
            Int row_scalar;
            row_scalar.Set32Bytes(shared_scalars[row_index].data());
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
                    ++it;
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
            uint64_t j = it->j;
            local_offset = local_start + i * m + j;
            if (local_offset < local_start + local_size) {
              candidate.Set(&base_scalar);
              candidate.Add(j);
              exact_candidate_checks.fetch_add(1, std::memory_order_relaxed);
              Point check = ec_mult(candidate);
              match = check.equals(search_target);
              if (match) {
                Int recovered;
                subtract_order_mod(recovered, candidate, active_shift);
                Point original_check = ec_mult(recovered);
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
              found.store(true);
            }
            break;
          }
          if (glv_orbit_mode) {
            bucket_checked = true;
          }
          ++it;
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
      "      --partner-table  exact-verify out-of-interval orbit branches (requires --glv-orbit)\n"
      "      --lambda-priority order local slices by the Lambda triangle, then fallback\n"
      "      --lambda-sectors N number of local numeric slices (default 16)\n"
      "      --cover-start V   external coverage interval start (Int/hex)\n"
      "      --cover-size V    external coverage interval width (Int/hex)\n"
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
      {"cover-start", required_argument, nullptr, 'c'},
      {"cover-size", required_argument, nullptr, 'C'},
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
  while ((c = getopt_long(argc, argv, "q:a:w:s:t:b:d:glL:PpkS:c:C:Dh", options, nullptr)) != -1) {
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
      case 'c': cover_start_text = optarg; break;
      case 'C': cover_size_text = optarg; break;
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

  curve.Init();
  beta.SetBase16("7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE");
  beta_g.ModMulK1(&beta, &curve.G.x);
  beta2.ModMulK1(&beta, &beta);
  glv_lambda.SetBase16("5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72");
  glv_lambda2.SetBase16("AC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE");
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

  // --canonical-k3 + --table-scalars: direct canonical pass. Every input
  // scalar is canonicalized first (two order multiplications per scalar,
  // no EC work) and equal canonical scalars - exactly one orbit - are
  // deduped BEFORE any point is computed, so build cost and peak memory
  // follow the stored table, not the raw input list. One point per orbit
  // instead of one point per input row; canonical_rows_merged counts input
  // rows dropped before the first scalar multiplication.
  if (!table_scalars.empty() && glv_orbit_mode && canonical_k3_mode) {
    size_t raw_input = table_scalars.size();
    std::vector<Int> canon;
    canon.reserve(raw_input);
    for (Int &value : table_scalars) {
      int rotation_unused;
      canon.push_back(canonical_orbit_scalar(value, rotation_unused));
    }
    std::sort(canon.begin(), canon.end(), [](const Int &a, const Int &b) {
      return const_cast<Int &>(a).IsLower(const_cast<Int *>(&b));
    });
    canon.erase(
        std::unique(canon.begin(), canon.end(),
                    [](const Int &a, const Int &b) {
                      return const_cast<Int &>(a).IsEqual(const_cast<Int *>(&b));
                    }),
        canon.end());
    canonical_rows_merged.fetch_add(raw_input - canon.size(),
                                    std::memory_order_relaxed);
    table_scalars.swap(canon);
  }

  uint64_t segment_count = range_size / segment_size + (range_size % segment_size != 0);
  uint64_t m = table_scalars.empty()
      ? ceil_sqrt_u64(std::min(segment_size, range_size))
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
  if (!table_scalars.empty()) {
    build_table_from_scalars(table_scalars);
  } else {
    build_shared_table(m);
  }
  double table_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - table_started).count();
  if (canonical_k3_mode) {
    std::cout << "  table rows stored: " << shared_baby_rows
              << " (canonical merge removed "
              << canonical_rows_merged.load() << ")\n";
  }

  auto search_started = std::chrono::steady_clock::now();
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
    Int source_start(&cover_start);
    uint64_t automatic_chunks = 0;
    while (!found.load(std::memory_order_relaxed) && source_start.IsLower(&cover_end)) {
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
                  << " size: " << chunk_size << "\n";
        std::free(source_hex);
      }
      run_target(shift, chunk_size, true);
      source_start.Add(chunk_size);
      ++automatic_chunks;
    }
    std::cout << "  automatic chunks attempted: " << automatic_chunks << "\n";
    if (!found.load(std::memory_order_relaxed)) {
      Int zero_shift;
      zero_shift.SetInt32(0);
      run_target(zero_shift, range_size, true);
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
  double search_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - search_started).count();
  double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();

  std::cout << "  segments completed: " << segments_done.load() << "\n"
            << "  raw baby rows built: " << raw_baby_rows.load() << "\n"
            << "  unique baby rows built: " << baby_rows.load() << "\n"
            << "  duplicate baby keys: "
            << (raw_baby_rows.load() - baby_rows.load()) << "\n"
            << "  giant steps: " << giant_steps.load() << "\n"
            << "  giant steps scheduled: " << giant_steps_scheduled.load() << "\n"
            << "  invariant hits: " << invariant_hits.load() << "\n"
            << "  orbit factors scheduled: " << orbit_factor_checks.load() << "\n"
            << "  exact candidate checks: " << exact_candidate_checks.load() << "\n"
            << "  EC scalar multiplications: " << ec_scalar_mults.load() << "\n"
            << "  EC point additions: " << ec_point_adds.load() << "\n"
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
    if (canonical_k3_mode) {
      std::cout << "  canonical rows merged: "
                << canonical_rows_merged.load() << "\n";
    }
    if (partner_table_mode) {
      std::cout << "  partner recoveries: "
                << partner_recoveries.load() << "\n";
    }
  }
  std::cout << "  giant steps planned: " << giant_steps_planned.load() << "\n"
            << "  table build seconds: " << table_seconds << "\n"
            << "  giant search seconds: " << search_seconds << "\n"
            << "  giant steps/second: "
            << (search_seconds > 0 ? giant_steps.load() / search_seconds : 0) << "\n"
            << "  total seconds: " << seconds << "\n";
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
