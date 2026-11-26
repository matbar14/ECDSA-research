/* Exact bounded Pollard-kangaroo benchmark for the local secp256k1 core.
 * This is a measurement harness, not a production key-search mode.
 */
#include <atomic>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <limits>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "secp256k1/SECP256k1.h"

namespace {
constexpr char GROUP_ORDER_HEX[] =
    "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";

Int scalar_mod_order(int64_t value) {
  Int out;
  if (value >= 0) {
    out.SetInt64(static_cast<uint64_t>(value));
    return out;
  }
  Int magnitude;
  magnitude.SetInt64(static_cast<uint64_t>(-(value + 1)) + 1);
  out.SetBase16(GROUP_ORDER_HEX);
  out.Sub(&magnitude);
  return out;
}

struct Key {
  unsigned char bytes[64];
  bool operator==(const Key &o) const { return std::memcmp(bytes, o.bytes, 64) == 0; }
};
struct KeyHash {
  size_t operator()(const Key &k) const {
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char b : k.bytes) { h ^= b; h *= 1099511628211ULL; }
    return static_cast<size_t>(h);
  }
};
struct Record { int kind; uint64_t distance; unsigned factor = 0; };
struct Result {
  bool found = false;
  uint64_t steps = 0, records = 0, scalar = 0, cross_collisions = 0;
  uint64_t orbit_candidates = 0;
};

struct GLV {
  Int beta;
  Int factors[6];

  void init(Secp256K1 &curve) {
    beta.SetBase16("7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE");
    factors[0].SetInt64(1);
    factors[1].SetBase16("5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72");
    factors[2].SetBase16("AC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE");
    for (unsigned i = 0; i < 3; ++i) {
      factors[i + 3].Set(&curve.order);
      factors[i + 3].Sub(&factors[i]);
      factors[i + 3].Mod(&curve.order);
    }
  }
};

Key key_of(Point &p) {
  p.Reduce();
  Key k{};
  p.x.Get32Bytes(k.bytes);
  p.y.Get32Bytes(k.bytes + 32);
  return k;
}

bool distinguished(Point &point, unsigned dp_bits) {
  if (dp_bits == 0) return true;
  for (unsigned bit = 0; bit < dp_bits; ++bit) {
    if (point.x.GetBit(bit) != 0) return false;
  }
  return true;
}

bool distinguished_key(const Key &key, unsigned dp_bits) {
  for (unsigned bit = 0; bit < dp_bits; ++bit) {
    const unsigned byte = 31u - bit / 8;
    const unsigned mask = 1u << (bit % 8);
    if ((key.bytes[byte] & mask) != 0) return false;
  }
  return true;
}

Point glv_transform(Secp256K1 &curve, GLV &glv, Point point, unsigned factor) {
  point.Reduce();
  const unsigned rotation = factor % 3;
  for (unsigned i = 0; i < rotation; ++i) {
    Int rotated_x;
    rotated_x.ModMulK1(&point.x, &glv.beta);
    point.x.Set(&rotated_x);
  }
  if (factor >= 3) point = curve.Negation(point);
  point.Reduce();
  return point;
}

struct OrbitKey {
  Key key;
  unsigned factor = 0;
};

OrbitKey canonical_orbit_key(Secp256K1 &curve, GLV &glv, Point &point) {
  OrbitKey best;
  bool have = false;
  for (unsigned factor = 0; factor < 6; ++factor) {
    Point transformed = glv_transform(curve, glv, point, factor);
    Key key = key_of(transformed);
    if (!have || std::memcmp(key.bytes, best.key.bytes, sizeof(key.bytes)) < 0) {
      best.key = key;
      best.factor = factor;
      have = true;
    }
  }
  return best;
}

unsigned inverse_factor(unsigned factor) {
  static constexpr unsigned inverse[6] = {0, 2, 1, 3, 5, 4};
  return inverse[factor];
}

bool point_matches_scalar(Secp256K1 &curve, Point &target, Int &scalar) {
  Point check = curve.ComputePublicKey(&scalar);
  check.Reduce();
  target.Reduce();
  return check.equals(target);
}

bool glv_transition_check(Secp256K1 &curve, GLV &glv) {
  Int scalar; scalar.SetInt64(7);
  Point source = curve.ComputePublicKey(&scalar);
  Point rotated = glv_transform(curve, glv, source, 1);
  Int expected_scalar;
  expected_scalar.Set(&glv.factors[1]);
  expected_scalar.ModMulK1order(&scalar);
  Point expected = curve.ComputePublicKey(&expected_scalar);
  expected.Reduce();
  if (!rotated.equals(expected)) return false;
  Point cycle = source;
  cycle = glv_transform(curve, glv, cycle, 1);
  cycle = glv_transform(curve, glv, cycle, 1);
  cycle = glv_transform(curve, glv, cycle, 1);
  cycle.Reduce();
  source.Reduce();
  return cycle.equals(source);
}

Result solve(Secp256K1 &curve, uint64_t start, uint64_t width,
             uint64_t target_scalar, unsigned dp_bits, unsigned jump_count,
             uint64_t seed, unsigned step_factor, bool glv_mode) {
  std::mt19937_64 rng(seed);
  GLV glv;
  if (glv_mode) {
    glv.init(curve);
    if (!glv_transition_check(curve, glv)) {
      std::fprintf(stderr, "[E] GLV state transition self-check failed\n");
      return Result{};
    }
  }
  const uint64_t jump_max = uint64_t{1} << (width == 0 ? 1 : (63 - __builtin_clzll(width)) / 2);
  std::vector<uint64_t> jumps(jump_count);
  std::vector<Point> jump_points;
  jump_points.reserve(jump_count);
  for (unsigned i = 0; i < jump_count; ++i) {
    jumps[i] = 1 + rng() % jump_max;
    Int s; s.SetInt64(jumps[i]);
    jump_points.push_back(curve.ComputePublicKey(&s));
  }

  Int target_int; target_int.SetInt64(target_scalar);
  Point target = curve.ComputePublicKey(&target_int);
  Int start_int; start_int.SetInt64(start);
  Point start_point = curve.ComputePublicKey(&start_int);
  Point neg_start = curve.Negation(start_point);
  Point relative_target = curve.Add(target, neg_start);
  relative_target.Reduce();
  Int width_int; width_int.SetInt64(width);
  Point tame = curve.ComputePublicKey(&width_int);
  Point wild = relative_target;
  wild.Reduce();

  // The walk is initialized in relative coordinates. Keep this assertion
  // close to the construction because a sign error makes the collision
  // equations appear to fail much later.
  {
    const int64_t relative_scalar = static_cast<int64_t>(target_scalar - start);
    Int check_scalar = scalar_mod_order(relative_scalar);
    Point expected_wild = curve.ComputePublicKey(&check_scalar);
    if (!wild.equals(expected_wild)) {
      std::fprintf(stderr, "[E] wild initialization scalar relation failed\n");
      return Result{};
    }
  }

  std::unordered_map<Key, Record, KeyHash> records;
  uint64_t tame_dist = 0, wild_dist = 0;
  Result result;
  const uint64_t root = width > 0
      ? (uint64_t{1} << ((63 - __builtin_clzll(width)) / 2)) : 1;
  const uint64_t cap = step_factor > 0 && root <= UINT64_MAX / step_factor
      ? step_factor * root : UINT64_MAX;
  const uint64_t end = start + width;

  auto verify = [&](uint64_t candidate) -> bool {
    if (candidate < start || candidate >= end) return false;
    Int s; s.SetInt64(candidate);
    Point check = curve.ComputePublicKey(&s);
    check.Reduce();
    target.Reduce();
    return check.equals(target);
  };
  auto verify_int = [&](Int &candidate) -> bool {
    Int low; low.SetInt64(start);
    Int high; high.SetInt64(end);
    if (candidate.IsLower(&low) || !candidate.IsLower(&high)) return false;
    return point_matches_scalar(curve, target, candidate);
  };
  auto orbit_candidate = [&](unsigned tame_factor, unsigned wild_factor,
                             uint64_t tame_distance,
                             uint64_t wild_distance) {
    Int ratio;
    ratio.Set(&glv.factors[inverse_factor(wild_factor)]);
    ratio.ModMulK1order(&glv.factors[tame_factor]);
    Int tame_scalar; tame_scalar.SetInt64(width + tame_distance);
    Int term;
    term.Set(&ratio);
    term.ModMulK1order(&tame_scalar);
    Int candidate; candidate.SetInt64(start);
    candidate.Mod(&curve.order);
    Int wild_offset; wild_offset.SetInt64(wild_distance);
    if (candidate.IsLower(&wild_offset)) candidate.Add(&curve.order);
    candidate.Sub(&wild_offset);
    candidate.Mod(&curve.order);
    candidate.Add(&term);
    candidate.Mod(&curve.order);
    return candidate;
  };
  auto visit = [&](Point &point, int kind, uint64_t distance) -> uint64_t {
    OrbitKey orbit;
    Key k;
    unsigned factor = 0;
    if (glv_mode) {
      orbit = canonical_orbit_key(curve, glv, point);
      k = orbit.key;
      factor = orbit.factor;
      if (!distinguished_key(k, dp_bits)) return 0;
    } else {
      k = key_of(point);
      if (!distinguished(point, dp_bits)) return 0;
    }
    auto same = records.find(k);
    if (same != records.end() && same->second.kind != kind) {
      ++result.cross_collisions;
      uint64_t td = kind == 0 ? distance : same->second.distance;
      uint64_t wd = kind == 1 ? distance : same->second.distance;
      if (glv_mode) {
        const unsigned tf = kind == 0 ? factor : same->second.factor;
        const unsigned wf = kind == 1 ? factor : same->second.factor;
        Int candidate = orbit_candidate(tf, wf, td, wd);
        ++result.orbit_candidates;
        if (verify_int(candidate)) return candidate.GetInt64();
      }
      const __int128 candidate = static_cast<__int128>(start) + width + td - wd;
      if (!glv_mode && candidate >= static_cast<__int128>(start) &&
          candidate < static_cast<__int128>(end) &&
          verify(static_cast<uint64_t>(candidate))) {
        return static_cast<uint64_t>(candidate);
      }
    }
    if (!glv_mode) {
      Point np = curve.Negation(point);
      Key nk = key_of(np);
      auto other = records.find(nk);
      if (other != records.end() && other->second.kind != kind) {
        ++result.cross_collisions;
        uint64_t td = kind == 0 ? distance : other->second.distance;
        uint64_t wd = kind == 1 ? distance : other->second.distance;
        // Opposite-point collision: relative target = -tame_dist-wild_dist.
        const __int128 candidate = static_cast<__int128>(start) - width - td - wd;
        if (candidate >= static_cast<__int128>(start) &&
            candidate < static_cast<__int128>(end) &&
            verify(static_cast<uint64_t>(candidate))) {
          return static_cast<uint64_t>(candidate);
        }
      }
    }
    records.emplace(k, Record{kind, distance, factor});
    return 0;
  };

  while (result.steps < cap) {
    for (int kind = 0; kind < 2; ++kind) {
      Point &point = kind == 0 ? tame : wild;
      unsigned index = static_cast<unsigned>(point.x.GetInt64() % jump_count);
      point = curve.Add(point, jump_points[index]);
      uint64_t &distance = kind == 0 ? tame_dist : wild_dist;
      distance += jumps[index];
      ++result.steps;
      uint64_t found = visit(point, kind, distance);
      if (found != 0) {
        result.found = true; result.scalar = found;
        result.records = records.size(); return result;
      }
    }
  }
  result.records = records.size();
  return result;
}
}

int main(int argc, char **argv) {
  unsigned bits = 20, trials = 16, threads = 1, dp_bits = 6, jumps = 64;
  unsigned step_factor = 1000;
  uint64_t start = uint64_t{1} << 40;
  uint64_t seed = 0x20261118ULL;
  bool glv_mode = false;
  static option opts[] = {
    {"bits", required_argument, nullptr, 'b'}, {"trials", required_argument, nullptr, 't'},
    {"threads", required_argument, nullptr, 'T'}, {"dp-bits", required_argument, nullptr, 'd'},
    {"jumps", required_argument, nullptr, 'j'}, {"start", required_argument, nullptr, 's'},
    {"seed", required_argument, nullptr, 'S'}, {"step-factor", required_argument, nullptr, 'f'},
    {"glv-orbit", no_argument, nullptr, 'g'},
    {nullptr, 0, nullptr, 0}};
  int c;
  while ((c = getopt_long(argc, argv, "b:t:T:d:j:s:S:f:g", opts, nullptr)) != -1) {
    try {
      switch (c) {
        case 'b': bits = static_cast<unsigned>(std::stoul(optarg)); break;
        case 't': trials = static_cast<unsigned>(std::stoul(optarg)); break;
        case 'T': threads = static_cast<unsigned>(std::stoul(optarg)); break;
        case 'd': dp_bits = static_cast<unsigned>(std::stoul(optarg)); break;
        case 'j': jumps = static_cast<unsigned>(std::stoul(optarg)); break;
        case 's': start = std::stoull(optarg, nullptr, 0); break;
        case 'S': seed = std::stoull(optarg, nullptr, 0); break;
        case 'f': step_factor = static_cast<unsigned>(std::stoul(optarg)); break;
        case 'g': glv_mode = true; break;
        default: return 2;
      }
    } catch (const std::exception &) {
      std::fprintf(stderr, "[E] invalid numeric option\n");
      return 2;
    }
  }
  if (bits == 0 || bits >= 63 || start == 0 || jumps == 0 ||
      threads == 0 || dp_bits > 256) {
    std::fprintf(stderr, "[E] require 1 <= bits < 63, start > 0, "
                 "jumps/threads > 0, dp-bits <= 256\n");
    return 2;
  }
  const uint64_t width = uint64_t{1} << bits;
  if (start > UINT64_MAX - width) {
    std::fprintf(stderr, "[E] start + width overflows uint64\n");
    return 2;
  }
  std::atomic<unsigned> next{0}; std::atomic<unsigned> solved{0};
  std::atomic<uint64_t> total_steps{0}, total_records{0};
  std::atomic<uint64_t> total_cross_collisions{0}, total_orbit_candidates{0};
  auto worker = [&](unsigned worker_id) {
    Secp256K1 curve; curve.Init();
    if (worker_id == 0) {
      Int a; a.SetInt64(start + 7); Point pa = curve.ComputePublicKey(&a);
      Int b; b.SetInt64(start); Point pb = curve.ComputePublicKey(&b);
      Point nb = curve.Negation(pb); Point diff = curve.Add(pa, nb); diff.Reduce();
      Int seven; seven.SetInt64(7); Point expected = curve.ComputePublicKey(&seven); expected.Reduce();
      if (!diff.equals(expected)) {
        std::fprintf(stderr, "[E] point subtraction self-check failed\n");
        return;
      }
    }
    for (;;) {
      unsigned trial = next.fetch_add(1); if (trial >= trials) break;
      std::mt19937_64 rng(seed + static_cast<uint64_t>(trial) * 1009);
      uint64_t target = start + rng() % width;
      Result r = solve(curve, start, width, target, dp_bits, jumps,
                       rng(), step_factor, glv_mode);
      if (r.found && r.scalar == target) solved.fetch_add(1);
      total_steps.fetch_add(r.steps); total_records.fetch_add(r.records);
      total_cross_collisions.fetch_add(r.cross_collisions);
      total_orbit_candidates.fetch_add(r.orbit_candidates);
    }
  };
  auto t0 = std::chrono::steady_clock::now();
  std::vector<std::thread> pool; for (unsigned i = 0; i < threads; ++i) pool.emplace_back(worker, i);
  for (auto &thread : pool) thread.join();
  double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("mode: %s\nwidth bits: %u\ntrials: %u\nthreads: %u\n"
              "solved: %u/%u\nmean steps: %.2f\nmean steps/sqrt(W): %.4f\n"
              "mean DP records: %.2f\nmean cross collisions: %.2f\n"
              "mean orbit candidate checks: %.2f\nwall seconds: %.6f\n"
              "steps/second: %.2f\n",
              glv_mode ? "glv-orbit-cpp" : "exact-negation-aware-cpp",
              bits, trials, threads, solved.load(), trials,
              trials ? double(total_steps.load()) / trials : 0.0,
              trials ? double(total_steps.load()) / trials / (uint64_t{1} << (bits / 2)) : 0.0,
              trials ? double(total_records.load()) / trials : 0.0,
              trials ? double(total_cross_collisions.load()) / trials : 0.0,
              trials ? double(total_orbit_candidates.load()) / trials : 0.0,
              seconds,
              seconds > 0 ? double(total_steps.load()) / seconds : 0.0);
  return solved == trials ? 0 : 1;
}
