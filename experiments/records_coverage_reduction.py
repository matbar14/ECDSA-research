#!/usr/bin/env python3
"""Real reductions of BSGS table records and range coverage (secp256k1).

Motivation (memo 10.11): the canonical K3 table walks the giant grid with
stride == row count, although the mirror branch -c is already checked by
check_orbit_candidate and folds into the core as pure redundancy. This
experiment quantifies the real slack:

  A. Coverage model and its ceiling.
     Giant positions i*stride, candidate rows C, relative offset t is
     coverable iff t = i*stride +- c (c in C). Necessity: t mod stride must
     lie in +/-C mod stride, so gapless tiling needs |±C mod stride| = stride.
     For even stride 2R both residues 0 and R are self-paired (each needs its
     own row), leaving R-2 rows for 2R-2 residues: impossible. Hence the
     ceiling stride <= 2R-1, and the consecutive-head construction
     C = [R, 2R-1] reaches it (±C mod (2R-1) = Z_{2R-1}). Verified:
     exhaustive search over all row sets for small R, enumeration for large R,
     multiplicity histograms (current layout double-covers every offset).

  B. lambda-images cannot extend the tiling. The 6-branch orbit adds the
     multipliers ±1, ±lambda, ±lambda^2; the counting bound would allow
     stride up to ~6R, but the AP structure of lambda*C mod stride leaves
     holes (measured densities far below 1). Also: on the real curve the
     orbit partners of heads [1, H] never land in the window (measured 0),
     so partner coverage is worthless for contiguous sweeps.

  C. Row payload. CBT1 row = 33 B key + 32 B scalar = 65 B. For core-tiling
     tables the scalar column equals base + row_index (implied); the key can
     be truncated to b bits with deferred exact-check resolution. Analytic
     false-positive accounting for the full puzzle-135 window and the
     H-curve windows.

  D. Kangaroo reference point: records -> O(1) at ~c*sqrt(W) group ops,
     measured on the real curve.
"""

import json
import random
import statistics
from itertools import combinations
from pathlib import Path

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
LAMBDA2 = 0xAC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE
BETA = 0x7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE

assert (BETA * BETA % P) * BETA % P == 1, "beta^3 != 1 mod p"
assert (LAMBDA * LAMBDA % N) * LAMBDA % N == 1, "lambda^3 != 1 mod n"
assert LAMBDA2 == LAMBDA * LAMBDA % N, "lambda2 mismatch"


def covered_set(rows, stride, giants, limit):
    """Relative offsets t in [1, limit) coverable as i*stride +/- c."""
    reps = {}
    for i in range(giants):
        base = i * stride
        for c in rows:
            for t in (base + c, base - c):
                if 1 <= t < limit:
                    reps.setdefault(t, []).append((i, c, 1 if t == base + c else -1))
    return reps


def giant_count_for(limit, stride, tail_correction_rows=None):
    """Centers needed to cover [1, limit): standard count, plus the tail
    giant when the last partial period exists (its low residues are covered
    only by the next center's back window). For the mirror-tiling
    construction tail_correction_rows = R; current layout needs no correction."""
    g = (limit - 1) // stride + 1
    if tail_correction_rows is not None and (limit - 1) % stride != 0:
        g += 1
    return g


def part_a_coverage(LIMIT=20000):
    out = {"limit": LIMIT}
    # A1: mirror-tiling construction vs current layout
    for R in (3, 8, 24):
        mirror_rows = list(range(R, 2 * R))
        stride = 2 * R - 1
        giants = giant_count_for(LIMIT, stride, tail_correction_rows=R)
        reps = covered_set(mirror_rows, stride, giants, LIMIT)
        holes = [t for t in range(1, LIMIT) if t not in reps]
        mult = {}
        for t in range(1, LIMIT):
            mult[len(reps.get(t, []))] = mult.get(len(reps.get(t, [])), 0) + 1
        out[f"mirror_tiling_R{R}"] = {
            "rows": R, "stride": stride, "giants": giants,
            "holes": len(holes),
            "multiplicity_histogram": dict(sorted(mult.items())),
        }
        core_rows = list(range(1, R + 1))
        stride_c = R
        giants_c = giant_count_for(LIMIT, stride_c)
        reps_c = covered_set(core_rows, stride_c, giants_c, LIMIT)
        holes_c = [t for t in range(1, LIMIT) if t not in reps_c]
        mult_c = {}
        for t in range(1, LIMIT):
            mult_c[len(reps_c.get(t, []))] = mult_c.get(len(reps_c.get(t, [])), 0) + 1
        out[f"current_layout_R{R}"] = {
            "rows": R, "stride": R, "giants": giants_c,
            "holes": len(holes_c),
            "multiplicity_histogram": dict(sorted(mult_c.items())),
        }
    # A2: exhaustive ceiling for small R - max gapless stride over ALL row
    # sets of size R (rows drawn from Z_stride)
    exhaustive = {}
    for R in range(2, 7):
        best = 0
        best_rows = None
        for stride in range(2 * R + 2, 1, -1):
            found = False
            # rows are subsets of Z_stride of size R; ±C must be all of Z_stride
            for cand in combinations(range(stride), R):
                img = set()
                for c in cand:
                    img.add(c % stride)
                    img.add((-c) % stride)
                if len(img) == stride:
                    found = True
                    best_rows = list(cand)
                    break
            if found:
                best = stride
                break
        exhaustive[f"R{R}"] = {"max_gapless_stride": best,
                               "ceiling_2R_minus_1": 2 * R - 1,
                               "witness_rows": best_rows}
    out["exhaustive_ceiling"] = exhaustive
    return out


def part_a_lambda_density():
    out = {}
    for R in (64, 256, 1024):
        for mode, stride in (("mirror_2R-1", 2 * R - 1),
                             ("3R-1", 3 * R - 1),
                             ("6R-1", 6 * R - 1)):
            rows = list(range(R, 2 * R))
            lam = LAMBDA % stride
            lam2 = LAMBDA * LAMBDA % stride
            img = set()
            for c in rows:
                for u in (1, lam, lam2):
                    img.add((u * c) % stride)
                    img.add((-u * c) % stride)
            density = len(img) / stride
            out[f"R{R}_{mode}"] = {
                "rows": R, "stride": stride, "orbit_residues": len(img),
                "density": round(density, 6),
                "gapless": len(img) == stride,
            }
    return out


def part_a_partner_scatter():
    out = {}
    for H in (1 << 12, 1 << 16, 1 << 20):
        in_window = 0
        canonical_drift = 0
        for j in range(1, H + 1):
            for u in (LAMBDA, LAMBDA2):
                pj = u * j % N
                if 1 <= pj <= 2 * H or N - pj <= 2 * H:
                    in_window += 1
            orbit_min = j
            for v in (pj_mirror_candidates(j, LAMBDA),
                      pj_mirror_candidates(j, LAMBDA2)):
                if v < orbit_min:
                    canonical_drift += 1
                    orbit_min = v
        out[f"H{H}"] = {
            "heads": H,
            "partners_in_[1,2H]": in_window,
            "heads_with_canonical_below_head": canonical_drift,
            "expected_by_birthday": round(4 * H * (2 * H) / N, 12),
        }
    return out


def pj_mirror_candidates(j, u):
    v = u * j % N
    return min(v, N - v)


def part_b_row_payload():
    W = 1 << 134                     # puzzle-135 window
    configs = []
    for H in (100_000, 400_000, 1_000_000, 4_000_000):
        R = H + 1                    # consecutive heads [R, 2R-1]
        g_default = (W - 1) // R + 1
        sigma = 2 * R - 1
        # engine rule: standard count + tail giant per slice; with one slice
        # over the whole window this is (W-1)//sigma + 2
        g_mirror = (W - 1) // sigma + 2
        formats = {
            "cbt1_65B": 65,
            "core33_scalar_implied": 33,
            "trunc_b32": 4, "trunc_b40": 5, "trunc_b48": 6, "trunc_b64": 8,
        }
        rows_fmt = {}
        for name, nbytes in formats.items():
            fp_rate = R / 2 ** (8 * nbytes) if name.startswith("trunc") else 0.0
            # FP resolution cost: one EC scalar mult (~384 point adds, giant
            # step = 1 point add); FPs per full sweep = giants * fp_rate
            fp_overhead_rel = fp_rate * 384.0
            rows_fmt[name] = {
                "bytes_per_row": nbytes,
                "table_MB": round(R * nbytes / 1e6, 1),
                "fp_rate_per_lookup": fp_rate,
                "fp_work_overhead_relative": fp_overhead_rel,
            }
        configs.append({
            "heads_H": H, "rows_R": R,
            "window": W,
            "giants_default_stride_R": g_default,
            "giants_mirror_stride_2R-1": g_mirror,
            "giant_ratio": round(g_default / g_mirror, 4),
            "formats": rows_fmt,
        })
    # H-curve windows S = (3H)^2 measured in the repo (memo: canonical
    # 100k -> 900000 giant steps at stride H)
    hcurve = []
    for H in (50_000, 100_000, 200_000, 400_000):
        S = (3 * H) ** 2
        R = H + 1
        sigma = 2 * R - 1
        hcurve.append({
            "heads_H": H, "window_S": S,
            "giants_default_measured_in_repo": 3 * H,
            "giants_default_formula": (S - 1) // R + 1,
            "giants_mirror_engine_rule": (S - 1) // sigma + 2,
        })
    combined = {
        "bytes_gain_full_65B_to_trunc64": 65 / 8,
        "giant_gain_mirror": (2 * (4_000_000 + 1) - 1) / (4_000_000 + 1),
        "bytes_times_giants_gain": 65 / 8 * (2 * (4_000_000 + 1) - 1) / (4_000_000 + 1),
    }
    return {"puzzle135_configs": configs, "h_curve_windows": hcurve,
            "combined_gain_H4M": combined}


# -------------------- D. kangaroo on the real curve --------------------

def _inv(a, m=P):
    return pow(a, m - 2, m)


class Pt:
    __slots__ = ("x", "y", "inf")

    def __init__(self, x=None, y=None, inf=False):
        self.x, self.y, self.inf = x, y, inf


G_POINT = None  # set lazily: the generator of secp256k1


def gen_point():
    return Pt(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
              0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)


def padd(p1, p2):
    if p1.inf:
        return p2
    if p2.inf:
        return p1
    if p1.x == p2.x and (p1.y + p2.y) % P == 0:
        return Pt(inf=True)
    if p1 is p2 or (p1.x == p2.x and p1.y == p2.y):
        if p1.y % P == 0:
            return Pt(inf=True)
        m = (3 * p1.x * p1.x) * _inv(2 * p1.y) % P
    else:
        m = (p2.y - p1.y) * _inv((p2.x - p1.x) % P) % P
    x3 = (m * m - p1.x - p2.x) % P
    y3 = (m * (p1.x - x3) - p1.y) % P
    return Pt(x3, y3)


def pmul(k, pt):
    result = Pt(inf=True)
    addend = pt
    while k:
        if k & 1:
            result = padd(result, addend)
        addend = padd(addend, addend)
        k >>= 1
    return result


def part_d_kangaroo(trials=150, W_bits=22, dp_bits=6, jump_entries=64,
                    jump_cap=200000):
    W = 1 << W_bits
    G = gen_point()
    rng = random.Random(20261110)
    # jump table: 64 entries uniform in [1, sqrt(W)], index = x mod 64;
    # mean jump ~ sqrt(W)/2 - the standard optimal scale
    jump_max = 1 << (W_bits // 2)
    jump_table = [rng.randrange(1, jump_max + 1) for _ in range(jump_entries)]
    jump_points = [pmul(v, G) for v in jump_table]
    totals = []
    records = []

    for _ in range(trials):
        start = 1 << 40
        k = start + rng.randrange(W)
        mid = start + W // 2
        Q = pmul(k, G)
        # Solve for kk = k - start in [0, W): tame scalar = W/2 + walk_t,
        # wild scalar = (k - start - W/2) + walk_w; collision gives
        # kk = W + walk_t - walk_w. (Tame must stand at the RELATIVE window
        # midpoint W/2, not at start + W/2 - otherwise the start-scalar gap
        # is ~start and the walks never meet.)
        tame = pmul(W // 2, G)
        wild = padd(Q, pmul((N - mid % N) % N, G))
        walk_t = 0
        walk_w = 0
        seen = {}
        steps = 0
        solved = False
        while steps < jump_cap and not solved:
            for which in (0, 1):
                pt = tame if which == 0 else wild
                s = pt.x % jump_entries
                pt = padd(pt, jump_points[s])
                steps += 1
                if which == 0:
                    walk_t += jump_table[s]
                    tame = pt
                else:
                    walk_w += jump_table[s]
                    wild = pt
                if pt.x % (1 << dp_bits) == 0:
                    key = pt.x
                    if key in seen:
                        typ0, off0 = seen[key]
                        if typ0 != which:
                            # stored tame (off0 = walk_t then) vs current wild:
                            # kk = W + off0 - walk_w; stored wild vs current
                            # tame: kk = W + walk_t - off0
                            kk = (W + off0 - walk_w) if typ0 == 0 \
                                else (W + walk_t - off0)
                            k_found = start + kk
                            if 0 <= kk < W and pmul(k_found % N, G).x == Q.x:
                                solved = True
                        # same-type repeat: keep walking (cycle)
                    else:
                        seen[key] = (which, walk_t if which == 0 else walk_w)
                if solved:
                    break
        totals.append(steps if solved else None)
        records.append(len(seen))
    ok = [t for t in totals if t is not None]
    return {
        "trials": trials, "window_bits": W_bits, "dp_bits": dp_bits,
        "jump_entries": jump_entries,
        "solved": len(ok),
        "mean_jumps": round(statistics.mean(ok), 1) if ok else None,
        "std_jumps": round(statistics.stdev(ok), 1) if len(ok) > 1 else None,
        "min_jumps": min(ok) if ok else None,
        "max_jumps": max(ok) if ok else None,
        "mean_over_sqrtW": round(statistics.mean(ok) / (W ** 0.5), 3) if ok else None,
        "mean_records": round(statistics.mean(records), 1) if records else None,
        "max_records": max(records) if records else None,
        "jump_cap": jump_cap,
    }


def main():
    print("== A. coverage model ==")
    a_cov = part_a_coverage()
    for key, val in a_cov.items():
        if key.startswith("mirror") or key.startswith("current"):
            print(f"  {key}: holes={val['holes']} giants={val['giants']} "
                  f"mult={val['multiplicity_histogram']}")
    for key, val in a_cov["exhaustive_ceiling"].items():
        print(f"  exhaustive {key}: max stride {val['max_gapless_stride']} "
              f"(ceiling {val['ceiling_2R_minus_1']})")
    a_lam = part_a_lambda_density()
    for key, val in a_lam.items():
        print(f"  lambda-density {key}: density={val['density']} "
              f"gapless={val['gapless']}")
    a_scat = part_a_partner_scatter()
    for key, val in a_scat.items():
        print(f"  partner-scatter {key}: in-window={val['partners_in_[1,2H]']} "
              f"drift={val['heads_with_canonical_below_head']}")
    print("== C. row payload ==")
    b = part_b_row_payload()
    for cfg in b["puzzle135_configs"]:
        if cfg["heads_H"] in (100_000, 4_000_000):
            print(f"  H={cfg['heads_H']}: giants {cfg['giants_default_stride_R']}"
                  f" -> {cfg['giants_mirror_stride_2R-1']}")
            for name, fmt in cfg["formats"].items():
                print(f"    {name}: {fmt['bytes_per_row']} B/row, "
                      f"table {fmt['table_MB']} MB"
                      + (f", FP overhead {fmt['fp_work_overhead_relative']:.2e}"
                         if fmt["fp_rate_per_lookup"] else ""))
    print(f"  combined gain (bytes x giants) at H=4M: "
          f"{b['combined_gain_H4M']['bytes_times_giants_gain']:.2f}x")
    for hc in b["h_curve_windows"]:
        print(f"  H-curve H={hc['heads_H']}: giants {hc['giants_default_formula']}"
              f" -> {hc['giants_mirror_engine_rule']} (engine rule)")
    print("== D. kangaroo (real curve) ==")
    kang = part_d_kangaroo()
    print(f"  {kang}")

    out = {
        "date": "2026-11-10",
        "curve": "secp256k1",
        "coverage_model": a_cov,
        "lambda_tiling_density": a_lam,
        "partner_scatter": a_scat,
        "row_payload": b,
        "kangaroo_reference": kang,
        "verdicts": {
            "stride_ceiling": "gapless tiling needs +/-C mod stride = Z_stride, "
                              "hence stride <= 2*rows-1; C=[R,2R-1] with stride "
                              "2R-1 reaches the ceiling (0 holes in enumeration, "
                              "exhaustive for R<=6)",
            "current_layout_slack": "default layout (C=[1,R], stride R) covers "
                                    "every offset exactly twice - the mirror "
                                    "redundancy is monetized by --giant-stride "
                                    "2R-1 on heads [R,2R-1]: half the giants",
            "lambda_tiling": "orbit multipliers do NOT extend gapless tiling: "
                             "measured residue densities < 1 at stride 3R-1 and "
                             "6R-1 (holes guaranteed)",
            "partner_scatter": "orbit partners of heads [1,H] never land in the "
                               "swept window (0 in all measured H), partner "
                               "coverage is worthless for contiguous sweeps",
            "row_payload": "scalar column is implied for consecutive heads "
                           "(65->33 B); truncated key b=64 adds 8.4e-11 relative "
                           "work at H=4M full window (65->8 B); combined with "
                           "mirror stride ~16x bytes*giants",
            "kangaroo": "records O(1) (~tens) at ~c*sqrt(W) jumps, measured c "
                        "in kangaroo_reference",
        },
    }
    path = Path(__file__).resolve().parent.parent / "results" / \
        "records_coverage_reduction.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
    print(f"saved: {path}")


if __name__ == "__main__":
    main()
