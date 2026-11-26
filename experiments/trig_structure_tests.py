#!/usr/bin/env python3
"""Sine-table analogy tests for the C6 / mu_6 structure of secp256k1.

Five checks, in the order of the memo plan:

1. rational_phi  - x(phi P) = beta*x(P) and x(-P) = x(P) verified against
                   real EC arithmetic, then orbit construction benchmarked
                   two ways: three scalar multiplications (EC) versus one
                   multiplication + two field multiplications (x-only
                   rational formulas).  The x-only path is what a baby
                   table should use for the orbit columns.
2. uniformity    - chi^2 of the 12-bit quantized angle over 2^22 window
                   scalars for L0 and L1, plus circular statistics
                   (circular mean / variance via arg(sum e^{i theta}));
                   arithmetic means of angles are never used.
3. orbit_angles  - angular distances between consecutive C6 orbit members
                   in the 12-bit metric; tests whether the metric is
                   phi-equivariant (spacing 1/6 of the turn) or maps the
                   orbit to 3 distinct angles spaced 1/3 (x-only view,
                   +/-P share x and hence the angle).
4. additivity    - 10^6 triples (P, Q, P+Q) built as [k1]G, [k2]G,
                   [k1+k2]G: chi^2 of Delta = a(P+Q)-a(P)-a(Q) mod 2^12
                   and conditional distribution of a(P+Q) inside fixed
                   (a(P), a(Q)) buckets.  Expected: uniform => the angle
                   additivity direction is closed as a negative result.
5. fft           - DFT of the fingerprint series along the scalar window
                   (2^20 consecutive scalars): spectral flatness, power in
                   low-frequency bands.  Flat spectrum = PRF-like metric,
                   no exploitable smoothness at scales > 2^k.

All angle values come from the verified cross-ratio Lambda (w1.py
convention), quantized as value * 2^bits // p.  Cached quantized triples
from angle_table_delta_sweep.py are reused where possible.
"""

import argparse
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from phase_signature_135bit_experiment import LAMBDA, LAMBDA2, lambda_x, BETA
from w1 import FIELD_PRIME

Q_BUILD_BITS = 18
ANGLE_BITS = 12
TURN = 1 << ANGLE_BITS
WORKDIR = "/tmp/angle_delta_sweep"
WINDOW_START = 1 << 134


def x_only_x(scalar):
    from coincurve import PrivateKey
    return int.from_bytes(
        PrivateKey(scalar.to_bytes(32, "big")).public_key.format()[1:33], "big")


def ec_x(scalar):
    from coincurve import PrivateKey
    return int.from_bytes(
        PrivateKey(scalar.to_bytes(32, "big")).public_key.format()[1:33], "big")


def section_rational_phi(trials):
    """Verify rational orbit formulas against EC, benchmark both paths."""
    from coincurve import PrivateKey, PublicKey
    rng = np.random.default_rng(20261126)
    # verification: x(phi P) == beta x(P), x(phi^2 P) == beta^2 x(P),
    # x(-phi^r P) == x(phi^r P)
    scalars = [int(rng.integers(1, 1 << 62)) * (1 << 62) + int(rng.integers(1, 1 << 62)) for _ in range(200)]
    for k in scalars:
        k = int(k)
        x0 = ec_x(k)
        x1 = ec_x((LAMBDA * k) % GROUP_ORDER)
        x2 = ec_x((LAMBDA2 * k) % GROUP_ORDER)
        assert x1 == BETA * x0 % FIELD_PRIME, k
        assert x2 == BETA2_of(x0), k
        assert ec_x((-k) % GROUP_ORDER) == x0, k
        assert ec_x((-LAMBDA * k) % GROUP_ORDER) == x1, k
        assert ec_x((-LAMBDA2 * k) % GROUP_ORDER) == x2, k
    print("rational_phi: x(phi P)=beta*x(P), x(phi^2 P)=beta^2*x(P), "
          "x(-P)=x(P) verified on 200 EC points (all three rotations)")

    # benchmark: orbit x-coordinates (P, phi P, phi^2 P) for `trials` scalars
    keys = [int(rng.integers(1, 1 << 62)) * (1 << 62) + int(rng.integers(1, 1 << 62)) for _ in range(trials)]
    lam_keys = [(LAMBDA * k) % GROUP_ORDER for k in keys]
    lam2_keys = [(LAMBDA2 * k) % GROUP_ORDER for k in keys]
    t0 = time.perf_counter()
    for k, lk, l2k in zip(keys, lam_keys, lam2_keys):
        x0 = ec_x(k)
        x1 = ec_x(lk)
        x2 = ec_x(l2k)
        _ = (x0, x1, x2)
    t_ec = time.perf_counter() - t0
    # x-only path: one multiplication + two field multiplications
    t0 = time.perf_counter()
    for k in keys:
        x0 = ec_x(k)
        x1 = BETA * x0 % FIELD_PRIME
        x2 = BETA * x1 % FIELD_PRIME
        _ = (x0, x1, x2)
    t_xonly = time.perf_counter() - t0
    result = {
        "trials": trials,
        "ec_path_us_per_orbit": t_ec / trials * 1e6,
        "xonly_path_us_per_orbit": t_xonly / trials * 1e6,
        "speedup": t_ec / t_xonly,
    }
    print(f"rational_phi: EC orbit build {result['ec_path_us_per_orbit']:.1f} us, "
          f"x-only {result['xonly_path_us_per_orbit']:.1f} us -> "
          f"{result['speedup']:.2f}x")
    return result


GROUP_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


def BETA2_of(x):
    return BETA * BETA * x % FIELD_PRIME


def load_angle_series(log2_w):
    path = os.path.join(WORKDIR, f"triple_L{log2_w}.v3.npz")
    with np.load(path) as data:
        q0 = (data["q0"].astype(np.uint64) >> (Q_BUILD_BITS - ANGLE_BITS))
        q1 = (data["q1"].astype(np.uint64) >> (Q_BUILD_BITS - ANGLE_BITS))
        q2 = (data["q2"].astype(np.uint64) >> (Q_BUILD_BITS - ANGLE_BITS))
    return q0, q1, q2


def section_uniformity(q0, q1):
    out = {}
    for name, q in (("L0", q0), ("L1", q1)):
        counts = np.bincount(q, minlength=TURN).astype(np.float64)
        n = counts.sum()
        expected = n / TURN
        chi2 = float(((counts - expected) ** 2 / expected).sum())
        dof = TURN - 1
        # circular statistics, never arithmetic means
        unit = np.exp(2j * np.pi * q.astype(np.float64) / TURN)
        z = unit.mean()
        circ_mean = float(np.angle(z) / (2 * np.pi) * TURN % TURN)
        circ_var = float(1 - np.abs(z))  # 0 = concentrated, 1 = uniform
        out[name] = {
            "chi2": chi2, "dof": dof, "chi2_reduced": chi2 / dof,
            "circular_mean_bucket": circ_mean,
            "circular_variance": circ_var,
            "max_bucket": int(counts.max()), "min_bucket": int(counts.min()),
        }
        print(f"uniformity {name}: chi2/dof={chi2/dof:.4f} "
              f"(uniform target 1.0), circular var={circ_var:.6f}, "
              f"circ mean={circ_mean:.1f}/{TURN}")
    return out


def section_orbit_angles(q0, q1, q2):
    """Angular distances between consecutive orbit members, 12-bit metric."""
    dists = []
    for a, b in ((q0, q1), (q1, q2), (q2, q0)):
        d = (b.astype(np.int64) - a.astype(np.int64)) % TURN
        d = np.minimum(d, TURN - d)  # undirected circular distance
        dists.append(d)
    d = np.concatenate(dists)
    sixth = TURN // 6
    third = TURN // 3
    hist6 = np.abs(d - sixth) <= 8
    hist3 = np.abs(d - third) <= 8
    result = {
        "samples": int(len(d)),
        "mean_dist": float(d.mean()),
        "expected_if_sixth": sixth,
        "expected_if_third": third,
        "frac_near_sixth": float(hist6.mean()),
        "frac_near_third": float(hist3.mean()),
        "uniform_ref_frac": (17) / TURN,  # +/-8 window under uniformity
    }
    print(f"orbit_angles: mean undirected distance {d.mean():.1f} "
          f"(1/6={sixth}, 1/3={third}); near-1/6: {result['frac_near_sixth']:.4f}, "
          f"near-1/3: {result['frac_near_third']:.4f}, "
          f"uniform would give {result['uniform_ref_frac']:.4f}")
    return result


def additivity_worker(job):
    seed, count = job
    from coincurve import PrivateKey
    rng = np.random.default_rng(seed)
    k1 = [int(rng.integers(1, 1 << 62)) * (1 << 62) + int(rng.integers(1, 1 << 62))
          for _ in range(count)]
    k2 = [int(rng.integers(1, 1 << 62)) * (1 << 62) + int(rng.integers(1, 1 << 62))
          for _ in range(count)]
    delta = np.empty(count, dtype=np.int64)
    ap = np.empty(count, dtype=np.int64)
    aq = np.empty(count, dtype=np.int64)
    for i in range(count):
        x_p = int.from_bytes(PrivateKey(k1[i].to_bytes(32, "big"))
                             .public_key.format()[1:33], "big")
        x_q = int.from_bytes(PrivateKey(k2[i].to_bytes(32, "big"))
                             .public_key.format()[1:33], "big")
        x_s = int.from_bytes(PrivateKey(((k1[i] + k2[i]) % GROUP_ORDER).to_bytes(32, "big"))
                             .public_key.format()[1:33], "big")
        ap[i] = lambda_x(x_p) * TURN // FIELD_PRIME
        aq[i] = lambda_x(x_q) * TURN // FIELD_PRIME
        delta[i] = ((lambda_x(x_s) * TURN // FIELD_PRIME) - ap[i] - aq[i]) % TURN
    return ap, aq, delta


def section_additivity(trials, workers):
    chunk = 25000
    jobs = [(1000 + i, min(chunk, trials - i * chunk))
            for i in range((trials + chunk - 1) // chunk)]
    aps, aqs, deltas = [], [], []
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for ap, aq, d in pool.map(additivity_worker, jobs):
            aps.append(ap)
            aqs.append(aq)
            deltas.append(d)
    delta = np.concatenate(deltas)
    ap = np.concatenate(aps)
    aq = np.concatenate(aqs)
    elapsed = time.perf_counter() - t0
    counts = np.bincount(delta, minlength=TURN).astype(np.float64)
    expected = len(delta) / TURN
    chi2 = float(((counts - expected) ** 2 / expected).sum())
    # conditional concentration: pair classes (aP//64, aQ//64), each ~244
    # samples at 1e6 triples; aggregate chi2 of a(P+Q) over 2^12 outcomes
    # inside every class, summed over the 4096 classes
    cls = (ap // 64).astype(np.int64) * 64 + (aq // 64)
    chi2_cond = 0.0
    dof_cond = 0
    for c in range(64 * 64):
        sel = delta[cls == c]
        if len(sel) < 20:
            continue
        cc = np.bincount(sel, minlength=TURN).astype(np.float64)
        exp_cc = len(sel) / TURN
        chi2_cond += float(((cc - exp_cc) ** 2 / exp_cc).sum())
        dof_cond += TURN - 1
    result = {
        "trials": int(len(delta)),
        "seconds": elapsed,
        "delta_chi2_reduced": chi2 / (TURN - 1),
        "delta_max_bucket_share": float(counts.max() / len(delta)),
        "cond_chi2_reduced": chi2_cond / dof_cond if dof_cond else None,
        "cond_dof": int(dof_cond),
    }
    print(f"additivity: {len(delta)} triples in {elapsed:.0f}s; "
          f"Delta chi2/dof={chi2/(TURN-1):.4f} (uniform 1.0); "
          f"max bucket share {counts.max()/len(delta):.2e} "
          f"(uniform {1.0/TURN:.2e}); conditional (aP//64,aQ//64) classes: "
          f"chi2/dof={chi2_cond/dof_cond:.4f} over {dof_cond} dof")
    return result


def section_fft(q0):
    series = q0.astype(np.float64) / TURN * 2 * np.pi - np.pi
    n = len(series)
    spec = np.abs(np.fft.rfft(series))
    power = spec ** 2
    low_bins = n // 200 + 1                    # lowest 1% of frequencies
    low = power[1:low_bins].sum()
    result = {
        "n": n,
        "max_over_mean_spec": float(spec.max() / spec[1:].mean()),
        "low_freq_power_share": float(low / power[1:].sum()),
        "uniform_ref_low_share": (low_bins - 1) / (n // 2),
    }
    print(f"fft: n={n}, max/mean spectrum = {result['max_over_mean_spec']:.2f} "
          f"(white noise ~5), lowest-1% band power share = "
          f"{result['low_freq_power_share']:.4%} (uniform "
          f"{result['uniform_ref_low_share']:.4%})")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=1_000_000)
    parser.add_argument("--orbit-bench-trials", type=int, default=20000)
    parser.add_argument("--workers", type=int, default=48)
    parser.add_argument("--output",
                        default="../results/trig_structure_tests.json")
    args = parser.parse_args()

    results = {}
    results["rational_phi"] = section_rational_phi(args.orbit_bench_trials)
    q0, q1, q2 = load_angle_series(22)
    results["uniformity"] = section_uniformity(q0, q1)
    results["orbit_angles"] = section_orbit_angles(q0, q1, q2)
    results["additivity"] = section_additivity(args.trials, args.workers)
    q0_fft, _, _ = load_angle_series(20)
    results["fft"] = section_fft(q0_fft)

    with open(args.output, "w") as handle:
        json.dump(results, handle, indent=1)
    print(f"saved -> {args.output}")


if __name__ == "__main__":
    main()
