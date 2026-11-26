#!/usr/bin/env python3
"""Delta-sweep for the quantized two-angle fingerprint table.

The reference experiment `parallel_angle_table_recovery.py` measured one
point: W = 100000, 12 bits per angle, mean bucket 1.020, 100/100 exact.
This sweep turns that point into a curve over the grid

    (log2 W, delta)  with  delta = B - ceil(log2 W),  B = fingerprint width.

Physics behind the parameterization: with N records in a 2^B keyspace the
mean bucket occupancy is lambda = N / 2^B, so at fixed delta the whole
bucket statistics is scale-free.  The sweep therefore measures quality
isolines in (log2 W, delta), not (log2 W, B).

Two fingerprint arms share one build pass per W:

  arm "pair"   - the verified two-angle fingerprint, quantized to
                 (b1, b2) = (ceil(B/2), floor(B/2)) bits per angle;
  arm "orbit"  - GLV-orbit canonicalization: the key of a record is the
                 minimum over the three x-distinct orbit pairs
                 (L_r, L_{r+1 mod 3}), r = 0,1,2.  Each record then covers
                 up to 6 scalars (orbit partners lambda^r*k, +/-), the same
                 6->1 idea as the engine --glv-orbit canonicalization.

  arm "sig"    - phase-signature arm, the table-model test of the
                 "--phase-signature shifts the wall left by ~1.5 bits"
                 hypothesis.  (d01, d12) = (L1-L0, L2-L1) mod p are
                 deterministic functions of the same invariant triple, so
                 at fixed total width B the fingerprint is
                 (pair @ B-2 bits) | 2 signature bits: a re-partition of
                 the same budget, not new entropy.  If the wall stays put,
                 the engine's 6->2 gain must come from candidate routing,
                 not from fingerprint information.

The exact invariant values are computed once per scalar at full precision
and quantized down per cell (floor division nests: q_{b-k} = q_b >> k),
so every delta cell of a given W uses exactly the same underlying data.

Metrics per cell: occupied keys, mean/p99/max bucket over occupied keys and
over queries, P(exact single-candidate recovery), Poisson consistency ratio,
estimated deferred-verification time per query, packed and unpacked memory.
"""

import argparse
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from phase_signature_135bit_experiment import BETA, lambda_x
from w1 import FIELD_PRIME, GROUP_ORDER

BETA2 = BETA * BETA % FIELD_PRIME
Q_BUILD_BITS = 18  # quantization stored during the build pass (per angle)


def x_of_scalar(scalar):
    from coincurve import PrivateKey
    return int.from_bytes(
        PrivateKey(scalar.to_bytes(32, "big")).public_key.format(compressed=True)[1:33],
        "big",
    )


def worker(job):
    """Exact (L0, L1, L2) triple per scalar, quantized at Q_BUILD_BITS."""
    start, offset_start, count = job
    q = 1 << Q_BUILD_BITS
    q0 = np.empty(count, dtype=np.uint32)
    q1 = np.empty(count, dtype=np.uint32)
    q2 = np.empty(count, dtype=np.uint32)
    sig = np.empty(count, dtype=np.uint8)  # 2 bits: (s01 << 1) | s12
    for i in range(count):
        x = x_of_scalar(start + offset_start + i)
        l0 = lambda_x(x)
        l1 = lambda_x(BETA * x % FIELD_PRIME)
        l2 = lambda_x(BETA2 * x % FIELD_PRIME)
        q0[i] = l0 * q // FIELD_PRIME
        q1[i] = l1 * q // FIELD_PRIME
        q2[i] = l2 * q // FIELD_PRIME
        d01 = (l1 - l0) % FIELD_PRIME
        d12 = (l2 - l1) % FIELD_PRIME
        sig[i] = ((d01 * 2 // FIELD_PRIME) << 1) | (d12 * 2 // FIELD_PRIME)
    return offset_start, q0, q1, q2, sig


def build_pass(log2_w, workers, workdir):
    """One pass over the window; cached on disk as .npz."""
    path = os.path.join(workdir, f"triple_L{log2_w}.v3.npz")
    if os.path.exists(path):
        with np.load(path) as data:
            return data["q0"], data["q1"], data["q2"], data["sig"]
    w = 1 << log2_w
    start = 1 << 134  # same canonical window as the reference experiment
    q0 = np.empty(w, dtype=np.uint32)
    q1 = np.empty(w, dtype=np.uint32)
    q2 = np.empty(w, dtype=np.uint32)
    sig = np.empty(w, dtype=np.uint8)
    chunk = 1 << 16
    jobs = [(start, o, min(chunk, w - o)) for o in range(0, w, chunk)]
    started = time.perf_counter()
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for offset_start, a0, a1, a2, asig in pool.map(worker, jobs):
            q0[offset_start:offset_start + len(a0)] = a0
            q1[offset_start:offset_start + len(a1)] = a1
            q2[offset_start:offset_start + len(a2)] = a2
            sig[offset_start:offset_start + len(asig)] = asig
            done += 1
            if done % max(1, len(jobs) // 20) == 0 or done == len(jobs):
                rate = (offset_start + len(a0)) / (time.perf_counter() - started)
                print(f"  L={log2_w} chunks {done}/{len(jobs)} "
                      f"({rate/1000:.0f}k scalars/s)", flush=True)
    np.savez(path, q0=q0, q1=q1, q2=q2, sig=sig)
    return q0, q1, q2, sig


def pair_key(q_a, q_b, b1, b2):
    """Quantize stored 16-bit values down to (b1, b2) bits and pack."""
    s1 = Q_BUILD_BITS - b1
    s2 = Q_BUILD_BITS - b2
    return (q_a.astype(np.uint64) >> s1) << np.uint64(b2) | \
           (q_b.astype(np.uint64) >> s2)


def orbit_key(q0, q1, q2, b1, b2):
    """Canonical key: min over the three x-distinct orbit pairs."""
    k0 = pair_key(q0, q1, b1, b2)
    k1 = pair_key(q1, q2, b1, b2)
    k2 = pair_key(q2, q0, b1, b2)
    return np.minimum(np.minimum(k0, k1), k2)


def bucket_stats(keys, query_keys):
    """Run-length statistics of equal-key buckets + query-side metrics."""
    keys = np.sort(keys)
    boundaries = np.flatnonzero(np.diff(keys)) + 1
    starts = np.concatenate(([0], boundaries))
    ends = np.concatenate((boundaries, [len(keys)]))
    counts = (ends - starts).astype(np.float64)
    left = np.searchsorted(keys, query_keys, side="left")
    right = np.searchsorted(keys, query_keys, side="right")
    query_counts = (right - left).astype(np.float64)
    return {
        "keys": keys,
        "occupied": int(len(counts)),
        "mean_bucket_occupied": float(counts.mean()),
        "p99_bucket_occupied": float(np.percentile(counts, 99)),
        "max_bucket": int(counts.max()),
        "frac_singleton": float((counts == 1).mean()),
        "mean_bucket_query": float(query_counts.mean()),
        "p99_bucket_query": float(np.percentile(query_counts, 99)),
        "query_counts": query_counts,
    }


def run_cell(log2_w, delta, arm, triples, queries, rng, verify_cost_s):
    q0, q1, q2, sig = triples
    w = 1 << log2_w
    big_b = log2_w + delta
    if big_b < 2:
        return None
    b1 = (big_b + 1) // 2
    b2 = big_b // 2
    if arm == "pair":
        keys = pair_key(q1, q2, b1, b2)
        query_offsets = rng.choice(w, size=min(queries, w), replace=False)
        query_keys = pair_key(q1[query_offsets], q2[query_offsets], b1, b2)
    elif arm == "orbit":
        keys = orbit_key(q0, q1, q2, b1, b2)
        query_offsets = rng.choice(w, size=min(queries, w), replace=False)
        query_keys = orbit_key(
            q0[query_offsets], q1[query_offsets], q2[query_offsets], b1, b2)
    elif arm == "sig":
        # B-2 bits of the pair + 2 phase-signature bits (same total budget).
        b1 = (big_b - 1) // 2  # pair part gets B-2 bits total
        b2 = (big_b - 2) // 2
        sbits = sig.astype(np.uint64) & np.uint64(3)
        keys = (pair_key(q1, q2, b1, b2) << np.uint64(2)) | sbits
        query_offsets = rng.choice(w, size=min(queries, w), replace=False)
        query_keys = ((pair_key(q1[query_offsets], q2[query_offsets], b1, b2)
                       << np.uint64(2))
                      | (sig[query_offsets].astype(np.uint64) & np.uint64(3)))
    else:
        raise ValueError(f"unknown arm {arm}")
    stats = bucket_stats(keys, query_keys)
    query_counts = stats.pop("query_counts")
    keys = stats.pop("keys")
    p_exact = float((query_counts == 1).mean())
    # Poisson consistency: occupied fraction and multiplicity spread.
    keyspace = float(1 << big_b)
    lam = w / keyspace
    predicted_occupied = keyspace * (1.0 - math.exp(-lam))
    poisson_ratio = stats["occupied"] / predicted_occupied
    packed_bytes = w * (big_b + log2_w) / 8.0
    if arm == "orbit":
        packed_bytes = w * (big_b + log2_w + 2) / 8.0  # + 2 bits orbit index
    unpacked_bytes = w * 16.0
    return {
        "arm": arm,
        "log2_w": log2_w,
        "w": int(w),
        "delta": delta,
        "B": big_b,
        "b1": b1,
        "b2": b2,
        "queries": int(len(query_counts)),
        "p_exact_single": p_exact,
        "mean_bucket_query": stats["mean_bucket_query"],
        "p99_bucket_query": stats["p99_bucket_query"],
        "mean_bucket_occupied": stats["mean_bucket_occupied"],
        "p99_bucket_occupied": stats["p99_bucket_occupied"],
        "max_bucket": stats["max_bucket"],
        "occupied": stats["occupied"],
        "frac_singleton_occupied": stats["frac_singleton"],
        "poisson_occupied_ratio": poisson_ratio,
        "verify_s_per_query": float(stats["mean_bucket_query"] * verify_cost_s),
        "verify_s_total_est": float(w * lam * verify_cost_s),
        "packed_bytes": packed_bytes,
        "packed_bytes_per_covered_scalar":
            packed_bytes / (6.0 * w) if arm == "orbit" else packed_bytes / w,
        "unpacked_bytes": unpacked_bytes,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log2-min", type=int, default=10)
    parser.add_argument("--log2-max", type=int, default=27)
    parser.add_argument("--delta-min", type=int, default=-2)
    parser.add_argument("--delta-max", type=int, default=6)
    parser.add_argument("--queries", type=int, default=200)
    parser.add_argument("--workers", type=int, default=52)
    parser.add_argument("--arms", default="pair,orbit,sig")
    parser.add_argument("--workdir", default="/tmp/angle_delta_sweep")
    parser.add_argument("--output", default="results/angle_table_delta_sweep.json")
    parser.add_argument("--seed", type=int, default=20261126)
    args = parser.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    # Deferred-verification cost per candidate: exact EC reconstruction +
    # full-precision invariant comparison, as in the reference script.
    from parallel_angle_table_recovery import angle_pair, point_from_scalar
    from w1 import DomeInvariantAnalyzer
    analyzer = DomeInvariantAnalyzer()
    t0 = time.perf_counter()
    n_verify = 20
    for i in range(n_verify):
        point = point_from_scalar((1 << 134) + 7919 * i)
        angle_pair(analyzer, point)
    verify_cost_s = (time.perf_counter() - t0) / n_verify
    print(f"deferred verify cost per candidate: {verify_cost_s*1e3:.3f} ms")

    arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    cells = []
    for log2_w in range(args.log2_min, args.log2_max + 1):
        triples = build_pass(log2_w, args.workers, args.workdir)
        for arm in arms:
            for delta in range(args.delta_min, args.delta_max + 1):
                rng = np.random.default_rng(args.seed + delta)
                cell = run_cell(log2_w, delta, arm, triples, args.queries,
                                rng, verify_cost_s)
                if cell is None:
                    continue
                cells.append(cell)
                print(f"L={log2_w:2d} d={delta:+d} {arm:5s} "
                      f"P_exact={cell['p_exact_single']:.3f} "
                      f"bucket_q={cell['mean_bucket_query']:.3f} "
                      f"p99_occ={cell['p99_bucket_occupied']:.0f} "
                      f"poisson={cell['poisson_occupied_ratio']:.4f}", flush=True)
        del triples

    with open(args.output, "w") as handle:
        json.dump({
            "description": __doc__,
            "seed": args.seed,
            "queries_per_cell": args.queries,
            "verify_cost_s_per_candidate": verify_cost_s,
            "canonical_window_start": 1 << 134,
            "cells": cells,
        }, handle, indent=1)
    print(f"saved {len(cells)} cells -> {args.output}")


if __name__ == "__main__":
    main()
