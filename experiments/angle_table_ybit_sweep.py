#!/usr/bin/env python3
"""Y-bit arm: the last independent fingerprint bit is sign(y).

The sig arm proved that bits derived from the invariant triple
(L0, L1, L2) add zero entropy beyond the stored high bits.  By the
mu_6-invariance argument the only scalar-dependent information outside
the x-derived triple is the sign of y.  This script adds exactly one
sign(y) bit to the pair fingerprint at a fixed budget

    key = (pair @ (B-1) bits) << 1 | sign(y),   B = L + delta

and checks the prediction: the y-bit is independent of the stored bits,
so ybit(delta) must reproduce pair(delta) exactly (bucket law
1 + 2^-delta, P(exact) = e^-2^-delta) - a full 1-bit gain over the
coarser pair@(B-1) part, and not one bit more.

Reuses the cached 18-bit quantized triples (triple_L*.v3.npz) built by
angle_table_delta_sweep.py; only the y parity is newly computed (the
compressed encoding already carries it, so the pass is derive-only).
"""

import argparse
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

Q_BUILD_BITS = 18


def y_parity_worker(job):
    start, offset_start, count = job
    from coincurve import PrivateKey
    y = np.empty(count, dtype=np.uint8)
    for i in range(count):
        k = start + offset_start + i
        y[i] = PrivateKey(k.to_bytes(32, "big")).public_key.format()[0] & 1
    return offset_start, y


def build_y_pass(log2_w, workers, workdir):
    path = os.path.join(workdir, f"yparity_L{log2_w}.v1.npz")
    if os.path.exists(path):
        with np.load(path) as data:
            return data["y"]
    w = 1 << log2_w
    start = 1 << 134
    yarr = np.empty(w, dtype=np.uint8)
    chunk = 1 << 16
    jobs = [(start, o, min(chunk, w - o)) for o in range(0, w, chunk)]
    started = time.perf_counter()
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for offset_start, part in pool.map(y_parity_worker, jobs):
            yarr[offset_start:offset_start + len(part)] = part
            done += 1
            if done % max(1, len(jobs) // 10) == 0 or done == len(jobs):
                rate = (offset_start + len(part)) / (time.perf_counter() - started)
                print(f"  y L={log2_w} chunks {done}/{len(jobs)} "
                      f"({rate/1000:.0f}k scalars/s)", flush=True)
    np.savez(path, y=yarr)
    return yarr


def pair_key(q_a, q_b, b1, b2):
    s1 = Q_BUILD_BITS - b1
    s2 = Q_BUILD_BITS - b2
    return (q_a.astype(np.uint64) >> s1) << np.uint64(b2) | \
           (q_b.astype(np.uint64) >> s2)


def bucket_stats(keys, query_keys):
    keys = np.sort(keys)
    boundaries = np.flatnonzero(np.diff(keys)) + 1
    starts = np.concatenate(([0], boundaries))
    ends = np.concatenate((boundaries, [len(keys)]))
    counts = (ends - starts).astype(np.float64)
    left = np.searchsorted(keys, query_keys, side="left")
    right = np.searchsorted(keys, query_keys, side="right")
    query_counts = (right - left).astype(np.float64)
    return counts, query_counts


def run_cell(log2_w, delta, triples, yarr, queries, rng, verify_cost_s):
    q0, q1, q2 = triples
    w = 1 << log2_w
    big_b = log2_w + delta
    b1 = (big_b // 2)          # pair part gets B-1 bits
    b2 = (big_b - 1) // 2
    keys = (pair_key(q1, q2, b1, b2) << np.uint64(1)) | yarr.astype(np.uint64)
    query_offsets = rng.choice(w, size=min(queries, w), replace=False)
    query_keys = ((pair_key(q1[query_offsets], q2[query_offsets], b1, b2)
                   << np.uint64(1)) | yarr[query_offsets].astype(np.uint64))
    counts, query_counts = bucket_stats(keys, query_keys)
    keyspace = float(1 << big_b)
    lam = w / keyspace
    predicted_occupied = keyspace * (1.0 - math.exp(-lam))
    return {
        "arm": "ybit",
        "log2_w": log2_w,
        "w": int(w),
        "delta": delta,
        "B": big_b,
        "b1": b1,
        "b2": b2,
        "queries": int(len(query_counts)),
        "p_exact_single": float((query_counts == 1).mean()),
        "mean_bucket_query": float(query_counts.mean()),
        "p99_bucket_query": float(np.percentile(query_counts, 99)),
        "mean_bucket_occupied": float(counts.mean()),
        "p99_bucket_occupied": float(np.percentile(counts, 99)),
        "max_bucket": int(counts.max()),
        "occupied": int(len(counts)),
        "frac_singleton_occupied": float((counts == 1).mean()),
        "poisson_occupied_ratio": float(len(counts) / predicted_occupied),
        "verify_s_per_query": float(query_counts.mean() * verify_cost_s),
        "packed_bytes": w * (big_b + log2_w) / 8.0,
        "unpacked_bytes": w * 16.0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log2-min", type=int, default=10)
    parser.add_argument("--log2-max", type=int, default=27)
    parser.add_argument("--delta-min", type=int, default=-2)
    parser.add_argument("--delta-max", type=int, default=6)
    parser.add_argument("--queries", type=int, default=200)
    parser.add_argument("--workers", type=int, default=56)
    parser.add_argument("--workdir", default="/tmp/angle_delta_sweep")
    parser.add_argument("--output", default="../results/angle_table_ybit_sweep.json")
    parser.add_argument("--seed", type=int, default=20261126)
    args = parser.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    verify_cost_s = 0.0007761108499835245  # measured in the delta sweep

    cells = []
    for log2_w in range(args.log2_min, args.log2_max + 1):
        qpath = os.path.join(args.workdir, f"triple_L{log2_w}.v3.npz")
        with np.load(qpath) as data:
            triples = (data["q0"], data["q1"], data["q2"])
        yarr = build_y_pass(log2_w, args.workers, args.workdir)
        for delta in range(args.delta_min, args.delta_max + 1):
            rng = np.random.default_rng(args.seed + delta)
            cell = run_cell(log2_w, delta, triples, yarr, args.queries,
                            rng, verify_cost_s)
            cells.append(cell)
            print(f"L={log2_w:2d} d={delta:+d} ybit "
                  f"P_exact={cell['p_exact_single']:.3f} "
                  f"bucket_q={cell['mean_bucket_query']:.3f} "
                  f"poisson={cell['poisson_occupied_ratio']:.4f}", flush=True)

    with open(args.output, "w") as handle:
        json.dump({
            "description": __doc__,
            "seed": args.seed,
            "queries_per_cell": args.queries,
            "cells": cells,
        }, handle, indent=1)
    print(f"saved {len(cells)} cells -> {args.output}")


if __name__ == "__main__":
    main()
