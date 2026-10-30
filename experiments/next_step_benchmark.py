#!/usr/bin/env python3
"""Compare direct angle tables, GLV deduplication, and BSGS on finite windows.

This is deliberately a bounded benchmark.  It measures actual rows in the
selected interval; it does not extrapolate a local table to the whole 2^135
space without reporting the scale.
"""

import argparse
import math
import random
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import (
    BASE_POINT,
    DomeInvariantAnalyzer,
    EllipticCurvePoint,
    ENDOMORPHISM_BETA,
    FIELD_PRIME,
)


def point_from_scalar(k):
    q = ec.derive_private_key(k, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(q.x, q.y)


def orbit_points(point):
    current = point
    result = []
    for _ in range(3):
        result.append(current)
        current = EllipticCurvePoint(
            ENDOMORPHISM_BETA * current.x % FIELD_PRIME, current.y
        )
    return result


def quantized_signature(analyzer, point, bits):
    values = [analyzer.compute_torus_invariant_lambda(q) for q in orbit_points(point)]
    return tuple(sorted(value * (1 << bits) // FIELD_PRIME for value in values))


def signature_worker(args):
    start, first, count, bits = args
    analyzer = DomeInvariantAnalyzer()
    result = []
    for offset in range(first, first + count):
        scalar = start + offset
        result.append((quantized_signature(analyzer, point_from_scalar(scalar), bits), offset))
    return result


def neg(point):
    if point.is_infinity:
        return point
    return EllipticCurvePoint(point.x, (-point.y) % FIELD_PRIME)


def point_key(point):
    return (point.x, point.y, point.is_infinity)


def bsgs_recover(range_start, range_size, target_offset):
    """Classical interval BSGS; returns recovered offset and table size."""
    m = math.isqrt(range_size - 1) + 1
    baby = {}
    current = EllipticCurvePoint(is_infinity=True)
    for j in range(m):
        baby[point_key(current)] = j
        current = current + BASE_POINT
    step = m * BASE_POINT
    target = point_from_scalar(range_start + target_offset) + neg(
        point_from_scalar(range_start)
    )
    giant = target
    for i in range(m + 1):
        j = baby.get(point_key(giant))
        if j is not None:
            offset = i * m + j
            if offset < range_size:
                return offset, m
        giant = giant + neg(step)
    return None, m


def run_one(range_start, range_size, bits, workers, chunk_size, seed):
    jobs = [
        (range_start, first, min(chunk_size, range_size - first), bits)
        for first in range(0, range_size, chunk_size)
    ]
    started = time.perf_counter()
    all_records = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for chunk in pool.map(signature_worker, jobs):
            all_records.extend(chunk)
    signature_seconds = time.perf_counter() - started

    direct_rows = range_size
    glv_rows = len({signature for signature, _ in all_records})
    target_offset = random.Random(seed).randrange(range_size)
    bsgs_started = time.perf_counter()
    recovered, bsgs_rows = bsgs_recover(range_start, range_size, target_offset)
    bsgs_seconds = time.perf_counter() - bsgs_started

    return {
        "range_size": range_size,
        "direct_rows": direct_rows,
        "glv_rows": glv_rows,
        "glv_ratio": direct_rows / glv_rows,
        "fingerprint_seconds": signature_seconds,
        "bsgs_rows": bsgs_rows,
        "bsgs_memory_bytes": bsgs_rows * (33 + 8),
        "bsgs_seconds": bsgs_seconds,
        "target_offset": target_offset,
        "bsgs_recovered": recovered,
        "bsgs_correct": recovered == target_offset,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-start", type=int, default=1 << 134)
    parser.add_argument("--sizes", type=int, nargs="+", default=[10000, 100000])
    parser.add_argument("--bits", type=int, default=12)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    for size in args.sizes:
        if size <= 0:
            parser.error("sizes must be positive")
        report = run_one(args.range_start, size, args.bits, args.workers,
                         args.chunk_size, args.seed)
        print(f"window={size}")
        print(f"  direct rows: {report['direct_rows']}")
        print(f"  GLV signature rows in window: {report['glv_rows']}")
        print(f"  actual GLV row reduction: {report['glv_ratio']:.6f}x")
        print(f"  signature build seconds: {report['fingerprint_seconds']:.2f}")
        print(f"  classical BSGS rows: {report['bsgs_rows']}")
        print(f"  BSGS estimated bytes: {report['bsgs_memory_bytes']}")
        print(f"  BSGS recovery seconds: {report['bsgs_seconds']:.3f}")
        print(f"  BSGS recovery correct: {report['bsgs_correct']}")


if __name__ == "__main__":
    main()
