#!/usr/bin/env python3
"""Bounded-interval BSGS with six-way GLV orbit canonicalization.

The table stores one canonical representative for each baby-step orbit.  On a
giant-step collision all six scalar multipliers are tested and the candidate
is verified in the original interval.  The benchmark reports actual
deduplication, not an assumed factor of six.
"""

import argparse
import math
import random
import time

from w1 import (
    BASE_POINT,
    ENDOMORPHISM_BETA,
    EllipticCurvePoint,
    FIELD_PRIME,
    GROUP_ORDER,
)


GLV_LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
GLV_FACTORS = (1, -1, GLV_LAMBDA, -GLV_LAMBDA,
               (GLV_LAMBDA * GLV_LAMBDA) % GROUP_ORDER,
               -(GLV_LAMBDA * GLV_LAMBDA) % GROUP_ORDER)


def point_key(point):
    if point.is_infinity:
        return (0, 0, 1)
    return (point.x, point.y, 0)


def neg(point):
    if point.is_infinity:
        return point
    return EllipticCurvePoint(point.x, (-point.y) % FIELD_PRIME)


def phi(point):
    if point.is_infinity:
        return point
    return EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)


def canonical_orbit_key(point):
    current = point
    keys = []
    for _ in range(3):
        keys.append(point_key(current))
        keys.append(point_key(neg(current)))
        current = phi(current)
    return min(keys)


def scalar_point(scalar):
    return scalar * BASE_POINT


def recover_glv(range_start, range_size, target_offset):
    m = math.isqrt(range_size - 1) + 1
    table = {}
    current = EllipticCurvePoint(is_infinity=True)
    started = time.perf_counter()
    for j in range(m):
        table.setdefault(canonical_orbit_key(current), []).append(j)
        current = current + BASE_POINT

    target = scalar_point(range_start + target_offset)
    base = scalar_point(range_start)
    giant = target + neg(base)
    giant_step = scalar_point(m)
    candidates = set()
    for i in range(m + 1):
        entries = table.get(canonical_orbit_key(giant), ())
        for j in entries:
            for factor in GLV_FACTORS:
                candidate = i * m + (factor * j) % GROUP_ORDER
                if 0 <= candidate < range_size:
                    candidate_point = scalar_point(range_start + candidate)
                    if point_key(candidate_point) == point_key(target):
                        candidates.add(candidate)
        giant = giant + neg(giant_step)

    return {
        "m": m,
        "raw_rows": m,
        "canonical_rows": len(table),
        "dedup_ratio": m / len(table),
        "candidates": sorted(candidates),
        "target_offset": target_offset,
        "build_query_seconds": time.perf_counter() - started,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-start", type=int, default=1 << 134)
    parser.add_argument("--sizes", type=int, nargs="+", default=[10000, 100000, 250000])
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()

    for size in args.sizes:
        if size <= 0:
            parser.error("sizes must be positive")
        target_offset = random.Random(args.seed + size).randrange(size)
        result = recover_glv(args.range_start, size, target_offset)
        print(f"window={size}")
        print(f"  target offset: {target_offset}")
        print(f"  raw baby rows: {result['raw_rows']}")
        print(f"  canonical GLV rows: {result['canonical_rows']}")
        print(f"  actual GLV table reduction: {result['dedup_ratio']:.6f}x")
        print(f"  candidates: {result['candidates']}")
        print(f"  recovery correct: {result['candidates'] == [target_offset]}")
        print(f"  seconds: {result['build_query_seconds']:.3f}")


if __name__ == "__main__":
    main()
