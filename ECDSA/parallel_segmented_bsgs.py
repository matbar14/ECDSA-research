#!/usr/bin/env python3
"""Parallel segmented-BSGS throughput and GLV-routing benchmark.

Every segment is an independent worker task.  It builds a local baby table and
searches the segment with giant steps.  With --glv, each worker additionally
routes the six transformed target points through the same segment; this
measures the real cost/benefit of orbit routing rather than assuming a 6x
reduction.  The benchmark covers the selected finite interval completely.
"""

import argparse
import math
import random
import time
from concurrent.futures import ProcessPoolExecutor

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import (
    BASE_POINT,
    ENDOMORPHISM_BETA,
    EllipticCurvePoint,
    FIELD_PRIME,
    GROUP_ORDER,
)


GLV_LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
GLV_FACTORS = (
    1,
    -1,
    GLV_LAMBDA,
    -GLV_LAMBDA,
    pow(GLV_LAMBDA, 2, GROUP_ORDER),
    -pow(GLV_LAMBDA, 2, GROUP_ORDER),
)


def point_key(point):
    return (point.x, point.y, point.is_infinity)


def neg(point):
    if point.is_infinity:
        return point
    return EllipticCurvePoint(point.x, (-point.y) % FIELD_PRIME)


def phi(point):
    if point.is_infinity:
        return point
    return EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)


def scalar_point(k):
    return k * BASE_POINT


def point_from_scalar(k):
    numbers = ec.derive_private_key(k, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(numbers.x, numbers.y)


def target_orbit(point):
    result = []
    current = point
    for _ in range(3):
        result.extend((current, neg(current)))
        current = phi(current)
    return result


def segment_worker(job):
    range_start, local_start, local_size, target_points = job
    m = math.isqrt(local_size - 1) + 1
    baby = {}
    current = EllipticCurvePoint(is_infinity=True)
    for j in range(m):
        baby.setdefault(point_key(current), j)
        current = current + BASE_POINT
    step = m * BASE_POINT
    anchor = scalar_point(range_start + local_start)
    results = []
    giant_steps = 0
    for target_index, target in enumerate(target_points):
        giant = target + neg(anchor)
        for i in range(m + 1):
            giant_steps += 1
            j = baby.get(point_key(giant))
            if j is not None:
                offset = i * m + j
                if offset < local_size:
                    results.append((target_index, local_start + offset))
            giant = giant + neg(step)
    return {
        "segment": local_start,
        "baby_rows": len(baby),
        "giant_steps": giant_steps,
        "matches": results,
    }


def run(range_start, range_size, segment_size, workers, glv, seed):
    target_offset = random.Random(seed).randrange(range_size)
    target_scalar = range_start + target_offset
    target = point_from_scalar(target_scalar)
    target_points = target_orbit(target) if glv else [target]
    jobs = []
    for local_start in range(0, range_size, segment_size):
        jobs.append((
            range_start,
            local_start,
            min(segment_size, range_size - local_start),
            target_points,
        ))

    started = time.perf_counter()
    reports = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for report in pool.map(segment_worker, jobs):
            reports.append(report)
    elapsed = time.perf_counter() - started
    matches = [
        (target_index, offset)
        for report in reports
        for target_index, offset in report["matches"]
    ]
    direct_matches = [offset for target_index, offset in matches if target_index == 0]
    glv_matches = sorted(set(offset for _, offset in matches))
    total_baby = sum(report["baby_rows"] for report in reports)
    total_giant = sum(report["giant_steps"] for report in reports)
    return {
        "target_offset": target_offset,
        "target_scalar": target_scalar,
        "segments": len(jobs),
        "total_baby_rows": total_baby,
        "peak_baby_rows": max(report["baby_rows"] for report in reports),
        "total_giant_steps": total_giant,
        "elapsed": elapsed,
        "giant_steps_per_second": total_giant / elapsed if elapsed else 0.0,
        "direct_matches": sorted(set(direct_matches)),
        "glv_matches": glv_matches,
        "correct": target_offset in direct_matches,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bit-length", type=int, default=135)
    parser.add_argument("--range-size", type=int, default=100000)
    parser.add_argument("--segment-size", type=int, default=10000)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--glv", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.bit_length <= 256 or args.range_size <= 0 or args.segment_size <= 0:
        parser.error("invalid range parameters")
    range_start = 1 << (args.bit_length - 1)
    if range_start + args.range_size > GROUP_ORDER:
        parser.error("range exceeds group order")
    report = run(range_start, args.range_size, args.segment_size,
                 args.workers, args.glv, args.seed)
    print("Parallel segmented BSGS throughput")
    print(f"  interval: [2^{args.bit_length - 1}, 2^{args.bit_length - 1} + {args.range_size})")
    print(f"  segments: {report['segments']}")
    print(f"  workers: {args.workers}")
    print(f"  GLV routing: {args.glv}")
    print(f"  target offset: {report['target_offset']}")
    print(f"  total baby rows: {report['total_baby_rows']}")
    print(f"  peak baby rows/worker: {report['peak_baby_rows']}")
    print(f"  total giant steps: {report['total_giant_steps']}")
    print(f"  elapsed seconds: {report['elapsed']:.3f}")
    print(f"  giant steps/second: {report['giant_steps_per_second']:.1f}")
    print(f"  direct matches: {report['direct_matches']}")
    print(f"  GLV-routed offsets: {report['glv_matches']}")
    print(f"  recovery correct: {report['correct']}")


if __name__ == "__main__":
    main()
