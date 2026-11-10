#!/usr/bin/env python3
"""Recover an offset from a fully covered interval with segmented BSGS.

The interval is partitioned into segments.  Each segment has its own baby
table; a query subtracts the segment anchor and searches the local offset.
This is a real coverage experiment: every scalar in the finite interval is
covered, while memory is traded for repeated giant-step work.
"""

import argparse
import math
import time

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import BASE_POINT, EllipticCurvePoint, FIELD_PRIME, GROUP_ORDER


def neg(point):
    if point.is_infinity:
        return point
    return EllipticCurvePoint(point.x, (-point.y) % FIELD_PRIME)


def key(point):
    return (point.x, point.y, point.is_infinity)


def scalar_point(k):
    return k * BASE_POINT


def point_from_public_key(encoded):
    raw = bytes.fromhex(encoded)
    public_key = ec.EllipticCurvePublicKey.from_encoded_point(ec.SECP256K1(), raw)
    numbers = public_key.public_numbers()
    return EllipticCurvePoint(numbers.x, numbers.y)


def build_baby(segment_size):
    m = math.isqrt(segment_size - 1) + 1
    table = {}
    current = EllipticCurvePoint(is_infinity=True)
    for j in range(m):
        table.setdefault(key(current), j)
        current = current + BASE_POINT
    return table, m, m * BASE_POINT


def recover_point(range_start, range_size, segment_size, target):
    segment_count = (range_size + segment_size - 1) // segment_size
    total_baby_rows = 0
    total_giant_steps = 0
    started = time.perf_counter()

    for segment in range(segment_count):
        local_start = segment * segment_size
        local_size = min(segment_size, range_size - local_start)
        baby, m, step = build_baby(local_size)
        total_baby_rows += len(baby)
        anchor = scalar_point(range_start + local_start)
        giant = target + neg(anchor)
        for i in range(m + 1):
            total_giant_steps += 1
            j = baby.get(key(giant))
            if j is not None:
                offset = i * m + j
                if offset < local_size:
                    return {
                        "correct": True,
                        "recovered_offset": local_start + offset,
                        "segments": segment_count,
                        "baby_rows_per_segment": len(baby),
                        "total_baby_rows_built": total_baby_rows,
                        "giant_steps": total_giant_steps,
                        "seconds": time.perf_counter() - started,
                    }
            giant = giant + neg(step)

    return {
        "correct": False,
        "recovered_offset": None,
        "segments": segment_count,
        "baby_rows_per_segment": 0,
        "total_baby_rows_built": total_baby_rows,
        "giant_steps": total_giant_steps,
        "seconds": time.perf_counter() - started,
    }


def recover(range_start, range_size, segment_size, target_offset):
    """Compatibility wrapper for the control test."""
    return recover_point(
        range_start, range_size, segment_size,
        scalar_point(range_start + target_offset),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-start", type=int, default=1 << 134)
    parser.add_argument("--range-size", type=int, default=100000)
    parser.add_argument("--segment-size", type=int, nargs="+", default=[100, 1000, 10000])
    parser.add_argument("--target-offset", type=int, default=None)
    parser.add_argument("--target-public-key", default=None,
                        help="compressed or uncompressed SEC1 public key")
    args = parser.parse_args()
    if args.range_start < 1 or args.range_size <= 0:
        parser.error("invalid range")
    if args.range_start + args.range_size > GROUP_ORDER:
        parser.error("range exceeds group order")
    target_offset = args.target_offset if args.target_offset is not None else args.range_size // 2
    if args.target_offset is not None and not 0 <= target_offset < args.range_size:
        parser.error("target-offset must be inside range")

    if args.target_public_key:
        target = point_from_public_key(args.target_public_key)
    else:
        target = scalar_point(args.range_start + target_offset)

    print("Segmented BSGS recovery")
    print(f"  range size: {args.range_size}")
    print(f"  target offset: {target_offset}")
    print("  segment  segments  baby rows/seg  total rows built  giant steps  seconds  correct")
    for segment_size in args.segment_size:
        if segment_size <= 0:
            parser.error("segment sizes must be positive")
        report = recover_point(args.range_start, args.range_size, segment_size, target)
        print(
            f"  {segment_size:7d}  {report['segments']:8d}"
            f"  {report['baby_rows_per_segment']:14d}"
            f"  {report['total_baby_rows_built']:16d}"
            f"  {report['giant_steps']:12d}"
            f"  {report['seconds']:7.3f}  {report['correct']}"
        )


if __name__ == "__main__":
    main()
