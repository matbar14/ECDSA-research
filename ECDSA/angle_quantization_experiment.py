#!/usr/bin/env python3
"""Measure scalar-range reduction from quantized first-radius angles.

For each scalar k in a known contiguous 256-bit range, compute P = kG and
the exact first-radius index a = Lambda(P).  Quantizing the corresponding
unit-circle angle into b bits is equivalent to

    q_b(k) = floor(2^b * a / p),  0 <= q_b < 2^b.

A precomputed table maps q_b to scalars in the range.  For every tested b the
script reports bucket sizes, the expected number of remaining candidates, and
the scalar span of each bucket.  The latter determines whether a cell could
provide a short contiguous interval for Pollard kangaroo.
"""

import math
import time

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import DomeInvariantAnalyzer, EllipticCurvePoint, FIELD_PRIME, GROUP_ORDER


RANGE_START = 1 << 255
RANGE_SIZE = 100000
QUANTIZATION_BITS = (8, 12, 16, 20)
PROGRESS_INTERVAL = 10000


def point_from_scalar(scalar):
    public_numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()
    ).public_key().public_numbers()
    return EllipticCurvePoint(public_numbers.x, public_numbers.y)


def quantized_angle(lambda_value, bit_count):
    return lambda_value * (1 << bit_count) // FIELD_PRIME


def run_experiment(range_start=RANGE_START, range_size=RANGE_SIZE, collect_buckets=False):
    if not isinstance(range_start, int) or not isinstance(range_size, int):
        raise ValueError("range_start and range_size must be integers")
    if range_start < 1 or range_size <= 0 or range_start + range_size > GROUP_ORDER:
        raise ValueError("range must be contained in the canonical scalar interval")

    analyzer = DomeInvariantAnalyzer()
    # cell -> [number of scalars, first offset, last offset].  The scan is in
    # scalar order, so updating only the last offset is sufficient.
    buckets = {bit_count: {} for bit_count in QUANTIZATION_BITS}
    start_time = time.perf_counter()

    for offset in range(range_size):
        scalar = range_start + offset
        lambda_value = analyzer.compute_torus_invariant_lambda(point_from_scalar(scalar))
        for bit_count, cells in buckets.items():
            cell = quantized_angle(lambda_value, bit_count)
            entry = cells.get(cell)
            if entry is None:
                cells[cell] = [1, offset, offset]
            else:
                entry[0] += 1
                entry[2] = offset
        if (offset + 1) % PROGRESS_INTERVAL == 0:
            print(f"  Completed: {offset + 1}/{range_size}")

    results = []
    for bit_count, cells in buckets.items():
        entries = tuple(cells.values())
        counts = tuple(entry[0] for entry in entries)
        spans = tuple(entry[2] - entry[1] + 1 for entry in entries)
        singleton_keys = sum(count for count in counts if count == 1)
        expected_candidates = math.fsum(count * count for count in counts) / range_size
        weighted_span = math.fsum(
            entry[0] * (entry[2] - entry[1] + 1) for entry in entries
        ) / range_size
        nonsingleton_spans = tuple(
            entry[2] - entry[1] + 1 for entry in entries if entry[0] > 1
        )
        results.append(
            {
                "bits": bit_count,
                "cells": 1 << bit_count,
                "nonempty_cells": len(cells),
                "minimum_bucket": min(counts),
                "maximum_bucket": max(counts),
                "singleton_keys": singleton_keys,
                "singleton_fraction": singleton_keys / range_size,
                "expected_candidates": expected_candidates,
                "weighted_span": weighted_span,
                "weighted_span_fraction": weighted_span / range_size,
                "maximum_span": max(spans),
                "maximum_span_fraction": max(spans) / range_size,
                "nonsingleton_cell_count": len(nonsingleton_spans),
                "mean_nonsingleton_span_fraction": (
                    math.fsum(nonsingleton_spans)
                    / len(nonsingleton_spans)
                    / range_size
                    if nonsingleton_spans
                    else 0.0
                ),
            }
        )

    report = {
        "range_start": range_start,
        "range_end": range_start + range_size - 1,
        "range_size": range_size,
        "results": results,
        "elapsed_seconds": time.perf_counter() - start_time,
    }
    if collect_buckets:
        report["buckets"] = buckets
    return report


def print_report(report):
    print("First-radius angle quantization experiment")
    print(f"  Scalar range start: 0x{report['range_start']:064x}")
    print(f"  Scalar range end:   0x{report['range_end']:064x}")
    print(f"  Scalars in table: {report['range_size']}")
    print("  bits  nonempty cells  bucket range  unique keys  candidates  mean span  max span")
    for result in report["results"]:
        print(
            f"  {result['bits']:4d}"
            f"  {result['nonempty_cells']:14d}"
            f"  {result['minimum_bucket']:4d}..{result['maximum_bucket']:<4d}"
            f"  {result['singleton_fraction']:10.4%}"
            f"  {result['expected_candidates']:10.6f}"
            f"  {result['weighted_span_fraction']:9.4%}"
            f"  {result['maximum_span_fraction']:8.4%}"
        )
        if result["nonsingleton_cell_count"]:
            print(
                "        non-singleton cells: "
                f"{result['nonsingleton_cell_count']}, "
                "mean span: "
                f"{result['mean_nonsingleton_span_fraction']:.4%}"
            )
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
