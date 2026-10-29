#!/usr/bin/env python3
"""Recover scalars from public points using a precomputed angle table.

The table covers one known contiguous scalar range.  For every k in that
range it stores the 20-bit first-angle cell and the exact Lambda(P) residue.
For a query public point P, recovery performs only:

    P -> Lambda(P) -> q_20(P) -> table bucket -> exact angle match -> k.

The validation query receives a public point and does not use its scalar until
after recovery has completed.  This is a bounded-range table experiment.
"""

from collections import Counter, defaultdict
import time

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import DomeInvariantAnalyzer, EllipticCurvePoint, FIELD_PRIME, GROUP_ORDER


RANGE_START = 1 << 255
RANGE_SIZE = 100000
ANGLE_BITS = 20
PROGRESS_INTERVAL = 10000


def point_from_scalar(scalar):
    public_numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()
    ).public_key().public_numbers()
    return EllipticCurvePoint(public_numbers.x, public_numbers.y)


def quantized_angle(lambda_value):
    return lambda_value * (1 << ANGLE_BITS) // FIELD_PRIME


def build_table(analyzer, range_start, range_size):
    """Build q_20 -> [(exact Lambda, k)] and retain public queries for validation."""
    table = defaultdict(list)
    queries = []
    for offset in range(range_size):
        scalar = range_start + offset
        point = point_from_scalar(scalar)
        lambda_value = analyzer.compute_torus_invariant_lambda(point)
        table[quantized_angle(lambda_value)].append((lambda_value, scalar))
        queries.append((point, scalar))
        if (offset + 1) % PROGRESS_INTERVAL == 0:
            print(f"  Table: {offset + 1}/{range_size}")
    return table, queries


def recover_scalar(analyzer, point, table):
    """Recover exact table matches using only the supplied public point."""
    lambda_value = analyzer.compute_torus_invariant_lambda(point)
    bucket = table[quantized_angle(lambda_value)]
    exact_matches = [scalar for stored_lambda, scalar in bucket if stored_lambda == lambda_value]
    return bucket, exact_matches


def run_experiment(range_start=RANGE_START, range_size=RANGE_SIZE):
    if range_start < 1 or range_size <= 0 or range_start + range_size > GROUP_ORDER:
        raise ValueError("range must be contained in the canonical scalar interval")

    analyzer = DomeInvariantAnalyzer()
    start_time = time.perf_counter()
    table, queries = build_table(analyzer, range_start, range_size)
    quantized_bucket_sizes = Counter()
    exact_match_sizes = Counter()
    recovered_count = 0

    for index, (point, expected_scalar) in enumerate(queries, start=1):
        bucket, exact_matches = recover_scalar(analyzer, point, table)
        quantized_bucket_sizes[len(bucket)] += 1
        exact_match_sizes[len(exact_matches)] += 1
        recovered_count += exact_matches == [expected_scalar]
        if index % PROGRESS_INTERVAL == 0:
            print(f"  Queries: {index}/{range_size}")

    return {
        "range_start": range_start,
        "range_end": range_start + range_size - 1,
        "range_size": range_size,
        "nonempty_cells": len(table),
        "recovered_count": recovered_count,
        "quantized_bucket_sizes": quantized_bucket_sizes,
        "exact_match_sizes": exact_match_sizes,
        "elapsed_seconds": time.perf_counter() - start_time,
    }


def weighted_average(counter):
    total = sum(counter.values())
    return sum(size * count for size, count in counter.items()) / total


def print_report(report):
    print("Bounded-range angle-table recovery experiment")
    print(f"  Scalar range start: 0x{report['range_start']:064x}")
    print(f"  Scalar range end:   0x{report['range_end']:064x}")
    print(f"  Table scalars: {report['range_size']}")
    print(f"  20-bit angle cells used: {report['nonempty_cells']}")
    print(f"  Exact k recoveries from public P: {report['recovered_count']}/{report['range_size']}")
    print(f"  Mean candidates after 20-bit cell: {weighted_average(report['quantized_bucket_sizes']):.6f}")
    print("  Exact-match candidate counts:")
    for size in sorted(report["exact_match_sizes"]):
        count = report["exact_match_sizes"][size]
        print(f"    {size} candidates: {count}")
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
