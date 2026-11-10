#!/usr/bin/env python3
"""Test table lookup by Lambda-triangle area on secp256k1.

This experiment does not assume that an area identifies a private scalar.  It
constructs a table containing pairs k and n-k.  Their public points are P and
-P, so their x coordinates, the Lambda orbit values, and the resulting phase
triangle area are equal.  A lookup therefore returns at least two candidates.

The finite-field Lambda computations are exact.  The Euclidean area is encoded
as the exact IEEE-754 binary64 result used by the selected phase model; equality
of the paired Lambda triples is checked separately and is exact.
"""

from collections import defaultdict
import random
import time

from lambda_circle_triangle_experiment import (
    GROUP_ORDER,
    lambda_orbit,
    point_from_scalar,
    triangle_area,
    unit_circle_vertices,
)
from w1 import DomeInvariantAnalyzer, EllipticCurvePoint


TABLE_KEY_COUNT = 100000
SEED = 20260922
PROGRESS_INTERVAL = 10000


def generate_256_bit_scalar(generator):
    return generator.randrange(1 << 255, GROUP_ORDER)


def area_record(analyzer, scalar, point):
    """Compute both the exact Lambda triple and the chosen numeric area key."""
    values = lambda_orbit(analyzer, point)
    if values is None:
        return None
    area = triangle_area(unit_circle_vertices(values))
    return {
        "scalar": scalar,
        "lambda_values": values,
        "area": area,
        "area_key": area.hex(),
    }


def run_experiment(table_key_count=TABLE_KEY_COUNT, seed=SEED):
    """Build a collision-aware area table and query every inserted key."""
    if not isinstance(table_key_count, int) or table_key_count <= 0:
        raise ValueError("table_key_count must be a positive integer")
    if table_key_count % 2:
        raise ValueError("table_key_count must be even so every key has a negation pair")

    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    table = defaultdict(list)
    records = []
    exact_lambda_pair_mismatches = 0
    area_pair_mismatches = 0
    start_time = time.perf_counter()

    for pair_index in range(1, table_key_count // 2 + 1):
        scalar = generate_256_bit_scalar(generator)
        point = point_from_scalar(scalar)
        negated_scalar = (-scalar) % GROUP_ORDER
        negated_point = EllipticCurvePoint(point.x, -point.y)

        positive_record = area_record(analyzer, scalar, point)
        negative_record = area_record(analyzer, negated_scalar, negated_point)
        if positive_record is None or negative_record is None:
            raise RuntimeError("Lambda is undefined for a generated table entry")

        exact_lambda_pair_mismatches += (
            positive_record["lambda_values"] != negative_record["lambda_values"]
        )
        area_pair_mismatches += (
            positive_record["area_key"] != negative_record["area_key"]
        )

        for record in (positive_record, negative_record):
            table[record["area_key"]].append(record["scalar"])
            records.append(record)

        if pair_index % (PROGRESS_INTERVAL // 2) == 0:
            print(f"  Completed: {pair_index * 2}/{table_key_count} table keys")

    lookup_successes = 0
    lookup_candidate_counts = defaultdict(int)
    for record in records:
        candidates = table[record["area_key"]]
        lookup_successes += record["scalar"] in candidates
        lookup_candidate_counts[len(candidates)] += 1

    collision_buckets = {
        key: candidates for key, candidates in table.items() if len(candidates) > 1
    }
    return {
        "table_key_count": table_key_count,
        "seed": seed,
        "exact_lambda_pair_mismatches": exact_lambda_pair_mismatches,
        "area_pair_mismatches": area_pair_mismatches,
        "distinct_area_keys": len(table),
        "collision_buckets": collision_buckets,
        "lookup_successes": lookup_successes,
        "lookup_candidate_counts": lookup_candidate_counts,
        "elapsed_seconds": time.perf_counter() - start_time,
    }


def print_report(report):
    print("Area lookup table experiment")
    print(f"  Table keys: {report['table_key_count']}")
    print(f"  Reproducibility seed: {report['seed']}")
    print(f"  Exact Lambda mismatches inside k, n-k pairs: {report['exact_lambda_pair_mismatches']}")
    print(f"  Area-key mismatches inside k, n-k pairs: {report['area_pair_mismatches']}")
    print(f"  Distinct area keys: {report['distinct_area_keys']}")
    print(f"  Area collision buckets: {len(report['collision_buckets'])}")
    print(f"  Original scalar found by lookup: {report['lookup_successes']}/{report['table_key_count']}")
    print("  Candidate count returned by an area lookup:")
    for candidate_count in sorted(report["lookup_candidate_counts"]):
        count = report["lookup_candidate_counts"][candidate_count]
        print(f"    {candidate_count} candidates: {count}")
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
