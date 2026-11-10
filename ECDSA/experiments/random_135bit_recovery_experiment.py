#!/usr/bin/env python3
"""Probe quantized-angle recovery for one random 135-bit scalar.

The full 135-bit interval is not enumerated. Instead, a random sparse table is
sampled across [2^134, 2^135), with the target scalar inserted as a control.
For each quantization size the script reports the bucket span, nearest
candidate, and exact (unquantized) matches.
"""

import argparse
import random
from collections import defaultdict

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import (
    DomeInvariantAnalyzer,
    EllipticCurvePoint,
    ENDOMORPHISM_BETA,
    FIELD_PRIME,
)


def point_from_scalar(k):
    q = ec.derive_private_key(k, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(q.x, q.y)


def angle_pair(analyzer, point):
    first = analyzer.compute_torus_invariant_lambda(point)
    second_point = EllipticCurvePoint(
        ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y
    )
    second = analyzer.compute_torus_invariant_lambda(second_point)
    return first, second


def quantized_pair(values, bits):
    return tuple(value * (1 << bits) // FIELD_PRIME for value in values)


def compressed_key(point):
    prefix = 0x02 | (point.y & 1)
    return f"{prefix:02x}{point.x:064x}"


def run(samples, bits, seed, include_target=True, key_mode="pair"):
    rng = random.Random(seed)
    start = 1 << 134
    width = 1 << 134
    target = rng.randrange(start, start + width)
    analyzer = DomeInvariantAnalyzer()
    target_point = point_from_scalar(target)
    target_exact = angle_pair(analyzer, target_point)
    rows = [target] if include_target else []
    rows.extend(rng.randrange(start, start + width) for _ in range(samples - len(rows)))
    report = []

    exact_cache = {target: target_exact}
    for scalar in rows:
        if scalar not in exact_cache:
            exact_cache[scalar] = angle_pair(analyzer, point_from_scalar(scalar))
    for bit_count in bits:
        buckets = defaultdict(list)
        table_values = [(scalar, exact_cache[scalar]) for scalar in rows]
        for scalar, values in table_values:
            key = (quantized_pair(values, bit_count) if key_mode == "pair"
                   else quantized_pair((values[0],), bit_count))
            buckets[key].append(scalar)
        target_key = (quantized_pair(target_exact, bit_count) if key_mode == "pair"
                      else quantized_pair((target_exact[0],), bit_count))
        bucket = buckets[target_key]
        exact_matches = [
            scalar for scalar in bucket if exact_cache[scalar] == target_exact
        ]
        distances = sorted(abs(scalar - target) for scalar in bucket)
        report.append({
            "bits": bit_count,
            "bucket_size": len(bucket),
            "bucket_span": max(bucket) - min(bucket) + 1 if bucket else 0,
            "bucket_span_fraction": ((max(bucket) - min(bucket) + 1) / width
                                      if bucket else 0.0),
            "nearest_distance": distances[0] if distances else None,
            "exact_matches": exact_matches,
        })
    return target, start, width, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=2000)
    parser.add_argument("--bits", type=int, nargs="+", default=[8, 12, 16, 20])
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--exclude-target", action="store_true",
                        help="do not insert target into the sparse table")
    parser.add_argument("--key", choices=("one", "pair"), default="pair",
                        help="quantized key: one angle or two-angle pair")
    args = parser.parse_args()
    target, start, width, report = run(
        args.samples, args.bits, args.seed, not args.exclude_target, args.key
    )
    print("Random 135-bit quantized-angle recovery")
    print(f"  target k: 0x{target:034x}")
    print(f"  target bit length: {target.bit_length()}")
    print(f"  target public key (compressed): {compressed_key(point_from_scalar(target))}")
    print(f"  sparse table rows: {args.samples}")
    print("  bits  bucket  span fraction  nearest distance  exact matches")
    for result in report:
        print(
            f"  {result['bits']:4d}  {result['bucket_size']:6d}"
            f"  {result['bucket_span_fraction']:13.6%}"
            f"  {str(result['nearest_distance']):>16s}"
            f"  {len(result['exact_matches'])}"
        )
    print(f"  target inserted in table: {not args.exclude_target}")
    print(f"  quantized key: {args.key}")
    print("  This measures bucket locality, not full 2^135 coverage.")


if __name__ == "__main__":
    main()
