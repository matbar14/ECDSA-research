#!/usr/bin/env python3
"""Screen the normalized relation between areas of P and 2P.

For P = kG this script computes the Lambda-triangle areas S(P) and S(2P)
using the phase model from lambda_circle_triangle_experiment.py.  It studies

    rho(k) = S(2P) / (S(P) + S(2P)).

For positive nonzero areas, rho is in (0, 1) and retains the ordinary ratio:

    S(2P) / S(P) = rho / (1-rho).

The first half of the samples learns the most frequent k mod 2^b in each rho
bin.  The independent second half tests whether this area relation predicts
low scalar bits better than the random baseline.  This is a screening
experiment, not a proof of absence of every possible relation.
"""

from collections import Counter, defaultdict
import math
import random
import time

from cryptography.hazmat.primitives.asymmetric import ec

from lambda_circle_triangle_experiment import (
    GROUP_ORDER,
    lambda_orbit,
    triangle_area,
    unit_circle_vertices,
)
from w1 import DomeInvariantAnalyzer, EllipticCurvePoint


SAMPLE_COUNT = 10000
SEED = 20260925
BIN_COUNT = 256
MODULI = (2, 4, 8, 16, 32, 64)
PROGRESS_INTERVAL = 2000


def generate_256_bit_scalar(generator):
    return generator.randrange(1 << 255, GROUP_ORDER)


def point_from_scalar(scalar):
    public_numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()
    ).public_key().public_numbers()
    return EllipticCurvePoint(public_numbers.x, public_numbers.y)


def area_from_point(analyzer, point):
    values = lambda_orbit(analyzer, point)
    if values is None:
        return None
    return triangle_area(unit_circle_vertices(values))


def normalized_area_ratio(first_area, second_area):
    """Return S(2P)/(S(P)+S(2P)), or None for a degenerate triangle."""
    denominator = first_area + second_area
    if first_area <= 0.0 or second_area <= 0.0 or denominator <= 0.0:
        return None
    return second_area / denominator


def ratio_bin(value):
    return min(BIN_COUNT - 1, int(value * BIN_COUNT))


def most_common_value(counter):
    return max(counter.items(), key=lambda item: (item[1], -item[0]))[0]


def run_experiment(sample_count=SAMPLE_COUNT, seed=SEED):
    if not isinstance(sample_count, int) or sample_count < 2 or sample_count % 2:
        raise ValueError("sample_count must be an even integer greater than one")

    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    training_count = sample_count // 2
    bin_target_counts = {modulus: defaultdict(Counter) for modulus in MODULI}
    global_target_counts = {modulus: Counter() for modulus in MODULI}
    undefined_count = 0
    ratios = []
    test_cases = []
    start_time = time.perf_counter()

    for sample_index in range(sample_count):
        scalar = generate_256_bit_scalar(generator)
        first_point = point_from_scalar(scalar)
        second_point = point_from_scalar((2 * scalar) % GROUP_ORDER)
        first_area = area_from_point(analyzer, first_point)
        second_area = area_from_point(analyzer, second_point)
        ratio = (
            normalized_area_ratio(first_area, second_area)
            if first_area is not None and second_area is not None
            else None
        )
        if ratio is None:
            undefined_count += 1
            continue

        ratios.append(ratio)
        feature = ratio_bin(ratio)
        if sample_index < training_count:
            for modulus in MODULI:
                target = scalar % modulus
                bin_target_counts[modulus][feature][target] += 1
                global_target_counts[modulus][target] += 1
        else:
            test_cases.append((scalar, feature))

        if (sample_index + 1) % PROGRESS_INTERVAL == 0:
            print(f"  Completed: {sample_index + 1}/{sample_count}")

    predictors = {
        modulus: {
            feature: most_common_value(target_counts)
            for feature, target_counts in bin_target_counts[modulus].items()
        }
        for modulus in MODULI
    }
    fallback_predictions = {
        modulus: most_common_value(global_target_counts[modulus])
        for modulus in MODULI
    }
    metrics = {
        modulus: {"correct": 0, "tested": 0, "seen_bin": 0}
        for modulus in MODULI
    }

    for scalar, feature in test_cases:
        for modulus in MODULI:
            prediction = predictors[modulus].get(feature, fallback_predictions[modulus])
            metrics[modulus]["tested"] += 1
            metrics[modulus]["correct"] += prediction == scalar % modulus
            metrics[modulus]["seen_bin"] += feature in predictors[modulus]

    results = []
    for modulus in MODULI:
        metric = metrics[modulus]
        results.append(
            {
                "modulus": modulus,
                "bits": modulus.bit_length() - 1,
                "accuracy": metric["correct"] / metric["tested"],
                "random_baseline": 1 / modulus,
                "seen_bin_coverage": metric["seen_bin"] / metric["tested"],
            }
        )

    return {
        "sample_count": sample_count,
        "seed": seed,
        "undefined_count": undefined_count,
        "ratio_min": min(ratios),
        "ratio_max": max(ratios),
        "ratio_mean": math.fsum(ratios) / len(ratios),
        "results": results,
        "elapsed_seconds": time.perf_counter() - start_time,
    }


def print_report(report):
    print("Area-ratio screening experiment")
    print(f"  256-bit scalars: {report['sample_count']}")
    print(f"  Reproducibility seed: {report['seed']}")
    print(f"  Undefined ratios: {report['undefined_count']}")
    print(f"  rho range: [{report['ratio_min']:.12f}, {report['ratio_max']:.12f}]")
    print(f"  rho mean: {report['ratio_mean']:.12f}")
    print("  modulus  bits  holdout accuracy  random baseline  seen-bin coverage")
    for result in report["results"]:
        print(
            f"  {result['modulus']:7d}  {result['bits']:4d}"
            f"  {result['accuracy']:16.4%}"
            f"  {result['random_baseline']:15.4%}"
            f"  {result['seen_bin_coverage']:17.4%}"
        )
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
