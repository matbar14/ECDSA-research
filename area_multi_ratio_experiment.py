#!/usr/bin/env python3
"""Screen joint area relations for P, 2P, 3P, and 5P.

For P = kG, define S(Q) as the Lambda-triangle area of Q under the explicit
phase model.  This experiment uses the three bounded area relations

    rho_m(k) = S(mP) / (S(P) + S(mP)),  m in {2, 3, 5}.

Each rho_m contains the same information as S(mP)/S(P) when both areas are
nonzero.  The tuple is quantized into a small joint grid, learned on one half
of the sample, and evaluated on a disjoint holdout half for k mod 2^b.
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


SAMPLE_COUNT = 6000
SEED = 20260926
MULTIPLIERS = (2, 3, 5)
BIN_COUNT = 8
MODULI = (2, 4, 8, 16, 32, 64)
PROGRESS_INTERVAL = 1000


def generate_256_bit_scalar(generator):
    return generator.randrange(1 << 255, GROUP_ORDER)


def point_from_scalar(scalar):
    public_numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()
    ).public_key().public_numbers()
    return EllipticCurvePoint(public_numbers.x, public_numbers.y)


def area_from_scalar(analyzer, scalar):
    values = lambda_orbit(analyzer, point_from_scalar(scalar))
    if values is None:
        return None
    return triangle_area(unit_circle_vertices(values))


def normalized_ratio(base_area, multiple_area):
    denominator = base_area + multiple_area
    if base_area <= 0.0 or multiple_area <= 0.0 or denominator <= 0.0:
        return None
    return multiple_area / denominator


def feature_from_ratios(ratios):
    return tuple(min(BIN_COUNT - 1, int(value * BIN_COUNT)) for value in ratios)


def most_common_value(counter):
    return max(counter.items(), key=lambda item: (item[1], -item[0]))[0]


def run_experiment(sample_count=SAMPLE_COUNT, seed=SEED):
    if not isinstance(sample_count, int) or sample_count < 2 or sample_count % 2:
        raise ValueError("sample_count must be an even integer greater than one")

    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    training_count = sample_count // 2
    feature_targets = {modulus: defaultdict(Counter) for modulus in MODULI}
    global_targets = {modulus: Counter() for modulus in MODULI}
    test_cases = []
    ratio_sums = [0.0 for _ in MULTIPLIERS]
    undefined_count = 0
    start_time = time.perf_counter()

    for sample_index in range(sample_count):
        scalar = generate_256_bit_scalar(generator)
        base_area = area_from_scalar(analyzer, scalar)
        multiple_areas = [
            area_from_scalar(analyzer, scalar * multiplier % GROUP_ORDER)
            for multiplier in MULTIPLIERS
        ]
        if base_area is None or any(area is None for area in multiple_areas):
            undefined_count += 1
            continue

        ratios = tuple(
            normalized_ratio(base_area, multiple_area)
            for multiple_area in multiple_areas
        )
        if any(ratio is None for ratio in ratios):
            undefined_count += 1
            continue

        for index, ratio in enumerate(ratios):
            ratio_sums[index] += ratio
        feature = feature_from_ratios(ratios)
        if sample_index < training_count:
            for modulus in MODULI:
                target = scalar % modulus
                feature_targets[modulus][feature][target] += 1
                global_targets[modulus][target] += 1
        else:
            test_cases.append((scalar, feature))

        if (sample_index + 1) % PROGRESS_INTERVAL == 0:
            print(f"  Completed: {sample_index + 1}/{sample_count}")

    predictors = {
        modulus: {
            feature: most_common_value(target_counts)
            for feature, target_counts in feature_targets[modulus].items()
        }
        for modulus in MODULI
    }
    fallback_predictions = {
        modulus: most_common_value(global_targets[modulus]) for modulus in MODULI
    }
    metrics = {
        modulus: {"correct": 0, "tested": 0, "seen_feature": 0}
        for modulus in MODULI
    }
    for scalar, feature in test_cases:
        for modulus in MODULI:
            prediction = predictors[modulus].get(feature, fallback_predictions[modulus])
            metrics[modulus]["tested"] += 1
            metrics[modulus]["correct"] += prediction == scalar % modulus
            metrics[modulus]["seen_feature"] += feature in predictors[modulus]

    results = []
    for modulus in MODULI:
        metric = metrics[modulus]
        results.append(
            {
                "modulus": modulus,
                "bits": modulus.bit_length() - 1,
                "accuracy": metric["correct"] / metric["tested"],
                "random_baseline": 1 / modulus,
                "feature_coverage": metric["seen_feature"] / metric["tested"],
            }
        )

    processed = sample_count - undefined_count
    return {
        "sample_count": sample_count,
        "seed": seed,
        "undefined_count": undefined_count,
        "ratio_means": tuple(value / processed for value in ratio_sums),
        "training_feature_count": len(feature_targets[MODULI[0]]),
        "results": results,
        "elapsed_seconds": time.perf_counter() - start_time,
    }


def print_report(report):
    print("Joint area-ratio screening experiment")
    print(f"  256-bit scalars: {report['sample_count']}")
    print(f"  Reproducibility seed: {report['seed']}")
    print(f"  Multipliers: {MULTIPLIERS}")
    print(f"  Joint bins: {BIN_COUNT}^{len(MULTIPLIERS)}")
    print(f"  Undefined ratio vectors: {report['undefined_count']}")
    print(f"  Training features observed: {report['training_feature_count']}")
    print("  Mean rho values:")
    for multiplier, value in zip(MULTIPLIERS, report["ratio_means"]):
        print(f"    rho_{multiplier}: {value:.12f}")
    print("  modulus  bits  holdout accuracy  random baseline  feature coverage")
    for result in report["results"]:
        print(
            f"  {result['modulus']:7d}  {result['bits']:4d}"
            f"  {result['accuracy']:16.4%}"
            f"  {result['random_baseline']:15.4%}"
            f"  {result['feature_coverage']:16.4%}"
        )
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
