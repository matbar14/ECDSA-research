#!/usr/bin/env python3
"""Test whether Lambda-orbit residues predict low residues of a scalar k.

For each 256-bit scalar k, the experiment computes P = kG and the ordered
triple (Lambda(P), Lambda(phi(P)), Lambda(phi^2(P))). It then tests whether
the triple modulo 2^b predicts k modulo 2^b on a disjoint holdout sample.

This is a statistical screening experiment. A held-out accuracy near the
random baseline means that these residue features do not provide an observed
range reduction for Pollard kangaroo.
"""

from collections import Counter, defaultdict
import random
import time

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import (
    DomeInvariantAnalyzer,
    ENDOMORPHISM_BETA,
    EllipticCurvePoint,
    FIELD_PRIME,
    GROUP_ORDER,
)


SAMPLE_COUNT = 100000
SEED = 20260916
MODULI = (2, 4, 8, 16, 32, 64)
PROGRESS_INTERVAL = 10000


def generate_256_bit_scalar(generator):
    minimum_scalar = 1 << 255
    return generator.randrange(minimum_scalar, GROUP_ORDER)


def point_from_scalar(scalar):
    public_numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(public_numbers.x, public_numbers.y)


def lambda_orbit(analyzer, point):
    phi_point = EllipticCurvePoint(
        ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y
    )
    phi2_point = EllipticCurvePoint(
        ENDOMORPHISM_BETA * phi_point.x % FIELD_PRIME, point.y
    )
    values = tuple(
        analyzer.compute_torus_invariant_lambda(orbit_point)
        for orbit_point in (point, phi_point, phi2_point)
    )
    return values


def most_common_value(counter):
    return max(counter.items(), key=lambda item: (item[1], -item[0]))[0]


def run_experiment(sample_count=SAMPLE_COUNT, seed=SEED):
    if not isinstance(sample_count, int) or sample_count < 2 or sample_count % 2:
        raise ValueError("sample_count must be an even integer greater than one")

    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    training_count = sample_count // 2
    feature_counts = {modulus: defaultdict(Counter) for modulus in MODULI}
    target_counts = {modulus: Counter() for modulus in MODULI}
    undefined_lambda_count = 0
    start_time = time.perf_counter()

    for sample_index in range(training_count):
        scalar = generate_256_bit_scalar(generator)
        values = lambda_orbit(analyzer, point_from_scalar(scalar))
        if any(value is None for value in values):
            undefined_lambda_count += 1
            continue
        for modulus in MODULI:
            target = scalar % modulus
            feature = tuple(value % modulus for value in values)
            feature_counts[modulus][feature][target] += 1
            target_counts[modulus][target] += 1

    predictors = {
        modulus: {
            feature: most_common_value(targets)
            for feature, targets in feature_counts[modulus].items()
        }
        for modulus in MODULI
    }
    fallback_predictions = {
        modulus: most_common_value(target_counts[modulus]) for modulus in MODULI
    }
    metrics = {
        modulus: {"correct": 0, "seen_feature": 0, "test_count": 0}
        for modulus in MODULI
    }

    for sample_index in range(training_count, sample_count):
        scalar = generate_256_bit_scalar(generator)
        values = lambda_orbit(analyzer, point_from_scalar(scalar))
        if any(value is None for value in values):
            undefined_lambda_count += 1
            continue
        for modulus in MODULI:
            target = scalar % modulus
            feature = tuple(value % modulus for value in values)
            prediction = predictors[modulus].get(feature, fallback_predictions[modulus])
            metrics[modulus]["test_count"] += 1
            metrics[modulus]["correct"] += prediction == target
            metrics[modulus]["seen_feature"] += feature in predictors[modulus]

        completed = sample_index + 1
        if completed % PROGRESS_INTERVAL == 0:
            print(f"  Completed: {completed}/{sample_count}")

    results = []
    for modulus in MODULI:
        metric = metrics[modulus]
        test_count = metric["test_count"]
        results.append(
            {
                "modulus": modulus,
                "bits": modulus.bit_length() - 1,
                "accuracy": metric["correct"] / test_count,
                "random_baseline": 1 / modulus,
                "feature_coverage": metric["seen_feature"] / test_count,
            }
        )

    return {
        "sample_count": sample_count,
        "seed": seed,
        "undefined_lambda_count": undefined_lambda_count,
        "results": results,
        "elapsed_seconds": time.perf_counter() - start_time,
    }


def print_report(report):
    print("Invariant leakage experiment")
    print(f"  256-bit scalars: {report['sample_count']}")
    print(f"  Reproducibility seed: {report['seed']}")
    print(f"  Undefined Lambda values: {report['undefined_lambda_count']}")
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
