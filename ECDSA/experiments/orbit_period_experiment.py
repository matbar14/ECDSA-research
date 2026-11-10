#!/usr/bin/env python3
"""Compare candidate turn-count functions with the endomorphism orbit period.

Observed value:
    d_observed(P) = min {m > 0 : phi^m(P) = P}

Candidate functions:
    d_constant(k) = 3
    d_lambda(k) = 1 + (Lambda(kG) mod TURN_MODULUS)
"""

import secrets
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
TURN_MODULUS = 256


def generate_256_bit_scalar():
    """Return a canonical secp256k1 scalar of exactly 256 bits."""
    minimum = 1 << 255
    return minimum + secrets.randbelow(GROUP_ORDER - minimum)


def endomorphism(point):
    return EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)


def observed_orbit_period(analyzer, projection):
    """Find the first positive m <= 3 for which phi^m(P) = P."""
    orbit_point = projection
    for exponent in range(1, 4):
        orbit_point = endomorphism(orbit_point)
        if analyzer.points_equal(orbit_point, projection):
            return exponent
    return None


def lambda_turn_count(analyzer, projection, turn_modulus):
    lambda_value = analyzer.compute_torus_invariant_lambda(projection)
    if lambda_value is None:
        return None
    return 1 + lambda_value % turn_modulus


def run_experiment(sample_count=SAMPLE_COUNT, turn_modulus=TURN_MODULUS):
    if not isinstance(sample_count, int) or sample_count <= 0:
        raise ValueError("sample_count must be a positive integer")
    if not isinstance(turn_modulus, int) or turn_modulus <= 0:
        raise ValueError("turn_modulus must be a positive integer")

    analyzer = DomeInvariantAnalyzer()
    counts = {
        "invalid_curve_point": 0,
        "undefined_period": 0,
        "period_not_three": 0,
        "constant_three_matches": 0,
        "lambda_undefined": 0,
        "lambda_matches": 0,
    }
    start_time = time.perf_counter()

    for sample_index in range(1, sample_count + 1):
        scalar = generate_256_bit_scalar()
        public_numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
        projection = EllipticCurvePoint(public_numbers.x, public_numbers.y)

        if not analyzer.is_valid_curve_point(projection):
            counts["invalid_curve_point"] += 1
            continue

        observed_turns = observed_orbit_period(analyzer, projection)
        if observed_turns is None:
            counts["undefined_period"] += 1
            continue
        if observed_turns != 3:
            counts["period_not_three"] += 1
        if observed_turns == 3:
            counts["constant_three_matches"] += 1

        predicted_turns = lambda_turn_count(analyzer, projection, turn_modulus)
        if predicted_turns is None:
            counts["lambda_undefined"] += 1
        elif predicted_turns == observed_turns:
            counts["lambda_matches"] += 1

        if sample_index % 10000 == 0:
            print(f"  Completed: {sample_index}/{sample_count}")

    elapsed_seconds = time.perf_counter() - start_time
    return {
        "sample_count": sample_count,
        "turn_modulus": turn_modulus,
        "counts": counts,
        "elapsed_seconds": elapsed_seconds,
    }


def print_report(result):
    counts = result["counts"]
    defined_periods = result["sample_count"] - counts["invalid_curve_point"] - counts["undefined_period"]
    defined_lambdas = defined_periods - counts["lambda_undefined"]

    print("Endomorphism orbit-period experiment")
    print(f"  256-bit scalars: {result['sample_count']}")
    print(f"  TURN_MODULUS: {result['turn_modulus']}")
    print(f"  Invalid curve points: {counts['invalid_curve_point']}")
    print(f"  Undefined orbit periods: {counts['undefined_period']}")
    print(f"  Observed periods not equal to 3: {counts['period_not_three']}")
    print(f"  d_constant(k) = 3 matches: {counts['constant_three_matches']}/{defined_periods}")
    print(f"  Undefined Lambda values: {counts['lambda_undefined']}")
    print(f"  d_lambda(k) matches: {counts['lambda_matches']}/{defined_lambdas}")
    print(f"  Elapsed time: {result['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
