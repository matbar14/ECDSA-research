#!/usr/bin/env python3
"""Experiment: derive a torus turn count d from a 256-bit height k.

Candidate relation:

    d(k) = 1 + (Lambda(kG) mod TURN_MODULUS)

The mapping is deterministic wherever Lambda(kG) is defined. Consequently,
the same height k always reconstructs the same turn count d.
"""

from collections import Counter
from fractions import Fraction
import secrets

from w1 import DomeInvariantAnalyzer, GROUP_ORDER


SAMPLE_COUNT = 1000
TURN_MODULUS = 256


def generate_256_bit_scalar():
    """Return a canonical secp256k1 scalar with bit length exactly 256."""
    minimum = 1 << 255
    return minimum + secrets.randbelow(GROUP_ORDER - minimum)


def derive_turn_count(analyzer, height, turn_modulus):
    """Evaluate d(k) and return d together with P and Lambda(P)."""
    projection = height * analyzer.base_point
    lambda_value = analyzer.compute_torus_invariant_lambda(projection)
    if lambda_value is None:
        return None, projection, None
    return 1 + lambda_value % turn_modulus, projection, lambda_value


def run_experiment(sample_count=SAMPLE_COUNT, turn_modulus=TURN_MODULUS):
    if not isinstance(sample_count, int) or sample_count <= 0:
        raise ValueError("sample_count must be a positive integer")
    if not isinstance(turn_modulus, int) or turn_modulus <= 0:
        raise ValueError("turn_modulus must be a positive integer")

    analyzer = DomeInvariantAnalyzer()
    undefined_lambda_count = 0
    recovery_failures = 0
    phase_failures = 0
    turn_counts = Counter()

    for _ in range(sample_count):
        height = generate_256_bit_scalar()
        turn_count, _, lambda_value = derive_turn_count(analyzer, height, turn_modulus)

        if lambda_value is None:
            undefined_lambda_count += 1
            continue

        recovered_turn_count, _, _ = derive_turn_count(analyzer, height, turn_modulus)
        if recovered_turn_count != turn_count:
            recovery_failures += 1

        if Fraction(turn_count * height, height) != turn_count:
            phase_failures += 1

        turn_counts[turn_count] += 1

    return {
        "sample_count": sample_count,
        "turn_modulus": turn_modulus,
        "undefined_lambda_count": undefined_lambda_count,
        "recovery_failures": recovery_failures,
        "phase_failures": phase_failures,
        "distinct_turn_counts": len(turn_counts),
        "minimum_turn_count": min(turn_counts) if turn_counts else None,
        "maximum_turn_count": max(turn_counts) if turn_counts else None,
    }


def print_report(result):
    print("d = f(k) experiment")
    print("  f(k) = 1 + (Lambda(kG) mod TURN_MODULUS)")
    print(f"  256-bit heights: {result['sample_count']}")
    print(f"  TURN_MODULUS: {result['turn_modulus']}")
    print(f"  Undefined Lambda values: {result['undefined_lambda_count']}")
    print(f"  d recovery failures: {result['recovery_failures']}")
    print(f"  Endpoint phase failures: {result['phase_failures']}")
    print(f"  Distinct d values: {result['distinct_turn_counts']}")
    print(
        "  d range: "
        f"{result['minimum_turn_count']}..{result['maximum_turn_count']}"
    )


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
