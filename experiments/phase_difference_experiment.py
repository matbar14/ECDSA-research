#!/usr/bin/env python3
"""Exact phase-difference experiment for the ordered Lambda orbit.

For the ordered Lambda triple (L0, L1, L2), define exact phase differences
over F_p:

    d01 = L1-L0,  d12 = L2-L1,  d20 = L0-L2.

Because L0+L1+L2 = 0 and p != 3, the pair (d01, d12) reconstructs the whole
ordered triple:

    L0 = -(2*d01+d12)/3,
    L1 = L0+d01,
    L2 = L1+d12.

The script compares candidate points obtained from L0 alone with candidates
that also match the complete ordered phase-difference signature.  It uses the
exact Lambda inversion already implemented in 2edcl.py; no decimal angle or
floating-point pi is used in the comparison.
"""

from collections import Counter
import importlib
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


SAMPLE_COUNT = 10000
SEED = 20261010
PROGRESS_INTERVAL = 1000
INVERSE_THREE = pow(3, -1, FIELD_PRIME)


def generate_256_bit_scalar(generator):
    return generator.randrange(1 << 255, GROUP_ORDER)


def point_from_scalar(scalar):
    public_numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()
    ).public_key().public_numbers()
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
    if any(value is None for value in values):
        return None
    return values


def phase_signature(values):
    first, second, third = values
    return (
        (second - first) % FIELD_PRIME,
        (third - second) % FIELD_PRIME,
    )


def reconstruct_lambda_triple(signature):
    """Recover ordered Lambda values from d01 and d12 exactly in F_p."""
    difference_01, difference_12 = signature
    first = -(2 * difference_01 + difference_12) * INVERSE_THREE % FIELD_PRIME
    second = (first + difference_01) % FIELD_PRIME
    third = (second + difference_12) % FIELD_PRIME
    return first, second, third


def w1_point(candidate):
    return EllipticCurvePoint(candidate.x, candidate.y)


def run_experiment(sample_count=SAMPLE_COUNT, seed=SEED):
    if not isinstance(sample_count, int) or sample_count <= 0:
        raise ValueError("sample_count must be a positive integer")

    inversion = importlib.import_module("2edcl")
    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    counters = Counter()
    start_time = time.perf_counter()

    for sample_index in range(1, sample_count + 1):
        scalar = generate_256_bit_scalar(generator)
        point = point_from_scalar(scalar)
        values = lambda_orbit(analyzer, point)
        if values is None:
            counters["undefined_lambda"] += 1
            continue

        signature = phase_signature(values)
        reconstructed_values = reconstruct_lambda_triple(signature)
        if reconstructed_values != values:
            counters["signature_reconstruction_failure"] += 1
            continue

        lambda_candidates = inversion.reconstruct_points_from_lambda(values[0])
        signature_candidates = []
        for candidate in lambda_candidates:
            candidate_values = lambda_orbit(analyzer, w1_point(candidate))
            if candidate_values is not None and phase_signature(candidate_values) == signature:
                signature_candidates.append(candidate)

        counters["processed"] += 1
        counters[f"lambda_candidates_{len(lambda_candidates)}"] += 1
        counters[f"signature_candidates_{len(signature_candidates)}"] += 1
        counters["original_in_signature_candidates"] += any(
            candidate.x == point.x and candidate.y == point.y
            for candidate in signature_candidates
        )
        counters["negated_original_in_signature_candidates"] += any(
            candidate.x == point.x and candidate.y == (-point.y) % FIELD_PRIME
            for candidate in signature_candidates
        )

        if sample_index % PROGRESS_INTERVAL == 0:
            print(f"  Completed: {sample_index}/{sample_count}")

    counters["elapsed_seconds"] = time.perf_counter() - start_time
    return counters


def print_report(counters):
    processed = counters["processed"]
    print("Exact phase-difference experiment")
    print(f"  Processed 256-bit scalars: {processed}")
    print(f"  Undefined Lambda triples: {counters['undefined_lambda']}")
    print(f"  Phase-signature reconstruction failures: {counters['signature_reconstruction_failure']}")
    print(f"  Original P in signature candidates: {counters['original_in_signature_candidates']}/{processed}")
    print(f"  Negated P in signature candidates: {counters['negated_original_in_signature_candidates']}/{processed}")
    print("  Candidates from Lambda(P) alone:")
    for count in range(5):
        if counters[f"lambda_candidates_{count}"]:
            print(f"    {count} points: {counters[f'lambda_candidates_{count}']}")
    print("  Candidates after ordered phase signature (d01, d12):")
    for count in range(5):
        if counters[f"signature_candidates_{count}"]:
            print(f"    {count} points: {counters[f'signature_candidates_{count}']}")
    print(f"  Elapsed time: {counters['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
