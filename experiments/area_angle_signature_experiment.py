#!/usr/bin/env python3
"""Collision experiment for the first-radius angle and triangle area.

For P = kG, the first radius is defined by the exact Lambda value from w1.py:

    theta(P) = Arg(exp(2*pi*i*Lambda(P)/p)).

The triangle area S(P) is built from the three ordered Lambda values at
P, phi(P), and phi^2(P).  The geometric signature is therefore

    (theta(P), S(P)).

For comparison without loss from floating-point angle representation, theta
is stored as the exact residue Lambda(P); the displayed angle is only its
unit-circle realization.  The experiment checks whether S distinguishes the
x-coordinate roots that remain when theta alone is fixed.
"""

from collections import Counter
import importlib
import math
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
SEED = 20261012
PROGRESS_INTERVAL = 1000


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


def first_radius_angle(lambda_value):
    """Return theta in radians for display; the residue remains the exact key."""
    return math.tau * (lambda_value / FIELD_PRIME)


def triangle_area(values):
    angles = tuple(first_radius_angle(value) for value in values)
    vertices = tuple(complex(math.cos(angle), math.sin(angle)) for angle in angles)
    first, second, third = vertices
    return abs(
        first.real * (second.imag - third.imag)
        + second.real * (third.imag - first.imag)
        + third.real * (first.imag - second.imag)
    ) / 2.0


def area_angle_signature(values):
    """Return exact first-radius residue and stable binary64 area encoding."""
    area = triangle_area(values)
    return values[0], area.hex(), area


def x_coordinates(points):
    return {point.x for point in points}


def run_experiment(sample_count=SAMPLE_COUNT, seed=SEED):
    if not isinstance(sample_count, int) or sample_count <= 0:
        raise ValueError("sample_count must be a positive integer")

    inversion = importlib.import_module("2edcl")
    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    counters = Counter()
    source_signatures = set()
    start_time = time.perf_counter()

    for sample_index in range(1, sample_count + 1):
        scalar = generate_256_bit_scalar(generator)
        point = point_from_scalar(scalar)
        values = lambda_orbit(analyzer, point)
        if values is None:
            counters["undefined_lambda"] += 1
            continue

        angle_residue, area_key, area = area_angle_signature(values)
        signature = angle_residue, area_key
        counters["repeated_source_signature"] += signature in source_signatures
        source_signatures.add(signature)

        lambda_candidates = inversion.reconstruct_points_from_lambda(angle_residue)
        lambda_xs = x_coordinates(lambda_candidates)
        signature_candidates = []
        candidate_values_by_x = {}
        for candidate in lambda_candidates:
            candidate_point = EllipticCurvePoint(candidate.x, candidate.y)
            candidate_values = lambda_orbit(analyzer, candidate_point)
            if candidate_values is None:
                continue
            candidate_values_by_x[candidate.x] = candidate_values
            candidate_angle, candidate_area_key, _ = area_angle_signature(candidate_values)
            if (candidate_angle, candidate_area_key) == signature:
                signature_candidates.append(candidate)
        signature_xs = x_coordinates(signature_candidates)

        counters["processed"] += 1
        counters[f"theta_x_candidates_{len(lambda_xs)}"] += 1
        counters[f"area_angle_x_candidates_{len(signature_xs)}"] += 1
        counters["original_x_retained"] += point.x in signature_xs
        counters["area_near_zero"] += area < 1e-12

        for candidate_x, candidate_values in candidate_values_by_x.items():
            if candidate_x == point.x:
                continue
            candidate_angle, candidate_area_key, _ = area_angle_signature(candidate_values)
            if (candidate_angle, candidate_area_key) == signature:
                counters["secondary_x_area_match"] += 1
                if sorted(candidate_values) == sorted(values):
                    counters["secondary_x_same_vertex_set"] += 1
                else:
                    counters["secondary_x_area_match_different_vertex_set"] += 1

        if sample_index % PROGRESS_INTERVAL == 0:
            print(f"  Completed: {sample_index}/{sample_count}")

    counters["elapsed_seconds"] = time.perf_counter() - start_time
    return counters


def print_report(counters):
    processed = counters["processed"]
    print("First-radius angle and triangle-area experiment")
    print(f"  Processed 256-bit scalars: {processed}")
    print(f"  Undefined Lambda triples: {counters['undefined_lambda']}")
    print(f"  Repeated source signatures: {counters['repeated_source_signature']}")
    print(f"  Original x retained by (theta, area): {counters['original_x_retained']}/{processed}")
    print(f"  Areas below 1e-12: {counters['area_near_zero']}")
    print(f"  Second x values with the same (theta, area): {counters['secondary_x_area_match']}")
    print(f"    Same unordered Lambda vertices: {counters['secondary_x_same_vertex_set']}")
    print(f"    Different vertices with same binary64 area: {counters['secondary_x_area_match_different_vertex_set']}")
    print("  x-coordinate candidates from theta alone:")
    for count in range(3):
        if counters[f"theta_x_candidates_{count}"]:
            print(f"    {count} x values: {counters[f'theta_x_candidates_{count}']}")
    print("  x-coordinate candidates from (theta, area):")
    for count in range(3):
        if counters[f"area_angle_x_candidates_{count}"]:
            print(f"    {count} x values: {counters[f'area_angle_x_candidates_{count}']}")
    print(f"  Elapsed time: {counters['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
