#!/usr/bin/env python3
"""Exact symmetric invariants of the three Lambda values on secp256k1.

Let z = x(P)/x(G), u = z^3, and let L0, L1, L2 be Lambda evaluated at
P, phi(P), phi^2(P), respectively.  Reducing with beta^2+beta+1 = 0 gives

    I2 = L0*L1 + L1*L2 + L2*L0 = J(u),
    I3 = L0*L1*L2 = J(u),
    J(u) = 27*u/(u-1)^2,
    D = product_{i<j}(Li-Lj)^2 = -J(u)^2*(4*J(u)+27).

Thus I2, I3, and D are algebraically dependent.  The inverse equation is

    J*u^2 - (2*J+27)*u + J = 0,

whose roots are u and u^-1.  The script verifies these identities on random
256-bit scalars and counts point candidates in the inverse invariant class.
"""

from collections import Counter
import random
import time

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import (
    BASE_POINT,
    DomeInvariantAnalyzer,
    ENDOMORPHISM_BETA,
    EllipticCurvePoint,
    FIELD_PRIME,
    GROUP_ORDER,
)


SAMPLE_COUNT = 100000
SEED = 20260929
PROGRESS_INTERVAL = 10000


def generate_256_bit_scalar(generator):
    return generator.randrange(1 << 255, GROUP_ORDER)


def point_from_scalar(scalar):
    public_numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()
    ).public_key().public_numbers()
    return EllipticCurvePoint(public_numbers.x, public_numbers.y)


def lambda_values(analyzer, point):
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


def symmetric_invariants(values):
    first, second, third = values
    i2 = (first * second + second * third + third * first) % FIELD_PRIME
    i3 = first * second * third % FIELD_PRIME
    discriminant = (
        (first - second) ** 2
        * (second - third) ** 2
        * (third - first) ** 2
    ) % FIELD_PRIME
    return i2, i3, discriminant


def invariant_j_from_u(u):
    denominator = (u - 1) % FIELD_PRIME
    if denominator == 0:
        return None
    return 27 * u * pow(denominator, -2, FIELD_PRIME) % FIELD_PRIME


def square_root_mod_prime(value):
    value %= FIELD_PRIME
    if value == 0:
        return 0
    if pow(value, (FIELD_PRIME - 1) // 2, FIELD_PRIME) != 1:
        return None
    return pow(value, (FIELD_PRIME + 1) // 4, FIELD_PRIME)


def inverse_u_roots(invariant_j):
    """Solve J*u^2 -(2J+27)*u + J = 0 over F_p."""
    if invariant_j == 0:
        return [0]
    coefficient_a = invariant_j
    coefficient_b = -(2 * invariant_j + 27) % FIELD_PRIME
    discriminant = (
        coefficient_b * coefficient_b - 4 * coefficient_a * invariant_j
    ) % FIELD_PRIME
    root = square_root_mod_prime(discriminant)
    if root is None:
        return []
    denominator_inverse = pow(2 * coefficient_a, -1, FIELD_PRIME)
    first = (-coefficient_b + root) * denominator_inverse % FIELD_PRIME
    second = (-coefficient_b - root) * denominator_inverse % FIELD_PRIME
    return [first] if first == second else [first, second]


def x_has_curve_points(x_coordinate):
    right_side = (pow(x_coordinate, 3, FIELD_PRIME) + 7) % FIELD_PRIME
    return pow(right_side, (FIELD_PRIME - 1) // 2, FIELD_PRIME) == 1


def inverse_polynomial_value(invariant_j, u):
    """Evaluate J*u^2 -(2J+27)*u + J modulo p."""
    return (
        invariant_j * u * u - (2 * invariant_j + 27) * u + invariant_j
    ) % FIELD_PRIME


def inverse_class_candidate_count(z):
    """Count affine points in the symmetric invariant class of z.

    The u branch always contributes the three x values z, beta*z, beta^2*z,
    each with two signs of y.  The u^-1 branch contributes another six points
    precisely when x(G)/z is the x-coordinate of a curve point.
    """
    inverse_z = pow(z, -1, FIELD_PRIME)
    inverse_x = inverse_z * BASE_POINT.x % FIELD_PRIME
    return 12 if x_has_curve_points(inverse_x) else 6


def run_experiment(sample_count=SAMPLE_COUNT, seed=SEED):
    if not isinstance(sample_count, int) or sample_count <= 0:
        raise ValueError("sample_count must be a positive integer")

    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    failures = Counter()
    candidate_counts = Counter()
    start_time = time.perf_counter()

    for sample_index in range(1, sample_count + 1):
        scalar = generate_256_bit_scalar(generator)
        point = point_from_scalar(scalar)
        values = lambda_values(analyzer, point)
        if values is None:
            failures["undefined_lambda"] += 1
            continue

        z = point.x * pow(BASE_POINT.x, -1, FIELD_PRIME) % FIELD_PRIME
        u = pow(z, 3, FIELD_PRIME)
        invariant_j = invariant_j_from_u(u)
        i2, i3, discriminant = symmetric_invariants(values)
        expected_discriminant = (
            -invariant_j * invariant_j * (4 * invariant_j + 27)
        ) % FIELD_PRIME
        inverse_u = pow(u, -1, FIELD_PRIME)

        failures["i2_not_i3"] += i2 != i3
        failures["i2_not_j"] += i2 != invariant_j
        failures["discriminant_not_j"] += discriminant != expected_discriminant
        failures["orbit_sum_not_zero"] += sum(values) % FIELD_PRIME != 0
        failures["inverse_u_residual_nonzero"] += (
            inverse_polynomial_value(invariant_j, u) != 0
        )
        failures["inverse_u_inverse_residual_nonzero"] += (
            inverse_polynomial_value(invariant_j, inverse_u) != 0
        )
        failures["inverse_u_roots_collide"] += u == inverse_u
        candidate_counts[inverse_class_candidate_count(z)] += 1

        if sample_index % PROGRESS_INTERVAL == 0:
            print(f"  Completed: {sample_index}/{sample_count}")

    return {
        "sample_count": sample_count,
        "seed": seed,
        "failures": failures,
        "candidate_counts": candidate_counts,
        "elapsed_seconds": time.perf_counter() - start_time,
    }


def print_report(report):
    print("Symmetric Lambda-orbit invariant experiment")
    print(f"  256-bit scalars: {report['sample_count']}")
    print(f"  Reproducibility seed: {report['seed']}")
    print("  Identity failures:")
    for name in (
        "undefined_lambda",
        "orbit_sum_not_zero",
        "i2_not_i3",
        "i2_not_j",
        "discriminant_not_j",
        "inverse_u_residual_nonzero",
        "inverse_u_inverse_residual_nonzero",
        "inverse_u_roots_collide",
    ):
        print(f"    {name}: {report['failures'][name]}")
    print("  Inverse symmetric-invariant point candidates:")
    for candidate_count in sorted(report["candidate_counts"]):
        count = report["candidate_counts"][candidate_count]
        print(f"    {candidate_count} points: {count}")
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
