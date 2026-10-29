#!/usr/bin/env python3
"""Measure what the quadratic Lambda inversion adds to the torus lookup.

One Lambda value yields up to two x-coordinates and two y-signs per x.  This
script checks whether the second x branch is a new GLV-orbit representative,
and whether the full ordered Lambda pair removes it.
"""

import argparse
import random
from collections import Counter

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import (
    BASE_POINT,
    DomeInvariantAnalyzer,
    EllipticCurvePoint,
    ENDOMORPHISM_BETA,
    FIELD_PRIME,
    GROUP_ORDER,
)


def point_from_scalar(k):
    q = ec.derive_private_key(k, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(q.x, q.y)


def lambda_coefficients(value):
    beta = ENDOMORPHISM_BETA
    xg = BASE_POINT.x
    beta2 = beta * beta % FIELD_PRIME
    a = -value * beta % FIELD_PRIME
    b = (value * xg * (1 + beta2) - xg * pow(1 - beta, 2, FIELD_PRIME)) % FIELD_PRIME
    c = -value * beta * xg * xg % FIELD_PRIME
    return a, b, c


def sqrt_mod(value):
    value %= FIELD_PRIME
    if value == 0:
        return 0
    if pow(value, (FIELD_PRIME - 1) // 2, FIELD_PRIME) != 1:
        return None
    root = pow(value, (FIELD_PRIME + 1) // 4, FIELD_PRIME)
    return root if root * root % FIELD_PRIME == value else None


def x_roots(value):
    a, b, c = lambda_coefficients(value)
    if a == 0:
        return [] if b == 0 else [(-c * pow(b, -1, FIELD_PRIME)) % FIELD_PRIME]
    disc = (b * b - 4 * a * c) % FIELD_PRIME
    root = sqrt_mod(disc)
    if root is None:
        return []
    inv = pow(2 * a, -1, FIELD_PRIME)
    first = (-b + root) * inv % FIELD_PRIME
    second = (-b - root) * inv % FIELD_PRIME
    return [first] if first == second else [first, second]


def curve_points_for_x(x):
    rhs = (x * x * x + 7) % FIELD_PRIME
    y = sqrt_mod(rhs)
    if y is None:
        return []
    return [EllipticCurvePoint(x, y), EllipticCurvePoint(x, -y % FIELD_PRIME)]


def orbit_xs(point):
    current = point
    result = set()
    for _ in range(3):
        result.add(current.x)
        current = EllipticCurvePoint(ENDOMORPHISM_BETA * current.x % FIELD_PRIME, current.y)
    return result


def ordered_pair(analyzer, point):
    first = analyzer.compute_torus_invariant_lambda(point)
    second_point = EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)
    second = analyzer.compute_torus_invariant_lambda(second_point)
    return first, second


def run(samples, seed):
    rng = random.Random(seed)
    analyzer = DomeInvariantAnalyzer()
    counts = Counter()
    for _ in range(samples):
        scalar = rng.randrange(1, GROUP_ORDER)
        point = point_from_scalar(scalar)
        first = analyzer.compute_torus_invariant_lambda(point)
        roots = x_roots(first)
        candidates = [candidate for root in roots for candidate in curve_points_for_x(root)]
        valid_lambda = [candidate for candidate in candidates
                        if analyzer.compute_torus_invariant_lambda(candidate) == first]
        target_pair = ordered_pair(analyzer, point)
        pair_candidates = [candidate for candidate in valid_lambda
                           if ordered_pair(analyzer, candidate) == target_pair]
        counts[f"x_roots_{len(roots)}"] += 1
        counts[f"lambda_points_{len(valid_lambda)}"] += 1
        counts[f"ordered_pair_points_{len(pair_candidates)}"] += 1
        counts["second_root_in_glv_x"] += len(roots) > 1 and roots[1] in orbit_xs(point)
        counts["original_in_pair"] += any(c.x == point.x and c.y == point.y for c in pair_candidates)
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    counts = run(args.samples, args.seed)
    print("Quadratic Lambda-branch experiment")
    print(f"  samples: {args.samples}")
    for key in sorted(counts):
        print(f"  {key}: {counts[key]}")
    print("  Note: branches recover public-point candidates; they do not recover the scalar k.")


if __name__ == "__main__":
    main()
