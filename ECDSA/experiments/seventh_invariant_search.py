#!/usr/bin/env python3
"""Search simple candidate features for a seventh, independent orbit key.

For the full C6 action on y^2=x^3+7, an invariant must be unchanged by
phi(x,y)=(beta*x,y) and i(x,y)=(x,-y).  The expected quotient generator is
u=x^3.  This script checks candidate symmetric/polynomial features and tests
whether they are constant on C6 orbits and whether they vary independently of
u on a finite sample.
"""

import argparse
import random
from collections import Counter

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import (
    DomeInvariantAnalyzer,
    EllipticCurvePoint,
    ENDOMORPHISM_BETA,
    FIELD_PRIME,
    GROUP_ORDER,
)


def point_from_scalar(k):
    q = ec.derive_private_key(k, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(q.x, q.y)


def orbit(point):
    result = []
    current = point
    for _ in range(3):
        result.extend((current, EllipticCurvePoint(current.x, -current.y % FIELD_PRIME)))
        current = EllipticCurvePoint(ENDOMORPHISM_BETA * current.x % FIELD_PRIME, current.y)
    return result


def candidates(point):
    x = point.x
    y = point.y
    x2 = x * x % FIELD_PRIME
    x3 = x * x2 % FIELD_PRIME
    y2 = y * y % FIELD_PRIME
    # All are cheap polynomial/rational-style descriptors.  The curve relation
    # makes y^2 = x^3+7, so several are deliberately expected to collapse.
    return {
        "x3": x3,
        "x6": x3 * x3 % FIELD_PRIME,
        "y2": y2,
        "x3_plus_7": (x3 + 7) % FIELD_PRIME,
        "x3_minus_y2": (x3 - y2) % FIELD_PRIME,
        "x3_times_y2": x3 * y2 % FIELD_PRIME,
        "x6_plus_y2": (x3 * x3 + y2) % FIELD_PRIME,
        "x3_over_y2": x3 * pow(y2, -1, FIELD_PRIME) % FIELD_PRIME,
        "x3_over_x3_plus_7": x3 * pow((x3 + 7) % FIELD_PRIME, -1, FIELD_PRIME) % FIELD_PRIME,
    }


def run(samples, seed):
    rng = random.Random(seed)
    orbit_failures = Counter()
    values_by_u = {name: {} for name in candidates(point_from_scalar(1))}
    seen_pairs = {name: set() for name in values_by_u}
    for _ in range(samples):
        point = point_from_scalar(rng.randrange(1, GROUP_ORDER))
        points = orbit(point)
        per_point = [candidates(q) for q in points]
        u = per_point[0]["x3"]
        for name in values_by_u:
            values = {entry[name] for entry in per_point}
            orbit_failures[name] += len(values) != 1
            value = next(iter(values))
            seen_pairs[name].add((u, value))
            values_by_u[name].setdefault(u, set()).add(value)
    independence = {
        name: sum(len(values) > 1 for values in mapping.values())
        for name, mapping in values_by_u.items()
    }
    return orbit_failures, independence, seen_pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    failures, independent, pairs = run(args.samples, args.seed)
    print("Seventh orbit-invariant search")
    print(f"  samples: {args.samples}")
    print("  candidate: orbit failures; u->multiple-values cases")
    for name in failures:
        print(f"  {name}: {failures[name]}; {independent[name]}")
    print("  Interpretation: a full C6 invariant with zero failures and zero independent cases is a function of u=x^3 on the sample.")


if __name__ == "__main__":
    main()
