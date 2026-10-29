#!/usr/bin/env python3
"""Search candidate orbit markers beyond K3.

The script distinguishes full C6 invariants from orientation markers.  A full
invariant must be constant on {+/-P,+/-phi(P),+/-phi^2(P)}; an orientation
marker is allowed to take several values there, but should label the six
positions consistently when combined with K3.
"""

import argparse
import math
import random
from collections import Counter, defaultdict

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


def neg(point):
    return EllipticCurvePoint(point.x, (-point.y) % FIELD_PRIME)


def c6_orbit(point):
    result = []
    current = point
    for _ in range(3):
        result.extend((current, neg(current)))
        current = EllipticCurvePoint(ENDOMORPHISM_BETA * current.x % FIELD_PRIME, current.y)
    return result


def lambdas(analyzer, point):
    result = []
    current = point
    for _ in range(3):
        result.append(analyzer.compute_torus_invariant_lambda(current))
        current = EllipticCurvePoint(ENDOMORPHISM_BETA * current.x % FIELD_PRIME, current.y)
    return tuple(result)


def candidate_values(analyzer, point):
    x = point.x
    y = point.y
    l0, l1, l2 = lambdas(analyzer, point)
    u = pow(x, 3, FIELD_PRIME)
    d01 = (l1 - l0) % FIELD_PRIME
    d12 = (l2 - l1) % FIELD_PRIME
    j = (l0*l1 + l1*l2 + l2*l0) % FIELD_PRIME
    inv_j = None if j == 0 else pow(j, -1, FIELD_PRIME)
    # Full invariants and orientation markers.  All arithmetic is exact in F_p.
    return {
        "K3": u,
        "J": j,
        "L0": l0,
        "d01": d01,
        "d12": d12,
        "L0_over_J": None if inv_j is None else l0 * inv_j % FIELD_PRIME,
        "d01_over_J": None if inv_j is None else d01 * inv_j % FIELD_PRIME,
        "x_over_xG": x * pow(BASE_POINT.x, -1, FIELD_PRIME) % FIELD_PRIME,
        "y_over_yG": y * pow(BASE_POINT.y, -1, FIELD_PRIME) % FIELD_PRIME,
        "xy_over_xGyG": (x * y) * pow(BASE_POINT.x * BASE_POINT.y % FIELD_PRIME, -1, FIELD_PRIME) % FIELD_PRIME,
    }


def canonical_signature(values):
    """Canonicalize a candidate over the six positions, if it is scalar-valued."""
    return min(values)


def run(samples, seed):
    rng = random.Random(seed)
    analyzer = DomeInvariantAnalyzer()
    names = tuple(candidate_values(analyzer, point_from_scalar(1)))
    orbit_failures = Counter()
    orbit_distinct = {name: Counter() for name in names}
    random_values = {name: set() for name in names}
    k3_marker_pairs = {name: set() for name in names}

    for _ in range(samples):
        point = point_from_scalar(rng.randrange(1, GROUP_ORDER))
        orbit = c6_orbit(point)
        rows = [candidate_values(analyzer, q) for q in orbit]
        for name in names:
            values = [row[name] for row in rows]
            distinct = len(set(values))
            orbit_distinct[name][distinct] += 1
            orbit_failures[name] += distinct != 1
            random_values[name].add(values[0])
            k3_marker_pairs[name].add((rows[0]["K3"], values[0]))

    report = {}
    for name in names:
        report[name] = {
            "orbit_failures": orbit_failures[name],
            "orbit_distinct_histogram": dict(sorted(orbit_distinct[name].items())),
            "distinct_random_values": len(random_values[name]),
            "distinct_k3_marker_pairs": len(k3_marker_pairs[name]),
        }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    report = run(args.samples, args.seed)
    print("Candidate invariant/marker search")
    print(f"  independent orbits: {args.samples}")
    print("  name              orbit failures  orbit values  random values  (K3,marker) pairs")
    for name, result in report.items():
        print(
            f"  {name:16s} {result['orbit_failures']:14d}"
            f"  {str(result['orbit_distinct_histogram']):12s}"
            f"  {result['distinct_random_values']:13d}"
            f"  {result['distinct_k3_marker_pairs']:17d}"
        )
    print("  Full C6 invariant: orbit failures must be 0.")
    print("  Useful orientation marker: orbit values should distinguish positions while (K3,marker) stays collision-free.")


if __name__ == "__main__":
    main()
