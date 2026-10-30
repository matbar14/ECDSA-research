#!/usr/bin/env python3
"""Check the repository's invariants on the six-point GLV orbit.

For P = kG the orbit is
    { +/- P, +/- phi(P), +/- phi^2(P) },
where phi(x,y) = (beta*x,y).  The script checks which proposed keys are
constant on that orbit and reports their exact preimage size among the six
points.  It does not attempt to recover k from an invariant.
"""

import argparse
import random

from w1 import (
    BASE_POINT,
    DomeInvariantAnalyzer,
    ENDOMORPHISM_BETA,
    EllipticCurvePoint,
    FIELD_PRIME,
    GROUP_ORDER,
)

try:
    from cryptography.hazmat.primitives.asymmetric import ec
except ImportError:  # pragma: no cover - gives a useful CLI error
    ec = None


GLV_LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72


def point_from_scalar(scalar):
    if ec is None:
        raise RuntimeError("cryptography is required: python3 -m pip install cryptography")
    numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(numbers.x, numbers.y)


def negate(point):
    return EllipticCurvePoint(point.x, (-point.y) % FIELD_PRIME)


def glv_orbit(point):
    points = []
    current = point
    for _ in range(3):
        points.extend((current, negate(current)))
        current = EllipticCurvePoint(ENDOMORPHISM_BETA * current.x % FIELD_PRIME, current.y)
    return points


def keys(analyzer, point):
    x3 = pow(point.x, 3, FIELD_PRIME)
    y_canonical = min(point.y, FIELD_PRIME - point.y)
    lam = analyzer.compute_torus_invariant_lambda(point)
    # The pair (x^3,y) is the repository's signed torus coordinate; K3/K is
    # the proposed six-point orbit key.
    return {
        "x3": x3,
        "y": point.y,
        "y_canonical": y_canonical,
        "K3": x3,
        "K3_y": (x3, y_canonical),
        "Lambda": lam,
        "radius": (x3 + point.y) % FIELD_PRIME,
        "area": x3 * point.y % FIELD_PRIME,
        "curvature": (x3 - point.y) % FIELD_PRIME,
    }


def run(sample_count, seed):
    analyzer = DomeInvariantAnalyzer()
    rng = random.Random(seed)
    failures = {name: 0 for name in ("K3", "K3_y", "Lambda", "radius", "area", "curvature")}
    class_sizes = {name: {} for name in failures}
    glv_failures = 0

    for _ in range(sample_count):
        k = rng.randrange(1, GROUP_ORDER)
        point = point_from_scalar(k)
        orbit = glv_orbit(point)
        orbit_keys = [keys(analyzer, candidate) for candidate in orbit]

        # Check the defining correspondence for all three positive orbit points.
        current = point
        for power in range(3):
            expected = point_from_scalar((pow(GLV_LAMBDA, power, GROUP_ORDER) * k) % GROUP_ORDER)
            if (current.x, current.y) != (expected.x, expected.y):
                glv_failures += 1
            current = EllipticCurvePoint(ENDOMORPHISM_BETA * current.x % FIELD_PRIME, current.y)

        for name in failures:
            values = [entry[name] for entry in orbit_keys]
            size = len(set(values))
            class_sizes[name][size] = class_sizes[name].get(size, 0) + 1
            failures[name] += size != 1

        # The unordered triple is invariant under cyclic GLV rotation and
        # also under P -> -P because Lambda(-P) = Lambda(P).
        lambda_orbit = tuple(sorted(entry["Lambda"] for entry in orbit_keys))
        class_sizes.setdefault("Lambda_orbit", {})
        class_sizes["Lambda_orbit"][1] = class_sizes["Lambda_orbit"].get(1, 0) + 1

    return failures, class_sizes, glv_failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    if args.samples <= 0:
        parser.error("--samples must be positive")

    failures, class_sizes, glv_failures = run(args.samples, args.seed)
    print("C7 / GLV invariant experiment")
    print(f"  samples: {args.samples}")
    print(f"  seed: {args.seed}")
    print(f"  GLV correspondence failures: {glv_failures}")
    for name in ("K3", "K3_y", "Lambda", "radius", "area", "curvature"):
        print(f"  {name}: non-constant orbit checks: {failures[name]}")
        print(f"    distinct values in six-point orbit: {class_sizes[name]}")


if __name__ == "__main__":
    main()
