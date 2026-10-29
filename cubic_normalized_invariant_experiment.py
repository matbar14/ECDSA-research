#!/usr/bin/env python3
"""Test whether cubic normalized Lambda invariants add information.

For L0+L1+L2=0 the repository derives I2=I3=J and
D=-J^2*(4J+27).  This script verifies the dependency on sampled points and
counts distinct exact keys.  Ratios such as I3/I2 are projective normalizations
and are expected to be constant, hence they cannot identify k.
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


def lambda_triple(analyzer, point):
    values = []
    current = point
    for _ in range(3):
        values.append(analyzer.compute_torus_invariant_lambda(current))
        current = EllipticCurvePoint(
            ENDOMORPHISM_BETA * current.x % FIELD_PRIME, current.y
        )
    return tuple(values)


def invariants(values):
    l0, l1, l2 = values
    i2 = (l0*l1 + l1*l2 + l2*l0) % FIELD_PRIME
    i3 = l0*l1*l2 % FIELD_PRIME
    d = ((l0-l1)**2 * (l1-l2)**2 * (l2-l0)**2) % FIELD_PRIME
    ratio_i3_i2 = None if i2 == 0 else i3 * pow(i2, -1, FIELD_PRIME) % FIELD_PRIME
    ratio_d_i2sq = None if i2 == 0 else d * pow(i2*i2 % FIELD_PRIME, -1, FIELD_PRIME) % FIELD_PRIME
    return i2, i3, d, ratio_i3_i2, ratio_d_i2sq


def run(samples, seed):
    rng = random.Random(seed)
    analyzer = DomeInvariantAnalyzer()
    rows = []
    failures = Counter()
    for _ in range(samples):
        scalar = rng.randrange(1, GROUP_ORDER)
        values = lambda_triple(analyzer, point_from_scalar(scalar))
        i2, i3, d, r13, rd = invariants(values)
        expected_d = (-i2*i2*(4*i2 + 27)) % FIELD_PRIME
        failures["sum"] += sum(values) % FIELD_PRIME != 0
        failures["i2_i3"] += i2 != i3
        failures["d_relation"] += d != expected_d
        rows.append((scalar, i2, i3, d, r13, rd))
    return rows, failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    rows, failures = run(args.samples, args.seed)
    keys = {
        "J=I2": len({row[1] for row in rows}),
        "(I2,I3,D)": len({row[1:4] for row in rows}),
        "I3/I2": len({row[4] for row in rows}),
        "D/I2^2": len({row[5] for row in rows}),
        "(J,D/I2^2)": len({(row[1], row[5]) for row in rows}),
    }
    print("Cubic normalized Lambda invariant experiment")
    print(f"  samples: {args.samples}")
    print(f"  identity failures: {dict(failures)}")
    for name, count in keys.items():
        print(f"  distinct {name}: {count}")
    print(f"  distinct normalized ratio I3/I2 values: {sorted({row[4] for row in rows})[:5]}")
    print("  Interpretation: normalized ratios discard J; they are not extra coordinates.")


if __name__ == "__main__":
    main()
