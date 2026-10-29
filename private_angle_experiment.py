#!/usr/bin/env python3
"""Direct sampling experiment: private scalar -> P -> Lambda -> unit angle.

For each distinct 256-bit private scalar k, the script builds P = kG and
computes Lambda(P) with DomeInvariantAnalyzer from w1.py.  The canonical field
residue a = Lambda(P) is the exact index of a point on the unit circle:

    theta = 2*pi*a/p.

The experiment records exact residue collisions, not equality of floating
point angles.  It also bins residues into equal arcs to inspect their sampled
distribution around the circle.
"""

import math
import random
import time

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import DomeInvariantAnalyzer, EllipticCurvePoint, FIELD_PRIME, GROUP_ORDER


SAMPLE_COUNT = 100000
SEED = 20261014
ARC_COUNT = 256
PROGRESS_INTERVAL = 10000


def generate_unique_256_bit_scalar(generator, used_scalars):
    while True:
        scalar = generator.randrange(1 << 255, GROUP_ORDER)
        if scalar not in used_scalars:
            used_scalars.add(scalar)
            return scalar


def point_from_scalar(scalar):
    public_numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()
    ).public_key().public_numbers()
    return EllipticCurvePoint(public_numbers.x, public_numbers.y)


def angle_radians(lambda_value):
    return math.tau * (lambda_value / FIELD_PRIME)


def run_experiment(sample_count=SAMPLE_COUNT, seed=SEED, collect_residues=False):
    if not isinstance(sample_count, int) or sample_count <= 0:
        raise ValueError("sample_count must be a positive integer")

    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    used_scalars = set()
    first_scalar_by_angle = {}
    arc_counts = [0] * ARC_COUNT
    boundary_angle_count = 0
    collision_pairs = []
    start_time = time.perf_counter()

    for sample_index in range(1, sample_count + 1):
        scalar = generate_unique_256_bit_scalar(generator, used_scalars)
        point = point_from_scalar(scalar)
        lambda_value = analyzer.compute_torus_invariant_lambda(point)
        if not 0 <= lambda_value < FIELD_PRIME:
            raise RuntimeError("Lambda residue is outside the canonical field interval")

        theta = angle_radians(lambda_value)
        if not 0.0 <= theta <= math.tau:
            raise RuntimeError("Unit-circle angle is outside the expected interval")
        boundary_angle_count += theta == math.tau

        if lambda_value in first_scalar_by_angle:
            collision_pairs.append((first_scalar_by_angle[lambda_value], scalar))
        else:
            first_scalar_by_angle[lambda_value] = scalar
        arc_counts[lambda_value * ARC_COUNT // FIELD_PRIME] += 1

        if sample_index % PROGRESS_INTERVAL == 0:
            print(f"  Completed: {sample_index}/{sample_count}")

    expected = sample_count / ARC_COUNT
    chi_square = math.fsum(
        (count - expected) ** 2 / expected for count in arc_counts
    )
    result = {
        "sample_count": sample_count,
        "seed": seed,
        "distinct_private_scalars": len(used_scalars),
        "distinct_angle_indices": len(first_scalar_by_angle),
        "angle_collision_count": len(collision_pairs),
        "boundary_angle_count": boundary_angle_count,
        "arc_count": ARC_COUNT,
        "arc_minimum": min(arc_counts),
        "arc_maximum": max(arc_counts),
        "arc_chi_square": chi_square,
        "elapsed_seconds": time.perf_counter() - start_time,
    }
    if collect_residues:
        result["angle_residues"] = set(first_scalar_by_angle)
        result["arc_counts"] = arc_counts
    return result


def print_report(report):
    print("Private-scalar to unit-angle experiment")
    print(f"  Distinct 256-bit private scalars: {report['distinct_private_scalars']}")
    print(f"  Reproducibility seed: {report['seed']}")
    print(f"  Distinct exact angle indices Lambda(P): {report['distinct_angle_indices']}")
    print(f"  Exact angle collisions: {report['angle_collision_count']}")
    print(f"  Angles rounded to 2*pi by binary64: {report['boundary_angle_count']}")
    print(f"  Equal-circle arcs: {report['arc_count']}")
    print(f"  Arc occupancy range: {report['arc_minimum']}..{report['arc_maximum']}")
    print(f"  Pearson chi-square for arc counts: {report['arc_chi_square']:.6f}")
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
