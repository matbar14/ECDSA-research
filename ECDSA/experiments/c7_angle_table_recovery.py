#!/usr/bin/env python3
"""Build and query a GLV-compressed table of exact angle signatures.

The signature is the unordered exact triple

    sort(Lambda(P), Lambda(phi(P)), Lambda(phi^2(P))).

It is unchanged for all six points {+/-P, +/-phi(P), +/-phi^2(P)}.  Thus a
table row stores one GLV orbit, not one oriented scalar.  A query public point
recovers the orbit's six scalar candidates; an external interval can reduce
that list to one candidate when it contains only one orbit member.
"""

import argparse
import random
from collections import defaultdict

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import (
    DomeInvariantAnalyzer,
    EllipticCurvePoint,
    ENDOMORPHISM_BETA,
    FIELD_PRIME,
    GROUP_ORDER,
)

GLV_LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72


def point_from_scalar(k):
    q = ec.derive_private_key(k, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(q.x, q.y)


def orbit_scalars(k):
    return {
        (sign * pow(GLV_LAMBDA, power, GROUP_ORDER) * k) % GROUP_ORDER
        for power in range(3)
        for sign in (1, -1)
    }


def mirror_classes(k):
    """Canonical representatives of the three GLV classes modulo negation."""
    return {min(scalar, GROUP_ORDER - scalar) for scalar in orbit_scalars(k)}


def orbit_points(point):
    points = []
    current = point
    for _ in range(3):
        points += [current, EllipticCurvePoint(current.x, (-current.y) % FIELD_PRIME)]
        current = EllipticCurvePoint(ENDOMORPHISM_BETA * current.x % FIELD_PRIME, current.y)
    return points


def angle_signature(analyzer, point, angle_bits=None):
    values = []
    # orbit_points interleaves +Q,-Q; use the three positive representatives.
    for candidate in orbit_points(point)[::2]:
        value = analyzer.compute_torus_invariant_lambda(candidate)
        if value is None:
            raise ValueError("Lambda is undefined on an orbit point")
        values.append(
            value if angle_bits is None
            else value * (1 << angle_bits) // FIELD_PRIME
        )
    return tuple(sorted(values))


def in_interval(scalar, start, size):
    return start <= scalar < start + size


def run(orbit_count, seed, interval_start=None, interval_size=None,
        source_start=None, source_size=None, mirror_canonical=False,
        angle_bits=None):
    analyzer = DomeInvariantAnalyzer()
    rng = random.Random(seed)
    table = defaultdict(set)
    expected_orbits = []

    for _ in range(orbit_count):
        if source_start is None:
            k = rng.randrange(1, GROUP_ORDER)
        else:
            k = source_start + rng.randrange(source_size)
        scalars = mirror_classes(k) if mirror_canonical else orbit_scalars(k)
        signature = angle_signature(analyzer, point_from_scalar(k), angle_bits)
        table[signature].update(scalars)
        expected_orbits.append(scalars)

    query_count = 0
    exact_orbit_hits = 0
    wrong_candidate_sets = 0
    filtered_unique = 0
    filtered_empty = 0
    filtered_multiple = 0
    for scalars in expected_orbits:
        for k in scalars:
            query_point = point_from_scalar(k if not mirror_canonical else min(k, GROUP_ORDER - k))
            recovered = table[angle_signature(analyzer, query_point, angle_bits)]
            query_count += 1
            exact_orbit_hits += recovered == scalars
            wrong_candidate_sets += recovered != scalars
            if interval_start is not None:
                filtered = [candidate for candidate in recovered
                            if in_interval(candidate, interval_start, interval_size)]
                if len(filtered) == 1:
                    filtered_unique += 1
                elif not filtered:
                    filtered_empty += 1
                else:
                    filtered_multiple += 1

    return {
        "orbits": orbit_count,
        "rows": len(table),
        "raw_scalars": orbit_count * 6,
        "queries": query_count,
        "exact_orbit_hits": exact_orbit_hits,
        "wrong_candidate_sets": wrong_candidate_sets,
        "candidate_sizes": sorted({len(values) for values in table.values()}),
        "stored_records": sum(len(values) for values in table.values()),
        "angle_bits": angle_bits,
        "interval_start": interval_start,
        "interval_size": interval_size,
        "filtered_unique": filtered_unique,
        "filtered_empty": filtered_empty,
        "filtered_multiple": filtered_multiple,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbits", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--interval-start", type=int, default=None,
                        help="inclusive scalar interval start for filtering")
    parser.add_argument("--interval-size", type=int, default=None,
                        help="interval width used with --interval-start")
    parser.add_argument("--source-start", type=int, default=None,
                        help="generate table/query seeds from this interval")
    parser.add_argument("--source-size", type=int, default=None,
                        help="width of source interval used with --source-start")
    parser.add_argument("--mirror-canonical", action="store_true",
                        help="store one canonical scalar for each pair {k,n-k}")
    parser.add_argument("--angle-bits", type=int, default=None,
                        help="quantize each of the three angles to this many bits")
    args = parser.parse_args()
    if args.orbits <= 0:
        parser.error("--orbits must be positive")
    if (args.interval_start is None) != (args.interval_size is None):
        parser.error("--interval-start and --interval-size must be supplied together")
    if args.interval_size is not None and args.interval_size <= 0:
        parser.error("--interval-size must be positive")
    if (args.source_start is None) != (args.source_size is None):
        parser.error("--source-start and --source-size must be supplied together")
    if args.source_size is not None and args.source_size <= 0:
        parser.error("--source-size must be positive")
    if args.source_start is not None and (
        args.source_start < 1 or args.source_start + args.source_size > GROUP_ORDER
    ):
        parser.error("source interval must be inside [1, n)")
    if args.angle_bits is not None and not 1 <= args.angle_bits <= 256:
        parser.error("--angle-bits must be in 1..256")
    report = run(args.orbits, args.seed, args.interval_start, args.interval_size,
                 args.source_start, args.source_size, args.mirror_canonical,
                 args.angle_bits)
    print("C7 GLV-compressed exact-angle table")
    print(f"  orbit rows: {report['rows']}")
    print(f"  raw scalar records (six-point baseline): {report['raw_scalars']}")
    print(f"  stored scalar records: {report['stored_records']}")
    print(f"  compression factor by stored records: {report['raw_scalars'] / report['stored_records']:.1f}x")
    print(f"  queries: {report['queries']}")
    print(f"  exact six-candidate orbit recoveries: {report['exact_orbit_hits']}/{report['queries']}")
    print(f"  wrong candidate sets: {report['wrong_candidate_sets']}")
    print(f"  candidate-set sizes in table: {report['candidate_sizes']}")
    if report["angle_bits"] is not None:
        print(f"  quantized angle bits per Lambda: {report['angle_bits']}")
    if args.mirror_canonical:
        print("  mirror mode: one representative min(k, n-k) per sign pair")
    if report["interval_start"] is not None:
        print("  interval filtering:")
        print(f"    unique candidate: {report['filtered_unique']}/{report['queries']}")
        print(f"    no candidate: {report['filtered_empty']}/{report['queries']}")
        print(f"    multiple candidates: {report['filtered_multiple']}/{report['queries']}")
    if args.source_start is not None:
        print(f"  source interval: [{args.source_start}, {args.source_start + args.source_size})")
    print("  Note: the signature recovers an orbit, not an oriented k.")


if __name__ == "__main__":
    main()
