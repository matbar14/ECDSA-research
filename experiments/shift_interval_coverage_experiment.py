#!/usr/bin/env python3
"""Exhaustively test coverage obtained from translated scalar windows.

For Q = kG and Q_t = Q + tG = (k+t)G, a covered window S catches k
exactly when k belongs to S-t.  This experiment enumerates a finite source
range, checks a precomputed shift set, and verifies the point identity on a
sample of successful translations.
"""

import argparse
import json
import random
from pathlib import Path

from w1 import BASE_POINT, FIELD_PRIME, GROUP_ORDER, EllipticCurvePoint
from w1_circle_uniformity_experiment import point_from_scalar


def neg(point):
    if point.is_infinity:
        return point
    return EllipticCurvePoint(point.x, (-point.y) % FIELD_PRIME)


def shifts_for_tiling(start, width, covered_start, covered_width):
    """Return shifts whose translated windows tile the source interval."""
    count = (width + covered_width - 1) // covered_width
    return [covered_start - (start + i * covered_width) for i in range(count)]


def run(start, width, covered_start, covered_width, shifts, checks, seed):
    if start < 1 or width <= 0 or start + width > GROUP_ORDER:
        raise ValueError("source interval must lie inside the scalar group")
    if covered_width <= 0:
        raise ValueError("covered width must be positive")

    hits_by_k = []
    selected = []
    for k in range(start, start + width):
        hits = [t for t in shifts
                if covered_start <= (k + t) % GROUP_ORDER < covered_start + covered_width]
        hits_by_k.append(len(hits))
        selected.append(hits[0] if hits else None)

    rng = random.Random(seed)
    candidates = [i for i, t in enumerate(selected) if t is not None]
    sample_indices = candidates[:checks]
    if len(candidates) > checks:
        sample_indices = rng.sample(candidates, checks)
    point_failures = 0
    recovery_failures = 0
    for index in sample_indices:
        k = start + index
        t = selected[index]
        q = point_from_scalar(k)
        shifted = q + (t * BASE_POINT if t >= 0 else (-t) * neg(BASE_POINT))
        shifted_scalar = (k + t) % GROUP_ORDER
        expected = point_from_scalar(shifted_scalar)
        if (shifted.x, shifted.y, shifted.is_infinity) != (expected.x, expected.y, expected.is_infinity):
            point_failures += 1
        recovered = (shifted_scalar - t) % GROUP_ORDER
        recovered_point = point_from_scalar(recovered)
        same_point = (recovered_point.x, recovered_point.y, recovered_point.is_infinity) == (
            q.x, q.y, q.is_infinity
        )
        if recovered != k or not same_point:
            recovery_failures += 1

    uncovered = sum(value == 0 for value in hits_by_k)
    report = {
        "source_start": start,
        "source_width": width,
        "covered_start": covered_start,
        "covered_width": covered_width,
        "shift_count": len(shifts),
        "shifts": shifts,
        "covered_keys": width - uncovered,
        "uncovered_keys": uncovered,
        "coverage_fraction": (width - uncovered) / width,
        "minimum_hits_per_key": min(hits_by_k),
        "maximum_hits_per_key": max(hits_by_k),
        "point_identity_failures": point_failures,
        "recovery_failures": recovery_failures,
        "verification_samples": len(sample_indices),
        "seed": seed,
    }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda value: int(value, 0), default=1 << 12)
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--covered-start", type=lambda value: int(value, 0), default=1)
    parser.add_argument("--covered-width", type=int, default=100)
    parser.add_argument("--shifts", type=int, nargs="*", default=None)
    parser.add_argument("--checks", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--json", type=Path, default=Path("results/shift_interval_coverage.json"))
    args = parser.parse_args()
    shifts = args.shifts
    if shifts is None:
        shifts = shifts_for_tiling(args.start, args.width, args.covered_start, args.covered_width)
    report = run(args.start, args.width, args.covered_start, args.covered_width,
                 shifts, args.checks, args.seed)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"source=[{args.start}, {args.start + args.width})")
    print(f"covered=[{args.covered_start}, {args.covered_start + args.covered_width})")
    print(f"shifts={report['shift_count']}, coverage={report['coverage_fraction']:.3%}, "
          f"hits/key={report['minimum_hits_per_key']}..{report['maximum_hits_per_key']}")
    print(f"point identity failures={report['point_identity_failures']}, "
          f"recovery failures={report['recovery_failures']}")
    print(f"report={args.json}")


if __name__ == "__main__":
    main()
