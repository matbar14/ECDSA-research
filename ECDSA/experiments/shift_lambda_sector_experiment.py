#!/usr/bin/env python3
"""Measure Lambda-sector hits after translating public targets by known shifts.

The source interval is tiled into translated copies of a smaller covered
window.  For every source scalar we evaluate the Lambda triangle of each
translated public point and compare its vertices with the scalar sector inside
the covered window.  The experiment reports exact scalar coverage separately
from the invariant-based hit rate.
"""

import argparse
import json
import random
from pathlib import Path

from w1 import (
    BASE_POINT,
    ENDOMORPHISM_BETA,
    FIELD_PRIME,
    GROUP_ORDER,
    DomeInvariantAnalyzer,
    EllipticCurvePoint,
)
from w1_circle_uniformity_experiment import point_from_scalar


def neg(point):
    return EllipticCurvePoint(point.x, (-point.y) % FIELD_PRIME)


def triangle(point, analyzer):
    values = []
    for _ in range(3):
        values.append(analyzer.compute_torus_invariant_lambda(point))
        point = EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)
    return tuple(values)


def tiling_shifts(start, width, covered_start, covered_width):
    count = (width + covered_width - 1) // covered_width
    return [covered_start - (start + i * covered_width) for i in range(count)]


def run(start, width, covered_start, covered_width, sectors, shifts, checks, seed):
    if start < 1 or width <= 0 or start + width > GROUP_ORDER:
        raise ValueError("source interval must lie inside the scalar group")
    if covered_start < 1 or covered_width <= 0:
        raise ValueError("covered window must contain positive scalars")
    analyzer = DomeInvariantAnalyzer()
    triangle_hits = 0
    inside_keys = 0
    any_shift_hits = 0
    shift_hit_counts = []
    invariant_hits_by_key = []
    point_identity_failures = 0
    point_check_indices = set(random.Random(seed).sample(range(width), min(checks, width)))

    for index, k in enumerate(range(start, start + width)):
        inside = []
        invariant_hit = False
        for t in shifts:
            shifted_scalar = (k + t) % GROUP_ORDER
            in_window = covered_start <= shifted_scalar < covered_start + covered_width
            if not in_window:
                continue
            inside.append(t)
            target_sector = ((shifted_scalar - covered_start) * sectors) // covered_width
            positions = tuple(value * sectors // FIELD_PRIME
                              for value in triangle(point_from_scalar(shifted_scalar), analyzer))
            invariant_hit |= any(position == target_sector for position in positions)
            if index in point_check_indices:
                q = point_from_scalar(k)
                shifted = q + (t * BASE_POINT if t >= 0 else (-t) * neg(BASE_POINT))
                expected = point_from_scalar(shifted_scalar)
                if (shifted.x, shifted.y, shifted.is_infinity) != (
                        expected.x, expected.y, expected.is_infinity):
                    point_identity_failures += 1
        inside_keys += bool(inside)
        any_shift_hits += invariant_hit
        shift_hit_counts.append(len(inside))
        invariant_hits_by_key.append(invariant_hit)

    report = {
        "source_start": start,
        "source_width": width,
        "covered_start": covered_start,
        "covered_width": covered_width,
        "sectors": sectors,
        "shift_count": len(shifts),
        "shifts": shifts,
        "scalar_coverage_fraction": inside_keys / width,
        "triangle_sector_hit_fraction": any_shift_hits / width,
        "triangle_hit_baseline_per_window": min(1.0, 3 / sectors),
        "minimum_inside_shifts_per_key": min(shift_hit_counts),
        "maximum_inside_shifts_per_key": max(shift_hit_counts),
        "point_identity_failures": point_identity_failures,
        "verification_samples": len(point_check_indices),
        "seed": seed,
    }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda value: int(value, 0), default=1 << 12)
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--covered-start", type=lambda value: int(value, 0), default=1)
    parser.add_argument("--covered-width", type=int, default=100)
    parser.add_argument("--sectors", type=int, default=16)
    parser.add_argument("--shifts", type=int, nargs="*", default=None)
    parser.add_argument("--checks", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--json", type=Path, default=Path("results/shift_lambda_sector.json"))
    args = parser.parse_args()
    shifts = args.shifts
    if shifts is None:
        shifts = tiling_shifts(args.start, args.width, args.covered_start, args.covered_width)
    report = run(args.start, args.width, args.covered_start, args.covered_width,
                 args.sectors, shifts, args.checks, args.seed)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"source=[{args.start}, {args.start + args.width})")
    print(f"covered=[{args.covered_start}, {args.covered_start + args.covered_width}), sectors={args.sectors}")
    print(f"shifts={report['shift_count']}, scalar coverage={report['scalar_coverage_fraction']:.3%}, "
          f"triangle hit={report['triangle_sector_hit_fraction']:.3%}, "
          f"baseline={report['triangle_hit_baseline_per_window']:.3%}")
    print(f"point identity failures={report['point_identity_failures']}")
    print(f"report={args.json}")


if __name__ == "__main__":
    main()
