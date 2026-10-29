#!/usr/bin/env python3
"""Measure public-key parity and public-point shifts for Lambda coverage."""

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from w1 import (
    BASE_POINT, DomeInvariantAnalyzer, ENDOMORPHISM_BETA,
    FIELD_PRIME, GROUP_ORDER, EllipticCurvePoint,
)
from w1_circle_uniformity_experiment import point_from_scalar


def lambda_triangle(point, analyzer):
    values = []
    for _ in range(3):
        values.append(analyzer.compute_torus_invariant_lambda(point))
        point = EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)
    return tuple(values)


def neg(point):
    return EllipticCurvePoint(point.x, (-point.y) % FIELD_PRIME)


def signed_delta(current, previous, width):
    value = current - previous
    if value > width // 2:
        value -= width
    if value < -width // 2:
        value += width
    return value


def run(start, width, sectors, shift_radius):
    analyzer = DomeInvariantAnalyzer()
    point = point_from_scalar(start)
    records = []
    parity_counts = Counter()
    parity_hits = Counter()
    parity_steps = defaultdict(list)
    negation_mismatches = 0
    previous_positions = None
    for offset in range(width):
        scalar = start + offset
        triangle = lambda_triangle(point, analyzer)
        projected = tuple(value * width // FIELD_PRIME for value in triangle)
        parity = int(point.y & 1)
        parity_counts[parity] += 1
        scalar_sector = offset * sectors // width
        vertex_sectors = tuple(position * sectors // width for position in projected)
        hits = sum(sector == scalar_sector for sector in vertex_sectors)
        parity_hits[parity] += hits > 0
        if previous_positions is not None:
            parity_steps[parity].extend(
                signed_delta(projected[i], previous_positions[i], width)
                for i in range(3)
            )
        previous_positions = projected
        if triangle != lambda_triangle(neg(point), analyzer):
            negation_mismatches += 1
        records.append({"scalar": scalar, "parity": parity, "hits": hits})
        point = point + BASE_POINT

    shift_target = start + width // 2
    shift_point = point_from_scalar(shift_target)
    shifts = []
    for shift in range(-shift_radius, shift_radius + 1):
        shifted_scalar = (shift_target + shift) % GROUP_ORDER
        shifted_point = shift_point
        if shift >= 0:
            for _ in range(shift):
                shifted_point = shifted_point + BASE_POINT
        else:
            inverse = EllipticCurvePoint(BASE_POINT.x, (-BASE_POINT.y) % FIELD_PRIME)
            for _ in range(-shift):
                shifted_point = shifted_point + inverse
        expected = point_from_scalar(shifted_scalar)
        shifts.append({
            "shift": shift,
            "scalar": shifted_scalar,
            "inside_window": start <= shifted_scalar < start + width,
            "point_matches_scalar_shift": analyzer.points_equal(shifted_point, expected),
            "compressed_parity": int(shifted_point.y & 1),
        })
    return {
        "start": start,
        "width": width,
        "sectors": sectors,
        "parity_counts": dict(parity_counts),
        "parity_hit_counts": dict(parity_hits),
        "parity_hit_rates": {
            str(parity): parity_hits[parity] / parity_counts[parity]
            for parity in parity_counts
        },
        "random_any_hit_baseline": 1 - (1 - 1 / sectors) ** 3,
        "negation_lambda_triangle_mismatches": negation_mismatches,
        "parity_step_direction": {
            str(parity): {
                "positive": sum(value > 0 for value in values) / len(values),
                "negative": sum(value < 0 for value in values) / len(values),
                "zero": sum(value == 0 for value in values) / len(values),
            }
            for parity, values in parity_steps.items()
        },
        "shift_target": shift_target,
        "shift_radius": shift_radius,
        "shift_count": len(shifts),
        "shift_inside_window_count": sum(item["inside_window"] for item in shifts),
        "shift_point_check_failures": sum(not item["point_matches_scalar_shift"] for item in shifts),
        "shifts": shifts,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda value: int(value, 0), default=1 << 134)
    parser.add_argument("--width", type=int, default=10000)
    parser.add_argument("--sectors", type=int, default=32)
    parser.add_argument("--shift-radius", type=int, default=100)
    parser.add_argument("--json", type=Path, default=Path("parity_shift_experiment.json"))
    args = parser.parse_args()
    report = run(args.start, args.width, args.sectors, args.shift_radius)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"range=[{args.start}, {args.start + args.width}), sectors={args.sectors}")
    print(f"parity counts={report['parity_counts']}")
    print(f"parity hit rates={report['parity_hit_rates']} "
          f"(random any-hit baseline={report['random_any_hit_baseline']:.3%})")
    print(f"Lambda(P) vs Lambda(-P) mismatches={report['negation_lambda_triangle_mismatches']}")
    print(f"parity step directions={report['parity_step_direction']}")
    print(f"shifts={report['shift_count']}, inside window={report['shift_inside_window_count']}, "
          f"point check failures={report['shift_point_check_failures']}")
    print(f"report={args.json}")


if __name__ == "__main__":
    main()
