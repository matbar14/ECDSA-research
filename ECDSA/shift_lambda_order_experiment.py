#!/usr/bin/env python3
"""Measure a lossless sector ordering based on Lambda vertices."""

import argparse
import json
from pathlib import Path

from w1 import ENDOMORPHISM_BETA, FIELD_PRIME, GROUP_ORDER, DomeInvariantAnalyzer, EllipticCurvePoint
from w1_circle_uniformity_experiment import point_from_scalar
from shift_lambda_sector_experiment import tiling_shifts


def triangle(point, analyzer):
    if point.is_infinity:
        return []
    values = []
    for _ in range(3):
        values.append(analyzer.compute_torus_invariant_lambda(point))
        point = EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)
    return values


def run(start, width, covered_start, covered_width, sectors, shifts):
    sector_width = (covered_width + sectors - 1) // sectors
    analyzer = DomeInvariantAnalyzer()
    ordered_total = 0
    full_total = 0
    baseline_ordered_total = 0
    recalled = 0
    for k in range(start, start + width):
        true_shift = None
        true_sector = None
        true_prediction = set()
        per_key_ordered = (len(shifts) - 1) * covered_width
        per_key_full = len(shifts) * covered_width
        for t in shifts:
            shifted = (k + t) % GROUP_ORDER
            sectors_seen = {
                value * sectors // FIELD_PRIME
                for value in triangle(
                    EllipticCurvePoint(is_infinity=True)
                    if shifted == 0 else point_from_scalar(shifted), analyzer
                )
            }
            if covered_start <= shifted < covered_start + covered_width:
                true_shift = t
                true_sector = ((shifted - covered_start) * sectors) // covered_width
                true_prediction = sectors_seen
        if true_shift is None:
            raise AssertionError("tiling shifts did not cover a source scalar")
        predicted = sorted(true_prediction)
        order = predicted + [sector for sector in range(sectors) if sector not in true_prediction]
        rank = order.index(true_sector) + 1
        per_key_ordered += min(rank * sector_width, covered_width)
        baseline_rank = true_sector + 1
        baseline_ordered_total += (len(shifts) - 1) * covered_width + baseline_rank * sector_width
        ordered_total += per_key_ordered
        full_total += per_key_full
        recalled += true_sector in true_prediction
    return {
        "source_start": start,
        "source_width": width,
        "covered_start": covered_start,
        "covered_width": covered_width,
        "sectors": sectors,
        "shift_count": len(shifts),
        "full_candidates": full_total,
        "ordered_candidates": ordered_total,
        "ordered_fraction_of_full": ordered_total / full_total,
        "ordered_savings": 1 - ordered_total / full_total,
        "canonical_baseline_candidates": baseline_ordered_total,
        "canonical_baseline_fraction": baseline_ordered_total / full_total,
        "prediction_recall": recalled / width,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda value: int(value, 0), default=1 << 12)
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--covered-start", type=lambda value: int(value, 0), default=1)
    parser.add_argument("--covered-width", type=int, default=100)
    parser.add_argument("--sectors", type=int, default=16)
    parser.add_argument("--json", type=Path, default=Path("shift_lambda_order.json"))
    args = parser.parse_args()
    shifts = tiling_shifts(args.start, args.width, args.covered_start, args.covered_width)
    report = run(args.start, args.width, args.covered_start, args.covered_width,
                 args.sectors, shifts)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"shifts={report['shift_count']}, sectors={args.sectors}")
    print(f"ordered candidates={report['ordered_fraction_of_full']:.3%} of full, "
          f"savings={report['ordered_savings']:.3%}")
    print(f"prediction recall={report['prediction_recall']:.3%}")
    print(f"report={args.json}")


if __name__ == "__main__":
    main()
