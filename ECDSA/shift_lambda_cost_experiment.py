#!/usr/bin/env python3
"""Estimate candidate reduction and recall of a Lambda-sector filter."""

import argparse
import json
from pathlib import Path

from w1 import ENDOMORPHISM_BETA, FIELD_PRIME, GROUP_ORDER, DomeInvariantAnalyzer, EllipticCurvePoint
from w1_circle_uniformity_experiment import point_from_scalar
from shift_lambda_sector_experiment import tiling_shifts


def triangle(point, analyzer):
    if point.is_infinity:
        return ()
    values = []
    for _ in range(3):
        values.append(analyzer.compute_torus_invariant_lambda(point))
        point = EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)
    return tuple(values)


def run(start, width, covered_start, covered_width, sectors, shifts):
    if start < 1 or width <= 0 or start + width > GROUP_ORDER:
        raise ValueError("source interval must lie inside the scalar group")
    sector_width = (covered_width + sectors - 1) // sectors
    analyzer = DomeInvariantAnalyzer()

    def point_for_scalar(scalar):
        return (EllipticCurvePoint(is_infinity=True)
                if scalar == 0 else point_from_scalar(scalar))
    total_full_candidates = 0
    total_vertex_candidates = 0
    true_shift_full_candidates = 0
    true_shift_vertex_candidates = 0
    recalled = 0
    for k in range(start, start + width):
        true_shift = None
        true_sector = None
        per_key_vertex = 0
        per_key_full = len(shifts) * covered_width
        for t in shifts:
            shifted_scalar = (k + t) % GROUP_ORDER
            sectors_seen = {
                value * sectors // FIELD_PRIME
                for value in triangle(point_for_scalar(shifted_scalar), analyzer)
            }
            per_key_vertex += len(sectors_seen) * sector_width
            if covered_start <= shifted_scalar < covered_start + covered_width:
                true_shift = t
                true_sector = ((shifted_scalar - covered_start) * sectors) // covered_width
                true_shift_vertex_candidates += len(sectors_seen) * sector_width
                true_shift_full_candidates += covered_width
                if true_sector in sectors_seen:
                    recalled += 1
        total_full_candidates += per_key_full
        total_vertex_candidates += per_key_vertex
        if true_shift is None:
            raise AssertionError("tiling shifts did not cover a source scalar")
    return {
        "source_start": start,
        "source_width": width,
        "covered_start": covered_start,
        "covered_width": covered_width,
        "sectors": sectors,
        "shift_count": len(shifts),
        "scalar_coverage_fraction": 1.0,
        "all_shift_full_candidates": total_full_candidates,
        "all_shift_vertex_candidates": total_vertex_candidates,
        "all_shift_candidate_fraction": total_vertex_candidates / total_full_candidates,
        "true_shift_full_candidates": true_shift_full_candidates,
        "true_shift_vertex_candidates": true_shift_vertex_candidates,
        "true_shift_candidate_fraction": true_shift_vertex_candidates / true_shift_full_candidates,
        "true_shift_recall": recalled / width,
        "missed_keys": width - recalled,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda value: int(value, 0), default=1 << 12)
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--covered-start", type=lambda value: int(value, 0), default=1)
    parser.add_argument("--covered-width", type=int, default=100)
    parser.add_argument("--sectors", type=int, default=16)
    parser.add_argument("--json", type=Path, default=Path("shift_lambda_cost.json"))
    args = parser.parse_args()
    shifts = tiling_shifts(args.start, args.width, args.covered_start, args.covered_width)
    report = run(args.start, args.width, args.covered_start, args.covered_width,
                 args.sectors, shifts)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"shifts={report['shift_count']}, sectors={args.sectors}")
    print(f"all-shift candidate fraction={report['all_shift_candidate_fraction']:.3%}")
    print(f"true-shift candidate fraction={report['true_shift_candidate_fraction']:.3%}")
    print(f"true-shift recall={report['true_shift_recall']:.3%}, missed={report['missed_keys']}")
    print(f"report={args.json}")


if __name__ == "__main__":
    main()
