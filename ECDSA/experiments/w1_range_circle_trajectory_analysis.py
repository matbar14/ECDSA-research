#!/usr/bin/env python3
"""Quantify Lambda-triangle motion on a circle matching a finite range."""

import argparse
import csv
import json
import math
import statistics
from collections import Counter
from pathlib import Path

from w1 import BASE_POINT, DomeInvariantAnalyzer, ENDOMORPHISM_BETA, FIELD_PRIME, GROUP_ORDER, EllipticCurvePoint
from w1_circle_uniformity_experiment import point_from_scalar


def values(point, analyzer):
    result = []
    for _ in range(3):
        result.append(analyzer.compute_torus_invariant_lambda(point))
        point = EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)
    return tuple(result)


def signed_delta(current, previous, width):
    delta = current - previous
    if delta > width // 2:
        delta -= width
    elif delta < -width // 2:
        delta += width
    return delta


def run(start, width, sectors, output_csv, output_json):
    if start < 1 or width < 4 or start + width > GROUP_ORDER:
        raise ValueError("use a finite interval contained in the scalar range")
    analyzer = DomeInvariantAnalyzer()
    point = point_from_scalar(start)
    positions = []
    for _ in range(width):
        triangle = values(point, analyzer)
        positions.append(tuple(value * width // FIELD_PRIME for value in triangle))
        point = point + BASE_POINT

    direction_rows = []
    for i in range(1, width):
        deltas = tuple(signed_delta(positions[i][j], positions[i - 1][j], width) for j in range(3))
        direction_rows.append(deltas)
    flattened = [delta for row in direction_rows for delta in row]
    hits = []
    for offset, triangle in enumerate(positions):
        scalar_sector = offset * sectors // width
        vertex_sectors = tuple(position * sectors // width for position in triangle)
        hits.append(sum(sector == scalar_sector for sector in vertex_sectors))
    sector_sequences = [tuple(position * sectors // width for position in triangle) for triangle in positions]
    return_counts = {}
    for lag in range(1, min(width // 2, 1001)):
        return_counts[lag] = sum(sector_sequences[i] == sector_sequences[i - lag]
                                 for i in range(lag, width)) / (width - lag)

    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["offset", "lambda0", "lambda1", "lambda2", "scalar_sector", "hit_count"])
        for offset, triangle in enumerate(positions):
            scalar_sector = offset * sectors // width
            writer.writerow([offset, *triangle, scalar_sector, hits[offset]])

    direction_summary = {
        "positive_fraction": sum(delta > 0 for delta in flattened) / len(flattened),
        "negative_fraction": sum(delta < 0 for delta in flattened) / len(flattened),
        "zero_fraction": sum(delta == 0 for delta in flattened) / len(flattened),
        "mean_signed_delta": statistics.fmean(flattened),
        "median_signed_delta": statistics.median(flattened),
        "per_vertex": [
            {
                "positive_fraction": sum(row[j] > 0 for row in direction_rows) / len(direction_rows),
                "negative_fraction": sum(row[j] < 0 for row in direction_rows) / len(direction_rows),
                "mean_delta": statistics.fmean(row[j] for row in direction_rows),
                "median_delta": statistics.median(row[j] for row in direction_rows),
            }
            for j in range(3)
        ],
    }
    hit_distribution = {str(value): hits.count(value) for value in range(4)}
    report = {
        "start": start, "width": width, "sectors": sectors,
        "circumference": width,
        "direction": direction_summary,
        "hit_distribution": hit_distribution,
        "hit_mean": statistics.fmean(hits),
        "hit_any_fraction": sum(value > 0 for value in hits) / len(hits),
        "hit_two_or_more_fraction": sum(value >= 2 for value in hits) / len(hits),
        "random_vertex_hit_baseline": 1 / sectors,
        "random_any_hit_baseline": 1 - (1 - 1 / sectors) ** 3,
        "random_two_or_more_baseline": 3 * (1 / sectors) ** 2 * (1 - 1 / sectors) + (1 / sectors) ** 3,
        "same_triangle_sector_fraction": sum(len(set(row)) == 1 for row in sector_sequences) / len(sector_sequences),
        "sector_return_rates": {str(k): v for k, v in return_counts.items() if k <= max(10, sectors * 3)},
        "csv": str(output_csv),
    }
    output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda value: int(value, 0), default=1 << 134)
    parser.add_argument("--width", type=int, default=10000)
    parser.add_argument("--sectors", type=int, default=32)
    parser.add_argument("--csv", type=Path, default=Path("results/w1_range_circle_trajectory.csv"))
    parser.add_argument("--json", type=Path, default=Path("results/w1_range_circle_trajectory.json"))
    args = parser.parse_args()
    report = run(args.start, args.width, args.sectors, args.csv, args.json)
    direction = report["direction"]
    print(f"range=[{args.start}, {args.start + args.width}), circumference={args.width}, sectors={args.sectors}")
    print(f"motion signs: +={direction['positive_fraction']:.3%}, "
          f"-={direction['negative_fraction']:.3%}, 0={direction['zero_fraction']:.3%}")
    for index, row in enumerate(direction["per_vertex"]):
        print(f"vertex {index}: +={row['positive_fraction']:.3%}, -={row['negative_fraction']:.3%}, "
              f"mean delta={row['mean_delta']:.3f}, median={row['median_delta']:.3f}")
    print(f"hit distribution={report['hit_distribution']}")
    print(f"any hit={report['hit_any_fraction']:.3%} (random {report['random_any_hit_baseline']:.3%}), "
          f"two+={report['hit_two_or_more_fraction']:.3%} (random {report['random_two_or_more_baseline']:.3%})")
    print(f"same triangle sector={report['same_triangle_sector_fraction']:.3%}")
    print(f"report={args.json}; csv={args.csv}")


if __name__ == "__main__":
    main()
