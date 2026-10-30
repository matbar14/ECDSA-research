#!/usr/bin/env python3
"""Trace Lambda triangles while subtracting G down to the zero point.

The trace starts at P=N*G and updates the point by P <- P-G. For every
nonzero scalar, it measures circular triangle collapse before and after a
stepwise rotation by -Lambda(G). The terminal 0*G point is recorded as the
point at infinity; Lambda is intentionally undefined there.
"""

import argparse
import csv
import json
import math
import statistics
from pathlib import Path

from w1 import (
    BASE_POINT, DomeInvariantAnalyzer, ENDOMORPHISM_BETA,
    FIELD_PRIME, EllipticCurvePoint,
)
from w1_circle_uniformity_experiment import point_from_scalar


def triangle_values(point, analyzer):
    if point.is_infinity:
        return None
    values = []
    for _ in range(3):
        values.append(analyzer.compute_torus_invariant_lambda(point))
        point = EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)
    return tuple(values)


def gaps(values, circumference):
    ordered = sorted(values)
    result = (
        ordered[1] - ordered[0],
        ordered[2] - ordered[1],
        circumference - ordered[2] + ordered[0],
    )
    return tuple(sorted(result))


def measure(values, circumference, sectors):
    arc_gaps = gaps(values, circumference)
    positions = tuple(value * sectors // circumference for value in values)
    return {
        "min_arc_fraction": (circumference - arc_gaps[2]) / circumference,
        "max_gap_fraction": arc_gaps[2] / circumference,
        "sector_count": len(set(positions)),
        "all_same_sector": len(set(positions)) == 1,
        "gaps": arc_gaps,
        "positions": positions,
    }


def rotate(values, rotation, circumference):
    return tuple((value + rotation) % circumference for value in values)


def summarize(records, key):
    valid = [record[key] for record in records if record[key] is not None]
    if not valid:
        return {}
    return {
        "count": len(valid),
        "mean_min_arc_fraction": statistics.fmean(row["min_arc_fraction"] for row in valid),
        "median_min_arc_fraction": statistics.median(row["min_arc_fraction"] for row in valid),
        "minimum_min_arc_fraction": min(row["min_arc_fraction"] for row in valid),
        "all_same_sector_fraction": statistics.fmean(row["all_same_sector"] for row in valid),
        "mean_sector_count": statistics.fmean(row["sector_count"] for row in valid),
    }


def run(max_scalar, sectors, output_csv, output_json):
    if max_scalar < 2:
        raise ValueError("max scalar must be at least 2")
    circumference = max_scalar
    analyzer = DomeInvariantAnalyzer()
    base_values = triangle_values(BASE_POINT, analyzer)
    # Rotate in the field first. Projecting Lambda(G)=1 to a small search
    # circle would round every one-step rotation to zero.
    base_rotation_field = base_values[0]
    point = point_from_scalar(max_scalar)
    negative_g = EllipticCurvePoint(BASE_POINT.x, (-BASE_POINT.y) % FIELD_PRIME)
    records = []
    for scalar in range(max_scalar, 0, -1):
        values = triangle_values(point, analyzer)
        projected = tuple(value * circumference // FIELD_PRIME for value in values)
        raw = measure(projected, circumference, sectors)
        # One inverse-G step corresponds to one opposite base rotation. At k=0
        # this rotation is zero, so the reference is anchored at the endpoint.
        rotation_field = (-scalar * base_rotation_field) % FIELD_PRIME
        rotated_projected = tuple(
            ((value + rotation_field) % FIELD_PRIME) * circumference // FIELD_PRIME
            for value in values
        )
        rotated = measure(rotated_projected, circumference, sectors)
        records.append({
            "scalar": scalar,
            "offset_from_zero": scalar,
            "rotation_field": rotation_field,
            "raw": raw,
            "rotated": rotated,
            "wrap_integer": sum(values) // FIELD_PRIME,
        })
        point = point + negative_g
    if not point.is_infinity:
        raise RuntimeError("descending G trace did not reach the zero point")

    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["scalar", "rotation_field", "wrap_integer", "raw_min_arc", "rotated_min_arc",
                         "raw_sector_count", "rotated_sector_count", "raw_all_same", "rotated_all_same"])
        for record in records:
            writer.writerow([
                record["scalar"], record["rotation_field"], record["wrap_integer"],
                record["raw"]["min_arc_fraction"], record["rotated"]["min_arc_fraction"],
                record["raw"]["sector_count"], record["rotated"]["sector_count"],
                record["raw"]["all_same_sector"], record["rotated"]["all_same_sector"],
            ])

    report = {
        "max_scalar": max_scalar,
        "circle_circumference": circumference,
        "sectors": sectors,
        "base_lambda_orbit": base_values,
        "terminal_point": "0*G = infinity; Lambda undefined",
        "raw_summary": summarize(records, "raw"),
        "rotated_summary": summarize(records, "rotated"),
        "wrap_integer_counts": {
            str(value): sum(record["wrap_integer"] == value for record in records)
            for value in sorted({record["wrap_integer"] for record in records})
        },
        "first_records": records[:10],
        "last_records": records[-10:],
        "csv": str(output_csv),
    }
    output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-scalar", type=int, default=10000)
    parser.add_argument("--sectors", type=int, default=32)
    parser.add_argument("--csv", type=Path, default=Path("results/triangle_descending_g.csv"))
    parser.add_argument("--json", type=Path, default=Path("results/triangle_descending_g.json"))
    args = parser.parse_args()
    report = run(args.max_scalar, args.sectors, args.csv, args.json)
    print(f"trace: {args.max_scalar}G -> 1G -> 0G, circumference={report['circle_circumference']}")
    print(f"base Lambda orbit: {report['base_lambda_orbit']}")
    print(f"wrap integers: {report['wrap_integer_counts']}")
    print("mode       count  mean arc  median arc  minimum arc  same-sector  mean sectors")
    for name in ("raw_summary", "rotated_summary"):
        row = report[name]
        print(f"{name:9s} {row['count']:6d} {row['mean_min_arc_fraction']:8.3%} "
              f"{row['median_min_arc_fraction']:10.3%} {row['minimum_min_arc_fraction']:11.3%} "
              f"{row['all_same_sector_fraction']:11.3%} {row['mean_sector_count']:12.3f}")
    print(f"csv={args.csv}; report={args.json}")


if __name__ == "__main__":
    main()
