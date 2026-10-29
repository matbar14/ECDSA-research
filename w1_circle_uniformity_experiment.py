#!/usr/bin/env python3
"""Measure how residues of w1.py invariants fill the unit circle.

The map a -> (cos(2*pi*a/p), sin(2*pi*a/p)) is a chosen visualization
convention for canonical residues a in [0,p). Sector counts use integer
arithmetic, so floating point rounding affects only the figure.
"""

import argparse
import json
import math
import random
from collections import Counter
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec
from scipy.stats import chi2

from w1 import (
    BASE_POINT, DomeInvariantAnalyzer, EllipticCurvePoint,
    ENDOMORPHISM_BETA, FIELD_PRIME, GROUP_ORDER,
)


NAMES = ("x3", "radius", "area", "curvature", "lambda_0", "lambda_1", "lambda_2")


def point_from_scalar(k):
    q = ec.derive_private_key(k, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(q.x, q.y)


def residues(point, analyzer):
    x3, y = analyzer.torus_coordinates(point)
    orbit = [point]
    for _ in range(2):
        previous = orbit[-1]
        orbit.append(EllipticCurvePoint(ENDOMORPHISM_BETA * previous.x % FIELD_PRIME, previous.y))
    lam = [analyzer.compute_torus_invariant_lambda(p) for p in orbit]
    return (
        x3, (x3 + y) % FIELD_PRIME, x3 * y % FIELD_PRIME,
        (x3 - y) % FIELD_PRIME, *lam,
    )


def sample(start, count, mode, seed, width):
    if start < 1 or count < 1 or start + width > GROUP_ORDER or count > width:
        raise ValueError("sample must be nonempty and lie within [1, group order)")
    analyzer = DomeInvariantAnalyzer()
    if mode == "consecutive":
        point = point_from_scalar(start)
        for _ in range(count):
            yield residues(point, analyzer)
            point = point + BASE_POINT
    else:
        rng = random.Random(seed)
        offsets = set()
        while len(offsets) < count:
            offsets.add(rng.randrange(width))
        for offset in offsets:
            yield residues(point_from_scalar(start + offset), analyzer)


def analyze(data, sectors):
    expected = len(data) / sectors
    result = {}
    for column, name in enumerate(NAMES):
        values = [row[column] for row in data]
        counts = Counter(value * sectors // FIELD_PRIME for value in values)
        bins = [counts[i] for i in range(sectors)]
        statistic = sum((v - expected) ** 2 / expected for v in bins)
        result[name] = {
            "occupied_sectors": len(counts),
            "min_count": min(bins),
            "max_count": max(bins),
            "chi_square": statistic,
            "degrees_of_freedom": sectors - 1,
            "chi_square_p_value": chi2.sf(statistic, sectors - 1),
            "counts": bins,
        }
    return result


def draw(data, result, sectors, mode, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(NAMES), 2, figsize=(13, 23), constrained_layout=True)
    theta_grid = [math.tau * i / 512 for i in range(513)]
    outline_x = [math.cos(t) for t in theta_grid]
    outline_y = [math.sin(t) for t in theta_grid]
    for column, name in enumerate(NAMES):
        values = [row[column] for row in data]
        angles = [math.tau * (v / FIELD_PRIME) for v in values]
        circle, histogram = axes[column]
        circle.plot(outline_x, outline_y, color="#26364f", linewidth=1)
        circle.scatter([math.cos(t) for t in angles], [math.sin(t) for t in angles],
                       s=3, alpha=0.08, color="#1680a8", rasterized=True)
        circle.set_aspect("equal")
        circle.set_xlim(-1.08, 1.08)
        circle.set_ylim(-1.08, 1.08)
        circle.set_title(f"{name}: unit circle")
        circle.axis("off")
        histogram.bar(range(sectors), result[name]["counts"], width=0.95, color="#1680a8")
        histogram.axhline(len(data) / sectors, color="#c43b34", linestyle="--", linewidth=1)
        histogram.set_title(f"{name}: p={result[name]['chi_square_p_value']:.3g}")
        histogram.set_xlabel("equal angular sector")
        histogram.set_ylabel("count")
    fig.suptitle(f"w1.py invariant residues on the unit circle: {mode}, n={len(data)}", fontsize=16)
    fig.savefig(output, dpi=170)
    plt.close(fig)


def run(start, count, width, sectors, mode, seed, figure, json_output):
    if sectors < 2 or count / sectors < 5:
        raise ValueError("need at least 5 expected observations per sector")
    data = list(sample(start, count, mode, seed, width))
    result = analyze(data, sectors)
    report = {
        "start": start, "sample_count": count, "sampling_width": width,
        "mode": mode, "seed": seed, "sectors": sectors,
        "circle_map": "a -> (cos(2*pi*a/p), sin(2*pi*a/p))",
        "invariants": result,
    }
    draw(data, result, sectors, mode, figure)
    json_output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda value: int(value, 0), default=1 << 134)
    parser.add_argument("--count", type=int, default=10000)
    parser.add_argument("--width", type=lambda value: int(value, 0), default=1 << 135)
    parser.add_argument("--sectors", type=int, default=32)
    parser.add_argument("--mode", choices=("consecutive", "random"), default="consecutive")
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--figure", type=Path, default=Path("w1_circle_uniformity.png"))
    parser.add_argument("--json", type=Path, default=Path("w1_circle_uniformity.json"))
    args = parser.parse_args()
    report = run(args.start, args.count, args.width, args.sectors,
                 args.mode, args.seed, args.figure, args.json)
    print(f"mode={args.mode}, samples={args.count}, sectors={args.sectors}, expected={args.count / args.sectors:.1f}/sector")
    print("invariant    occupied   min..max    chi2     p-value")
    for name, row in report["invariants"].items():
        print(f"{name:11s} {row['occupied_sectors']:2d}/{args.sectors:<2d} "
              f"{row['min_count']:4d}..{row['max_count']:<4d} "
              f"{row['chi_square']:8.2f} {row['chi_square_p_value']:.6g}")
    print(f"figure={args.figure}, report={args.json}")


if __name__ == "__main__":
    main()
