#!/usr/bin/env python3
"""Inspect the Lambda trajectory on the phase torus.

For consecutive scalars k this records
    (Lambda(kG)/p, Lambda(phi(kG))/p)
and writes an SVG scatter/trajectory plot.  It also measures quantized torus
cells as possible anchors: repeated cells, scalar spans, and short-span hit
fractions.  All collision checks use exact field residues; floating point is
used only for plotting and circular distances.
"""

import argparse
import math
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


def point_from_scalar(k):
    q = ec.derive_private_key(k, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(q.x, q.y)


def torus_values(analyzer, point):
    first = analyzer.compute_torus_invariant_lambda(point)
    second_point = EllipticCurvePoint(
        ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y
    )
    second = analyzer.compute_torus_invariant_lambda(second_point)
    return first, second


def quantize(values, bits):
    return tuple(value * (1 << bits) // FIELD_PRIME for value in values)


def circular_delta(left, right):
    delta = abs(left - right)
    return min(delta, 1.0 - delta)


def write_svg(points, output, title):
    width = 900
    height = 900
    margin = 60
    plot = width - 2 * margin
    circles = []
    for index, (u, v) in enumerate(points):
        x = margin + plot * u
        y = height - margin - plot * v
        color = "#d62728" if index == 0 else "#1f77b4"
        circles.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="1.7" fill="{color}"/>')
    polyline = " ".join(
        f"{margin + plot*u:.3f},{height-margin-plot*v:.3f}" for u, v in points
    )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{margin}" y="28" font-family="sans-serif" font-size="18">{title}</text>
<rect x="{margin}" y="{margin}" width="{plot}" height="{plot}" fill="#fafafa" stroke="#333"/>
<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{margin}" stroke="#ddd"/>
<polyline points="{polyline}" fill="none" stroke="#9ecae1" stroke-width="0.6" opacity="0.7"/>
{''.join(circles)}
<text x="{margin}" y="{height-18}" font-family="sans-serif" font-size="13">Lambda(kG) / p</text>
<text x="15" y="{margin+plot/2}" transform="rotate(-90 15 {margin+plot/2})" font-family="sans-serif" font-size="13">Lambda(phi(kG)) / p</text>
</svg>'''
    with open(output, "w", encoding="utf-8") as handle:
        handle.write(svg)


def run(range_start, range_size, bits, output, sample_mode, seed):
    if range_start < 1 or range_size <= 0 or range_start + range_size > GROUP_ORDER:
        raise ValueError("range must be contained in [1,n)")
    analyzer = DomeInvariantAnalyzer()
    scalars = list(range(range_start, range_start + range_size))
    if sample_mode == "random":
        rng = random.Random(seed)
        scalars = sorted(rng.sample(scalars, min(range_size, 10000)))

    exact = []
    plotted = []
    for scalar in scalars:
        values = torus_values(analyzer, point_from_scalar(scalar))
        exact.append((scalar, values))
        plotted.append((values[0] / FIELD_PRIME, values[1] / FIELD_PRIME))

    cells = defaultdict(list)
    for scalar, values in exact:
        cells[quantize(values, bits)].append(scalar)
    repeated = [members for members in cells.values() if len(members) > 1]
    spans = [max(members) - min(members) + 1 for members in repeated]
    weights = [len(members) for members in repeated]
    total_repeated = sum(weights)

    adjacent_steps = [
        math.hypot(circular_delta(left[0], right[0]), circular_delta(left[1], right[1]))
        for left, right in zip(plotted, plotted[1:])
    ]
    write_svg(plotted, output, f"Lambda phase-torus trajectory ({sample_mode}, {len(plotted)} points)")

    return {
        "points": len(plotted),
        "cells": len(cells),
        "repeated_cells": len(repeated),
        "mean_candidates": sum(len(members) for members in cells.values()) / len(cells),
        "min_span_fraction": (min(spans) / range_size if spans else 0.0),
        "mean_span_fraction": (sum(spans) / len(spans) / range_size if spans else 0.0),
        "max_span_fraction": (max(spans) / range_size if spans else 0.0),
        "good_1pct": (sum(weight for weight, span in zip(weights, spans) if span <= range_size * .01) / total_repeated if repeated else 0.0),
        "good_5pct": (sum(weight for weight, span in zip(weights, spans) if span <= range_size * .05) / total_repeated if repeated else 0.0),
        "step_mean": sum(adjacent_steps) / len(adjacent_steps) if adjacent_steps else 0.0,
        "step_min": min(adjacent_steps) if adjacent_steps else 0.0,
        "step_max": max(adjacent_steps) if adjacent_steps else 0.0,
        "svg": output,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-start", type=int, default=1 << 134)
    parser.add_argument("--range-size", type=int, default=2000)
    parser.add_argument("--bits", type=int, default=12)
    parser.add_argument("--mode", choices=("consecutive", "random"), default="consecutive")
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--output", default="figures/torus_trajectory.svg")
    args = parser.parse_args()
    report = run(args.range_start, args.range_size, args.bits, args.output, args.mode, args.seed)
    print("Lambda phase-torus trajectory")
    print(f"  mode: {args.mode}")
    print(f"  points: {report['points']}")
    print(f"  quantization bits/coordinate: {args.bits}")
    print(f"  occupied cells: {report['cells']}")
    print(f"  repeated cells: {report['repeated_cells']}")
    print(f"  mean candidates/cell: {report['mean_candidates']:.6f}")
    print(f"  repeated-cell span min/mean/max: {report['min_span_fraction']:.3%} / {report['mean_span_fraction']:.3%} / {report['max_span_fraction']:.3%}")
    print(f"  repeated-hit fraction span<=1%: {report['good_1pct']:.3%}")
    print(f"  repeated-hit fraction span<=5%: {report['good_5pct']:.3%}")
    print(f"  consecutive torus step mean/min/max: {report['step_mean']:.6f} / {report['step_min']:.6f} / {report['step_max']:.6f}")
    print(f"  svg: {report['svg']}")


if __name__ == "__main__":
    main()
