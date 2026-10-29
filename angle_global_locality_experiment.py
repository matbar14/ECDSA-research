#!/usr/bin/env python3
"""Measure anchor locality on a wide, non-contiguous scalar sample.

Unlike angle_anchor_experiment.py, this script samples k uniformly across a
wide interval. For every repeated quantized angle cell it measures the scalar
span and the number of offset bits needed to address that span.
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


def values(analyzer, point):
    first = analyzer.compute_torus_invariant_lambda(point)
    second_point = EllipticCurvePoint(
        ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y
    )
    second = analyzer.compute_torus_invariant_lambda(second_point)
    return first, second


def quantize(value, bits):
    return value * (1 << bits) // FIELD_PRIME


def summarize(cells, width):
    repeated = [offsets for offsets in cells.values() if len(offsets) > 1]
    if not repeated:
        return {
            "cells": len(cells), "repeated": 0, "mean_span_fraction": 0.0,
            "max_span_fraction": 0.0, "mean_offset_bits": 0.0,
            "max_offset_bits": 0, "min_span_fraction": 0.0,
            "weighted_good_1pct": 0.0, "weighted_good_5pct": 0.0,
            "weighted_good_10pct": 0.0,
        }
    spans = [max(x) - min(x) + 1 for x in repeated]
    bits = [max(1, (span - 1).bit_length()) for span in spans]
    weights = [len(offsets) for offsets in repeated]
    return {
        "cells": len(cells),
        "repeated": len(repeated),
        "mean_span_fraction": sum(spans) / len(spans) / width,
        "max_span_fraction": max(spans) / width,
        "min_span_fraction": min(spans) / width,
        "mean_offset_bits": sum(bits) / len(bits),
        "max_offset_bits": max(bits),
        "weighted_good_1pct": sum(w for w, span in zip(weights, spans)
                                   if span <= width * 0.01) / sum(weights),
        "weighted_good_5pct": sum(w for w, span in zip(weights, spans)
                                   if span <= width * 0.05) / sum(weights),
        "weighted_good_10pct": sum(w for w, span in zip(weights, spans)
                                    if span <= width * 0.10) / sum(weights),
    }


def run(sample_count, width, bits, seed):
    if sample_count <= 0 or width <= 0:
        raise ValueError("sample_count and width must be positive")
    if width >= GROUP_ORDER:
        raise ValueError("width must be smaller than group order")
    rng = random.Random(seed)
    start = 1 << 255
    analyzer = DomeInvariantAnalyzer()
    cells = {b: {"one_angle": defaultdict(list),
                 "two_angles": defaultdict(list)} for b in bits}

    for _ in range(sample_count):
        offset = rng.randrange(width)
        first, second = values(analyzer, point_from_scalar(start + offset))
        for b in bits:
            cells[b]["one_angle"][quantize(first, b)].append(offset)
            cells[b]["two_angles"][(quantize(first, b), quantize(second, b))].append(offset)

    return {
        b: {key: summarize(value, width) for key, value in groups.items()}
        for b, groups in cells.items()
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--width-bits", type=int, default=135)
    parser.add_argument("--bits", type=int, nargs="+", default=[8, 12, 16, 20])
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    report = run(args.samples, 1 << args.width_bits, args.bits, args.seed)
    print("Global angle-anchor locality experiment")
    print(f"  samples: {args.samples}")
    print(f"  sampled width: 2^{args.width_bits}")
    print("  bits  key          cells  repeated  min span   mean span   max span   good<=1%  good<=5%  good<=10%  mean offset bits  max bits")
    for b in args.bits:
        for key in ("one_angle", "two_angles"):
            r = report[b][key]
            print(
                f"  {b:4d}  {key:10s}  {r['cells']:6d}  {r['repeated']:8d}"
                f"  {r['min_span_fraction']:9.3%}  {r['mean_span_fraction']:9.3%}"
                f"  {r['max_span_fraction']:9.3%}"
                f"  {r['weighted_good_1pct']:8.3%}  {r['weighted_good_5pct']:8.3%}"
                f"  {r['weighted_good_10pct']:9.3%}"
                f"  {r['mean_offset_bits']:16.2f}  {r['max_offset_bits']:8d}"
            )


if __name__ == "__main__":
    main()
