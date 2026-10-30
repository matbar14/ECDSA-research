#!/usr/bin/env python3
"""Test whether quantized angles can recover a block anchor and offset.

The proposed scheme is k = anchor + offset.  For it to save scalar storage,
an angle cell must be local in scalar order: entries in one cell should lie in
one (or a few) anchor blocks.  This script measures that locality and also
tests a two-angle key (Lambda(P), Lambda(phi(P))).
"""

import argparse
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


def lambda_pair(analyzer, point):
    first = analyzer.compute_torus_invariant_lambda(point)
    second_point = EllipticCurvePoint(
        ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y
    )
    second = analyzer.compute_torus_invariant_lambda(second_point)
    return first, second


def quantize(value, bits):
    return value * (1 << bits) // FIELD_PRIME


def measure(cells, range_size, block_size):
    block_counts = [len({offset // block_size for offset in offsets})
                    for offsets in cells.values()]
    spans = [max(offsets) - min(offsets) + 1 for offsets in cells.values()]
    unique_block_cells = sum(count == 1 for count in block_counts)
    unique_offsets = sum(len(offsets) == 1 for offsets in cells.values())
    return {
        "cells": len(cells),
        "mean_candidates": sum(len(offsets) for offsets in cells.values()) / range_size,
        "unique_block_fraction": unique_block_cells / len(cells),
        "mean_blocks": sum(block_counts) / len(block_counts),
        "max_blocks": max(block_counts),
        "mean_span_fraction": sum(spans) / len(spans) / range_size,
        "max_span_fraction": max(spans) / range_size,
        "unique_offset_fraction": unique_offsets / range_size,
    }


def run(range_start, range_size, block_size, bits):
    if range_start < 1 or range_size <= 0 or range_start + range_size > GROUP_ORDER:
        raise ValueError("range must be contained in [1, n)")
    if block_size <= 0 or not bits or any(bit <= 0 for bit in bits):
        raise ValueError("block_size and bits must be positive")

    analyzer = DomeInvariantAnalyzer()
    first_cells = {b: defaultdict(list) for b in bits}
    pair_cells = {b: defaultdict(list) for b in bits}
    for offset in range(range_size):
        first, second = lambda_pair(analyzer, point_from_scalar(range_start + offset))
        for b in bits:
            first_cells[b][quantize(first, b)].append(offset)
            pair_cells[b][(quantize(first, b), quantize(second, b))].append(offset)

    return {
        "range_start": range_start,
        "range_size": range_size,
        "block_size": block_size,
        "results": {
            b: {"one_angle": measure(first_cells[b], range_size, block_size),
                "two_angles": measure(pair_cells[b], range_size, block_size)}
            for b in bits
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-size", type=int, default=10000)
    parser.add_argument("--block-size", type=int, default=100)
    parser.add_argument("--bits", type=int, nargs="+", default=[8, 12, 16, 20])
    args = parser.parse_args()
    report = run(1 << 255, args.range_size, args.block_size, args.bits)
    print("Angle anchor locality experiment")
    print(f"  range size: {report['range_size']}")
    print(f"  anchor block size: {report['block_size']}")
    print("  bits  key          cells  mean blocks  max blocks  mean span  unique offset")
    for bits in args.bits:
        for key_name, label in (("one_angle", "one angle"), ("two_angles", "two angles")):
            result = report["results"][bits][key_name]
            print(
                f"  {bits:4d}  {label:10s}  {result['cells']:5d}"
                f"  {result['mean_blocks']:11.3f}  {result['max_blocks']:10d}"
                f"  {result['mean_span_fraction']:9.3%}"
                f"  {result['unique_offset_fraction']:13.3%}"
            )


if __name__ == "__main__":
    main()
