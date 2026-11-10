#!/usr/bin/env python3
"""Search for a compact deterministic decoder from angle marker to offset.

The test is deliberately falsifiable.  It measures whether a quantized
two-angle fingerprint has a simple affine decoder, whether sorting markers
clusters scalar offsets, and whether consecutive scalars have predictable
marker transitions.  A high bucket hit rate alone is not enough: the decoder
must reproduce the offset without storing a hidden lookup table.
"""

import argparse
import math
from concurrent.futures import ProcessPoolExecutor

from parallel_angle_table_recovery import angle_pair, fingerprint, point_from_scalar
from w1 import DomeInvariantAnalyzer


def worker(args):
    start, first, count, bits = args
    analyzer = DomeInvariantAnalyzer()
    return [
        (fingerprint(angle_pair(analyzer, point_from_scalar(start + offset)), bits), offset)
        for offset in range(first, first + count)
    ]


def build(start, size, bits, workers, chunk):
    jobs = [(start, first, min(chunk, size - first), bits)
            for first in range(0, size, chunk)]
    records = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(worker, jobs):
            records.extend(result)
    return records


def affine_fit(records):
    """Least-squares fit offset ~= a*marker+b after integer normalization."""
    n = len(records)
    xs = [marker for marker, _ in records]
    ys = [offset for _, offset in records]
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    denominator = sum((x - mean_x) ** 2 for x in xs)
    slope = (sum((x - mean_x) * (y - mean_y) for x, y in records) /
             denominator if denominator else 0.0)
    intercept = mean_y - slope * mean_x
    errors = [abs((slope * x + intercept) - y) for x, y in records]
    return slope, intercept, sum(errors) / n, max(errors)


def sorted_marker_locality(records, size):
    ordered = sorted(records)
    deltas = [abs(right[1] - left[1]) for left, right in zip(ordered, ordered[1:])]
    return {
        "mean_delta_fraction": sum(deltas) / len(deltas) / size if deltas else 0.0,
        "median_delta_fraction": sorted(deltas)[len(deltas) // 2] / size if deltas else 0.0,
        "max_delta_fraction": max(deltas) / size if deltas else 0.0,
    }


def transition_stats(records, size):
    ordered = sorted(records, key=lambda item: item[1])
    deltas = [abs(right[0] - left[0]) for left, right in zip(ordered, ordered[1:])]
    marker_bits = max(1, max((marker.bit_length() for marker, _ in records), default=1))
    return {
        "mean_transition_fraction": sum(deltas) / len(deltas) / (1 << marker_bits) if deltas else 0.0,
        "median_transition_fraction": sorted(deltas)[len(deltas) // 2] / (1 << marker_bits) if deltas else 0.0,
        "marker_bits": marker_bits,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-size", type=int, default=10000)
    parser.add_argument("--bits", type=int, nargs="+", default=[8, 12, 16])
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--chunk-size", type=int, default=1000)
    args = parser.parse_args()
    start = 1 << 134
    print("Marker decoder search")
    print(f"  complete window: {args.range_size} offsets at [2^134,2^134+W)")
    print("  bits  unique  collisions  affine MAE  affine max  sorted median delta  transition median")
    for bits in args.bits:
        records = build(start, args.range_size, bits, args.workers, args.chunk_size)
        unique = len({marker for marker, _ in records})
        slope, intercept, mae, maximum = affine_fit(records)
        locality = sorted_marker_locality(records, args.range_size)
        transitions = transition_stats(records, args.range_size)
        print(
            f"  {bits:4d}  {unique:7d}  {len(records)-unique:10d}"
            f"  {mae/args.range_size:10.3%}  {maximum/args.range_size:10.3%}"
            f"  {locality['median_delta_fraction']:18.3%}"
            f"  {transitions['median_transition_fraction']:16.3%}"
        )
        print(f"       affine slope={slope:.6e}, intercept={intercept:.3f}")
    print("  A useful marker decoder must make affine/locality errors small without storing offset lists.")


if __name__ == "__main__":
    main()
