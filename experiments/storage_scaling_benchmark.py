#!/usr/bin/env python3
"""Benchmark structural compression of the working angle table.

It builds a complete finite window, groups records by quantized two-angle
fingerprint, and compares fixed-width records with simple bucket/delta/varint
encodings.  It also projects the measured bits/record to the full canonical
135-bit window [2^134,2^135).
"""

import argparse
import math
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

from parallel_angle_table_recovery import fingerprint, point_from_scalar, angle_pair
from w1 import DomeInvariantAnalyzer, GROUP_ORDER


def worker(args):
    start, first, count, bits = args
    analyzer = DomeInvariantAnalyzer()
    return [
        (fingerprint(angle_pair(analyzer, point_from_scalar(start + offset)), bits), offset)
        for offset in range(first, first + count)
    ]


def varint_bytes(value):
    total = 1
    while value >= 128:
        value >>= 7
        total += 1
    return total


def build(start, size, bits, workers, chunk):
    jobs = [(start, first, min(chunk, size - first), bits)
            for first in range(0, size, chunk)]
    records = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(worker, jobs):
            records.extend(result)
    return records


def summarize(records, range_size, bits):
    offset_bits = max(1, (range_size - 1).bit_length())
    grouped = defaultdict(list)
    for key, offset in records:
        grouped[key].append(offset)
    for offsets in grouped.values():
        offsets.sort()

    fixed_bits = 2 * bits + offset_bits
    fixed_bytes = (fixed_bits * len(records) + 7) // 8

    # Bucket format: fingerprint once, count as a varint, then sorted offsets
    # as first-offset plus positive deltas, each varint.
    bucket_bits = 0
    bucket_bytes = 0
    delta_bytes = 0
    for key, offsets in grouped.items():
        bucket_bytes += (2 * bits + 7) // 8 + varint_bytes(len(offsets))
        previous = 0
        for offset in offsets:
            delta = offset - previous
            delta_bytes += varint_bytes(delta)
            previous = offset
    bucket_bytes += delta_bytes

    return {
        "records": len(records),
        "buckets": len(grouped),
        "mean_bucket": len(records) / len(grouped),
        "max_bucket": max(map(len, grouped.values())),
        "fixed_bits": fixed_bits,
        "fixed_bytes": fixed_bytes,
        "bucket_bytes": bucket_bytes,
        "bucket_bits_per_record": bucket_bytes * 8 / len(records),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-size", type=int, default=100000)
    parser.add_argument("--bits", type=int, nargs="+", default=[8, 12, 16])
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--chunk-size", type=int, default=1000)
    args = parser.parse_args()
    start = 1 << 134
    if start + args.range_size > GROUP_ORDER:
        parser.error("range exceeds group order")

    print("Storage scaling benchmark")
    print(f"  complete finite window: {args.range_size} entries at [2^134,2^134+W)")
    print("  bits  buckets  mean bucket  max  fixed B/row  bucket B/row  projected full-window ZiB")
    full_classes = (1 << 134) / 6
    for bits in args.bits:
        records = build(start, args.range_size, bits, args.workers, args.chunk_size)
        result = summarize(records, args.range_size, bits)
        projected = full_classes * result["bucket_bits_per_record"] / 8 / 2**70
        print(
            f"  {bits:4d}  {result['buckets']:7d}  {result['mean_bucket']:11.3f}"
            f"  {result['max_bucket']:3d}  {result['fixed_bits']/8:11.3f}"
            f"  {result['bucket_bits_per_record']/8:12.3f}  {projected:18.3e}"
        )


if __name__ == "__main__":
    main()
