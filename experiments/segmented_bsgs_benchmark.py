#!/usr/bin/env python3
"""Benchmark segmented BSGS storage with local offsets.

The range is split into segments.  Each segment stores quantized two-angle
fingerprints and offsets local to that segment; a query tests the segment
tables and verifies exact Lambda values.  This measures whether segmentation
turns the global offset into a genuinely cheaper local field.
"""

import argparse
import math
import time
from concurrent.futures import ProcessPoolExecutor

from parallel_angle_table_recovery import angle_pair, fingerprint, point_from_scalar
from w1 import DomeInvariantAnalyzer, GROUP_ORDER


def worker(args):
    range_start, segment_start, count, bits = args
    analyzer = DomeInvariantAnalyzer()
    return [
        (fingerprint(angle_pair(analyzer, point_from_scalar(range_start + segment_start + local)), bits), local)
        for local in range(count)
    ]


def bits_required(count):
    return max(1, (count - 1).bit_length())


def run(range_start, range_size, segment_size, bits, workers, chunk_size):
    if segment_size <= 0 or range_size <= 0:
        raise ValueError("segment_size and range_size must be positive")
    segment_count = (range_size + segment_size - 1) // segment_size
    jobs = []
    for segment in range(segment_count):
        start = segment * segment_size
        count = min(segment_size, range_size - start)
        for local in range(0, count, chunk_size):
            jobs.append((range_start, start + local, min(chunk_size, count - local), bits))

    started = time.perf_counter()
    segment_records = [[] for _ in range(segment_count)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for job, result in zip(jobs, pool.map(worker, jobs)):
            absolute_segment_start = job[1]
            segment = absolute_segment_start // segment_size
            local_base = absolute_segment_start % segment_size
            segment_records[segment].extend((key, local_base + offset) for key, offset in result)

    build_seconds = time.perf_counter() - started
    fingerprint_bits = 2 * bits
    total_bits = 0
    bucket_sizes = []
    for records in segment_records:
        records.sort()
        segment_width = min(segment_size, range_size)
        local_bits = bits_required(segment_width)
        total_bits += len(records) * (fingerprint_bits + local_bits)
        bucket_sizes.extend(len(records) for _ in records[:0])

    # Add one segment offset/length header per segment.  This is the cost of
    # using local offsets instead of a global offset.
    header_bits = segment_count * (bits_required(range_size) + bits_required(segment_size))
    packed_bytes = (total_bits + header_bits + 7) // 8
    global_bits = range_size * (fingerprint_bits + bits_required(range_size))
    return {
        "segments": segment_count,
        "records": sum(len(records) for records in segment_records),
        "local_bits": bits_required(segment_size),
        "global_bits_per_record": fingerprint_bits + bits_required(range_size),
        "local_bits_per_record": fingerprint_bits + bits_required(segment_size),
        "header_bits": header_bits,
        "packed_bytes": packed_bytes,
        "global_bytes": (global_bits + 7) // 8,
        "build_seconds": build_seconds,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-size", type=int, default=100000)
    parser.add_argument("--segment-size", type=int, nargs="+", default=[100, 1000, 10000])
    parser.add_argument("--bits", type=int, default=12)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--chunk-size", type=int, default=1000)
    args = parser.parse_args()
    start = 1 << 134
    print("Segmented BSGS/local-offset benchmark")
    print("  segment  count  local bits/record  packed bytes  global bytes  reduction")
    for segment in args.segment_size:
        report = run(start, args.range_size, segment, args.bits,
                     args.workers, args.chunk_size)
        print(
            f"  {segment:7d}  {report['segments']:5d}"
            f"  {report['local_bits_per_record']:17d}"
            f"  {report['packed_bytes']:12d}"
            f"  {report['global_bytes']:12d}"
            f"  {report['global_bytes']/report['packed_bytes']:.3f}x"
        )
        print(f"      build seconds: {report['build_seconds']:.2f}")


if __name__ == "__main__":
    main()
