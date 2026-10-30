#!/usr/bin/env python3
"""Compare offset storage with marker/rank encodings on a complete window.

A marker can identify a bucket, but recovery also needs the member's position
inside that bucket or a decoder from marker to scalar offset.  This benchmark
measures the exact bit cost of those alternatives on real fingerprints.
"""

import argparse
import math
from collections import defaultdict

from storage_scaling_benchmark import build


def bits(n):
    return max(1, (n - 1).bit_length())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-size", type=int, default=100000)
    parser.add_argument("--angle-bits", type=int, default=12)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--chunk-size", type=int, default=1000)
    args = parser.parse_args()
    records = build(1 << 134, args.range_size, args.angle_bits,
                    args.workers, args.chunk_size)
    buckets = defaultdict(list)
    for fingerprint, offset in records:
        buckets[fingerprint].append(offset)
    bucket_count = len(buckets)
    fingerprint_bits = 2 * args.angle_bits
    offset_bits = bits(args.range_size)
    bucket_id_bits = bits(bucket_count)
    rank_bits = max(bits(len(values)) for values in buckets.values())

    # Three representations, excluding a shared header/index:
    # full offset; bucket marker + rank; bucket marker only (not decodable if
    # buckets contain multiple members).
    full_bits = len(records) * offset_bits
    marker_rank_bits = len(records) * (bucket_id_bits + rank_bits)
    marker_only_bits = len(records) * bucket_id_bits
    print("Marker encoding benchmark")
    print(f"  records: {len(records)}")
    print(f"  buckets: {bucket_count}")
    print(f"  fingerprint bits: {fingerprint_bits}")
    print(f"  offset bits: {offset_bits}")
    print(f"  bucket-id bits: {bucket_id_bits}")
    print(f"  max rank bits: {rank_bits}")
    print(f"  full offset storage: {full_bits/8:.1f} bytes")
    print(f"  marker + rank storage: {marker_rank_bits/8:.1f} bytes")
    print(f"  marker-only storage: {marker_only_bits/8:.1f} bytes")
    print("  marker-only decodable: yes only when every bucket has one member")
    print("  rule: a marker is useful only if it has a deterministic decoder to offset")


if __name__ == "__main__":
    main()
