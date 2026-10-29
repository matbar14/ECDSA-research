#!/usr/bin/env python3
"""Build the working angle table with bit-packed real records.

This uses the same two-angle fingerprint and exact deferred verification as
parallel_angle_table_recovery.py, but packs each fingerprint/offset pair into
the minimum fixed-width bit record for the selected finite window.
"""

import argparse
import bisect
import os
import struct
from concurrent.futures import ProcessPoolExecutor

from parallel_angle_table_recovery import (
    angle_pair,
    build_table,
    exact_pair,
    fingerprint,
    point_from_scalar,
)
from w1 import GROUP_ORDER, DomeInvariantAnalyzer


HEADER = struct.Struct("<8s32sQQQQ")
MAGIC = b"CATBLv1\0"


def bits_required(count):
    return max(1, (count - 1).bit_length())


def pack_records(records, fingerprint_bits, offset_bits):
    width = fingerprint_bits + offset_bits
    total_bits = width * len(records)
    payload = bytearray((total_bits + 7) // 8)
    cursor = 0
    mask = (1 << width) - 1
    for key, offset in records:
        value = ((key << offset_bits) | offset) & mask
        for bit in range(width):
            if value & (1 << bit):
                payload[(cursor + bit) >> 3] |= 1 << ((cursor + bit) & 7)
        cursor += width
    return bytes(payload), width


def unpack_records(payload, count, fingerprint_bits, offset_bits, width):
    records = []
    mask = (1 << width) - 1
    for index in range(count):
        cursor = index * width
        value = 0
        for bit in range(width):
            if payload[(cursor + bit) >> 3] & (1 << ((cursor + bit) & 7)):
                value |= 1 << bit
        value &= mask
        records.append((value >> offset_bits, value & ((1 << offset_bits) - 1)))
    return records


def write_table(path, range_start, range_size, bits, records):
    payload, width = pack_records(records, 2 * bits, bits_required(range_size))
    with open(path, "wb") as handle:
        handle.write(HEADER.pack(
            MAGIC, range_start.to_bytes(32, "big"), range_size, bits,
            len(records), width
        ))
        handle.write(payload)
    return len(payload) + HEADER.size, width


def read_table(path):
    with open(path, "rb") as handle:
        header = handle.read(HEADER.size)
        magic, start_bytes, range_size, bits, count, width = HEADER.unpack(header)
        if magic != MAGIC:
            raise ValueError("not a compact angle table")
        payload = handle.read()
    records = unpack_records(payload, count, 2 * bits, bits_required(range_size), width)
    return int.from_bytes(start_bytes, "big"), range_size, bits, records


def lookup(records, key):
    first = bisect.bisect_left(records, (key, -1))
    result = []
    for stored, offset in records[first:]:
        if stored != key:
            break
        result.append(offset)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bit-length", type=int, default=135)
    parser.add_argument("--range-size", type=int, default=100000)
    parser.add_argument("--bits", type=int, default=12)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--queries", type=int, default=100)
    parser.add_argument("--output", default="compact_angle_table.bin")
    args = parser.parse_args()
    range_start = 1 << (args.bit_length - 1)
    if range_start + args.range_size > GROUP_ORDER:
        parser.error("range must be contained in [1,n)")
    print("Compact parallel angle table")
    records, build_seconds, _ = build_table(
        range_start, args.range_size, args.bits, args.workers,
        args.chunk_size, args.output + ".temporary"
    )
    records.sort()
    file_bytes, width = write_table(
        args.output, range_start, args.range_size, args.bits, records
    )
    loaded_start, loaded_size, loaded_bits, loaded = read_table(args.output)
    analyzer = DomeInvariantAnalyzer()
    correct = 0
    bucket_sizes = []
    for offset in range(min(args.queries, args.range_size)):
        scalar = range_start + offset
        values = exact_pair(analyzer, scalar)
        bucket = lookup(loaded, fingerprint(values, loaded_bits))
        bucket_sizes.append(len(bucket))
        matches = [candidate for candidate in bucket
                   if exact_pair(analyzer, loaded_start + candidate) == values]
        correct += matches == [offset]
    print(f"  range start: 2^{args.bit_length - 1}")
    print(f"  records: {len(records)}")
    print(f"  packed record bits: {width}")
    print(f"  packed file bytes: {file_bytes}")
    print(f"  old 16-byte estimate: {args.range_size * 16 + 72}")
    print(f"  storage reduction: {(args.range_size * 16 + 72) / file_bytes:.3f}x")
    print(f"  build seconds: {build_seconds:.2f}")
    print(f"  mean bucket: {sum(bucket_sizes)/len(bucket_sizes):.3f}")
    print(f"  exact recoveries: {correct}/{len(bucket_sizes)}")
    print(f"  output: {args.output}")


if __name__ == "__main__":
    main()
