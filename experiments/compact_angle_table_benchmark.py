#!/usr/bin/env python3
"""Measure bit-packed storage for the working two-angle table scheme.

The lookup semantics are unchanged:
    public P -> exact (Lambda0,Lambda1) -> quantized fingerprint -> offsets
    -> exact verification.

This benchmark packs each (fingerprint, offset) pair into the minimum fixed
number of bits for the selected finite window and compares it with the current
16-byte record format. It is a storage benchmark, not a full 2^135 build.
"""

import argparse
import math
import os
import struct


def bits_required(value_count):
    return max(1, (value_count - 1).bit_length())


def pack_records(records, fingerprint_bits, offset_bits, output):
    records = list(records)
    record_bits = fingerprint_bits + offset_bits
    total_bits = len(records) * record_bits
    packed = bytearray((total_bits + 7) // 8)
    cursor = 0
    mask = (1 << record_bits) - 1
    for fingerprint, offset in records:
        value = ((fingerprint << offset_bits) | offset) & mask
        for bit in range(record_bits):
            if value & (1 << bit):
                packed[(cursor + bit) >> 3] |= 1 << ((cursor + bit) & 7)
        cursor += record_bits
    with open(output, "wb") as handle:
        handle.write(struct.pack("<QQQQ", len(records), fingerprint_bits,
                                 offset_bits, record_bits))
        handle.write(packed)
    return os.path.getsize(output), record_bits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-size", type=int, default=100000)
    parser.add_argument("--angle-bits", type=int, default=12)
    parser.add_argument("--records", type=int, default=None)
    parser.add_argument("--output", default="results/compact_angle_table.bin")
    args = parser.parse_args()
    if args.range_size <= 0 or not 1 <= args.angle_bits <= 32:
        parser.error("invalid range-size or angle-bits")
    records_count = args.records or args.range_size
    fingerprint_bits = 2 * args.angle_bits
    offset_bits = bits_required(args.range_size)
    compact_bytes, record_bits = pack_records(
        ((index % (1 << fingerprint_bits), index) for index in range(records_count)),
        fingerprint_bits, offset_bits, args.output
    )
    old_bytes = 24 + records_count * 16
    theoretical = (record_bits * records_count + 7) // 8 + 32
    print("Compact angle-table storage benchmark")
    print(f"  records: {records_count}")
    print(f"  fingerprint bits: {fingerprint_bits}")
    print(f"  offset bits: {offset_bits}")
    print(f"  packed record bits: {record_bits}")
    print(f"  packed file bytes: {compact_bytes}")
    print(f"  current 16-byte-record estimate: {old_bytes}")
    print(f"  storage reduction: {old_bytes / compact_bytes:.3f}x")
    print(f"  output: {args.output}")


if __name__ == "__main__":
    main()
