#!/usr/bin/env python3
"""Full bounded-range build/query test for a quantized angle table.

The table covers every scalar in a finite test range.  Workers store only a
packed quantized two-angle fingerprint and the range offset.  A query first
looks up the fingerprint and then recomputes exact Lambda values for the small
candidate bucket, modeling deferred verification without storing k in the
fingerprint itself.
"""

import argparse
import bisect
import os
import random
import struct
import time
from concurrent.futures import ProcessPoolExecutor

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import (
    DomeInvariantAnalyzer,
    EllipticCurvePoint,
    ENDOMORPHISM_BETA,
    FIELD_PRIME,
    GROUP_ORDER,
)


MAGIC = b"ATBLv2\0\0"
HEADER = struct.Struct("<8s32sQQQQ")
RECORD = struct.Struct("<QQ")  # quantized fingerprint, scalar offset


def point_from_scalar(k):
    q = ec.derive_private_key(k, ec.SECP256K1()).public_key().public_numbers()
    return EllipticCurvePoint(q.x, q.y)


def compressed_public_key(point):
    prefix = 0x02 | (point.y & 1)
    return f"{prefix:02x}{point.x:064x}"


def point_from_public_key(encoded):
    try:
        raw = bytes.fromhex(encoded)
    except ValueError as error:
        raise ValueError(
            "target public key must be a full hex string "
            "(compressed: 66 hex characters, starting with 02 or 03)"
        ) from error
    if len(raw) not in (33, 65) or raw[0] not in (2, 3, 4):
        raise ValueError(
            "target public key must be SEC1 compressed (33 bytes) or "
            "uncompressed (65 bytes) hex"
        )
    public_key = ec.EllipticCurvePublicKey.from_encoded_point(ec.SECP256K1(), raw)
    q = public_key.public_numbers()
    return EllipticCurvePoint(q.x, q.y)


def angle_pair(analyzer, point):
    first = analyzer.compute_torus_invariant_lambda(point)
    second_point = EllipticCurvePoint(
        ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y
    )
    second = analyzer.compute_torus_invariant_lambda(second_point)
    return first, second


def fingerprint(values, bits):
    q1 = values[0] * (1 << bits) // FIELD_PRIME
    q2 = values[1] * (1 << bits) // FIELD_PRIME
    return (q1 << bits) | q2


def worker(args):
    range_start, offset_start, count, bits = args
    analyzer = DomeInvariantAnalyzer()
    output = []
    for offset in range(offset_start, offset_start + count):
        point = point_from_scalar(range_start + offset)
        output.append((fingerprint(angle_pair(analyzer, point), bits), offset))
    return output


def exact_pair(analyzer, scalar):
    return angle_pair(analyzer, point_from_scalar(scalar))


def build_table(range_start, range_size, bits, workers, chunk_size, output_path):
    jobs = [
        (range_start, offset, min(chunk_size, range_size - offset), bits)
        for offset in range(0, range_size, chunk_size)
    ]
    records = []
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for index, chunk in enumerate(pool.map(worker, jobs), start=1):
            records.extend(chunk)
            if index % max(1, len(jobs) // 10) == 0 or index == len(jobs):
                print(f"  chunks: {index}/{len(jobs)}")

    records.sort()
    with open(output_path, "wb") as handle:
        handle.write(HEADER.pack(
            MAGIC, range_start.to_bytes(32, "big"), range_size, bits,
            len(records), RECORD.size
        ))
        for key, offset in records:
            handle.write(RECORD.pack(key, offset))

    elapsed = time.perf_counter() - started
    return records, elapsed, os.path.getsize(output_path)


def lookup(records, key):
    # The table is sorted by (fingerprint, offset); binary-search the equal-key run.
    first = bisect.bisect_left(records, (key, -1))
    matches = []
    for stored_key, offset in records[first:]:
        if stored_key == key:
            matches.append(offset)
        else:
            break
    return matches


def read_table(path):
    with open(path, "rb") as handle:
        header = handle.read(HEADER.size)
        magic, start_bytes, range_size, bits, record_count, record_size = HEADER.unpack(header)
        if magic != MAGIC or record_size != RECORD.size:
            raise ValueError("unsupported angle table format")
        payload = handle.read()
    if len(payload) != record_count * RECORD.size:
        raise ValueError("truncated angle table")
    records = [RECORD.unpack_from(payload, index)
               for index in range(0, len(payload), RECORD.size)]
    return {
        "range_start": int.from_bytes(start_bytes, "big"),
        "range_size": range_size,
        "bits": bits,
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--range-start", type=int, default=None)
    parser.add_argument("--bit-length", type=int, default=None,
                        help="use canonical bit-length window starting at 2^(bit_length-1)")
    parser.add_argument("--range-size", type=int, default=100000)
    parser.add_argument("--bits", type=int, default=16)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--queries", type=int, default=100)
    parser.add_argument("--target-public-key", default=None,
                        help="compressed/uncompressed public key to recover")
    parser.add_argument("--target-scalar", type=int, default=None,
                        help="test-only helper: derive the query public key")
    parser.add_argument("--output", default="/tmp/parallel_angle_table.bin")
    parser.add_argument("--load-table", default=None,
                        help="load an existing table instead of rebuilding it")
    args = parser.parse_args()

    if args.bit_length is not None and not 1 <= args.bit_length <= 256:
        parser.error("--bit-length must be in 1..256")
    if args.range_start is None:
        args.range_start = (1 << (args.bit_length - 1)) if args.bit_length else (1 << 255)
    if args.range_start < 1 or args.range_size <= 0:
        parser.error("invalid scalar range")
    if args.range_start + args.range_size > GROUP_ORDER:
        parser.error("range must be contained in [1, n)")
    if not 1 <= args.bits <= 32:
        parser.error("--bits must be in 1..32 (two coordinates fit in 64 bits)")
    if args.target_public_key and args.target_scalar is not None:
        parser.error("use either --target-public-key or --target-scalar")
    if args.workers <= 0 or args.chunk_size <= 0 or args.queries <= 0:
        parser.error("workers and chunk-size must be positive")

    print("Parallel bounded angle-table recovery")
    build_seconds = None
    if args.load_table:
        table_path = args.load_table
        table_data = read_table(table_path)
        table_start = table_data["range_start"]
        table_size = table_data["range_size"]
        table_bits = table_data["bits"]
        records = table_data["records"]
        file_bytes = os.path.getsize(table_path)
        print(f"  loaded table: {table_path}")
    else:
        table_path = args.output
        records, build_seconds, file_bytes = build_table(
            args.range_start, args.range_size, args.bits, args.workers,
            args.chunk_size, table_path
        )
        table_data = read_table(table_path)
        table_start = args.range_start
        table_size = args.range_size
        table_bits = args.bits
    print(f"  range start: 0x{table_start:064x}")
    print(f"  range size: {table_size}")
    print(f"  covered interval end: {table_start + table_size - 1:#066x}")
    canonical_space_bits = table_start.bit_length()
    print(f"  coverage fraction of 2^{canonical_space_bits}: "
          f"{table_size / (1 << canonical_space_bits):.6e}")
    print(f"  angle bits per coordinate: {table_bits}")
    if args.target_scalar is not None:
        if not table_start <= args.target_scalar < table_start + table_size:
            parser.error("target scalar must be inside the table range")
        query_point = point_from_scalar(args.target_scalar)
        public_key = compressed_public_key(query_point)
        query_offsets = [args.target_scalar - table_start]
    elif args.target_public_key:
        query_point = point_from_public_key(args.target_public_key)
        public_key = compressed_public_key(query_point)
        query_offsets = []
    else:
        rng = random.Random(args.seed)
        query_offsets = [rng.randrange(table_size) for _ in range(args.queries)]
        query_point = None
        public_key = None

    analyzer = DomeInvariantAnalyzer()
    bucket_sizes = []
    exact_counts = []
    correct = 0
    for target_offset in query_offsets:
        target_point = point_from_scalar(table_start + target_offset)
        target_values = angle_pair(analyzer, target_point)
        candidates = lookup(table_data["records"], fingerprint(target_values, table_bits))
        exact_matches = []
        for offset in candidates:
            candidate_point = point_from_scalar(table_start + offset)
            if (angle_pair(analyzer, candidate_point) == target_values
                    and candidate_point.x == target_point.x
                    and candidate_point.y == target_point.y):
                exact_matches.append(offset)
        bucket_sizes.append(len(candidates))
        exact_counts.append(len(exact_matches))
        correct += exact_matches == [target_offset]

    if args.target_public_key:
        target_values = angle_pair(analyzer, query_point)
        candidates = lookup(table_data["records"], fingerprint(target_values, table_bits))
        exact_matches = []
        for offset in candidates:
            candidate_point = point_from_scalar(table_start + offset)
            if (angle_pair(analyzer, candidate_point) == target_values
                    and candidate_point.x == query_point.x
                    and candidate_point.y == query_point.y):
                exact_matches.append(offset)
        bucket_sizes.append(len(candidates))
        exact_counts.append(len(exact_matches))
        correct = len(exact_matches) == 1
        query_offsets = exact_matches
        recovered_scalars = [table_start + offset for offset in exact_matches]
    else:
        recovered_scalars = []

    print(f"  table file: {table_path}")
    print(f"  table bytes: {file_bytes}")
    if build_seconds is not None:
        print(f"  build seconds: {build_seconds:.2f}")
    print(f"  queries: {len(bucket_sizes)}")
    if public_key:
        print(f"  query public key: {public_key}")
    if args.target_public_key:
        print(f"  recovered offsets: {query_offsets}")
        print("  recovered private scalars:")
        for scalar in recovered_scalars:
            print(f"    decimal: {scalar}")
            print(f"    hex:    0x{scalar:064x}")
    print(f"  mean fingerprint bucket: {sum(bucket_sizes) / len(bucket_sizes):.3f}")
    print(f"  maximum fingerprint bucket: {max(bucket_sizes)}")
    print(f"  exact one-candidate recoveries: {int(bool(correct))}/{len(bucket_sizes)}"
          if args.target_public_key else
          f"  exact one-candidate recoveries: {correct}/{len(bucket_sizes)}")
    print(f"  exact-match multiplicities: {sorted(set(exact_counts))}")


if __name__ == "__main__":
    main()
