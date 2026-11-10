#!/usr/bin/env python3
"""Generate a random bounded key and recover its offset with segmented BSGS."""

import argparse
import random

from cryptography.hazmat.primitives.asymmetric import ec

from segmented_bsgs_recovery import recover
from w1 import GROUP_ORDER


def compressed_public_key(scalar):
    numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bit-length", type=int, default=135)
    parser.add_argument("--range-size", type=int, default=100000)
    parser.add_argument("--segment-size", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--target-public-key", default=None,
                        help="recover this public key instead of generating a control target")
    args = parser.parse_args()
    range_start = 1 << (args.bit_length - 1)
    if range_start + args.range_size > GROUP_ORDER:
        parser.error("range exceeds group order")
    if args.target_public_key:
        from segmented_bsgs_recovery import point_from_public_key, recover_point
        target_point = point_from_public_key(args.target_public_key)
        target_offset = None
        target_scalar = None
        report = recover_point(range_start, args.range_size, args.segment_size, target_point)
    else:
        target_offset = random.Random(args.seed).randrange(args.range_size)
        target_scalar = range_start + target_offset
        report = recover(range_start, args.range_size, args.segment_size, target_offset)
    print("Random bounded segmented-BSGS recovery")
    print(f"  seed: {args.seed}")
    print(f"  interval: [2^{args.bit_length - 1}, 2^{args.bit_length - 1} + {args.range_size})")
    if target_scalar is not None:
        print(f"  target public key: {compressed_public_key(target_scalar)}")
        print(f"  target scalar (control): 0x{target_scalar:064x}")
        print(f"  target offset (control): {target_offset}")
    else:
        print(f"  target public key: {args.target_public_key}")
    print(f"  recovered offset: {report['recovered_offset']}")
    print(f"  recovered scalar: 0x{range_start + report['recovered_offset']:064x}"
          if report['recovered_offset'] is not None else "  recovered scalar: None")
    print(f"  segments: {report['segments']}")
    print(f"  total baby rows built: {report['total_baby_rows_built']}")
    print(f"  giant steps: {report['giant_steps']}")
    print(f"  seconds: {report['seconds']:.3f}")
    print(f"  recovery found a scalar in range: {report['correct']}")
    if target_offset is not None:
        print(f"  recovery correct: {report['correct'] and report['recovered_offset'] == target_offset}")


if __name__ == "__main__":
    main()
