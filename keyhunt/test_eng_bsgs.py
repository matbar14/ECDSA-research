#!/usr/bin/env python3
"""Black-box regression test for eng_bsgs.

The private scalar is generated in this harness and is never passed to the
binary. eng_bsgs receives only the public key and the interval.
"""

import argparse
import random
import re
import subprocess

from cryptography.hazmat.primitives.asymmetric import ec


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--range-size", type=int, default=3_000_000_000_000)
    parser.add_argument("--baby-rows", type=int, default=300_000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--glv-orbit", action="store_true")
    parser.add_argument("--phase-signature", action="store_true",
                        help="route GLV branches by the (d01,d12) signature")
    args = parser.parse_args()

    start = 1 << 134
    rng = random.Random(args.seed)
    offset = rng.randrange(args.range_size)
    scalar = start + offset
    numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
    public_key = f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"

    command = [
        args.binary,
        "--target", public_key,
        "--start", hex(start),
        "--size", str(args.range_size),
        "--baby-rows", str(args.baby_rows),
        "--threads", str(args.threads),
    ]
    if args.glv_orbit:
        command.append("--glv-orbit")
    if args.phase_signature:
        if not args.glv_orbit:
            raise SystemExit("--phase-signature requires --glv-orbit")
        command.append("--phase-signature")
    result = subprocess.run(command, text=True, capture_output=True)
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    match = re.search(r"recovered private key: 0x([0-9a-fA-F]{64})", result.stdout)
    recovered = int(match.group(1), 16) if match else None
    print(f"black-box target public key: {public_key}")
    print(f"black-box expected scalar: 0x{scalar:064x}")
    print(f"black-box target offset: {offset}")
    print(f"black-box recovery correct: {recovered == scalar}")
    if args.phase_signature:
        routes = re.search(r"signature routes: (\d+)", result.stdout)
        anomalies = re.search(r"signature anomalies \(6-branch fallback\): (\d+)",
                              result.stdout)
        print(f"black-box signature routes: {routes.group(1) if routes else '?'}")
        print(f"black-box signature anomalies: "
              f"{anomalies.group(1) if anomalies else '?'}")
    if result.returncode != 0 or recovered != scalar:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
