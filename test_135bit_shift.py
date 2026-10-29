#!/usr/bin/env python3
"""End-to-end 135-bit scalar test for eng_bsgs shifted targets."""

import argparse
import random
import re
import subprocess

from cryptography.hazmat.primitives.asymmetric import ec


def compressed_public_key(scalar):
    numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./keyhunt/eng_bsgs")
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--window-size", type=int, default=4096)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    lower = 1 << 134
    upper = 1 << 135
    scalar = rng.randrange(lower, upper)
    window_start = lower
    local_offset = rng.randrange(args.window_size)
    shift = window_start + local_offset - scalar
    public_key = compressed_public_key(scalar)
    command = [
        args.binary,
        "--target", public_key,
        "--start", hex(window_start),
        "--size", str(args.window_size),
        "--segment", str(args.window_size),
        "--baby-rows", str(max(32, int(args.window_size ** 0.5))),
        f"--shift={shift}",
        "--threads", "1",
    ]
    result = subprocess.run(command, text=True, capture_output=True)
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    match = re.search(r"recovered private key: 0x([0-9a-fA-F]{64})", result.stdout)
    recovered = int(match.group(1), 16) if match else None
    print(f"test scalar: 0x{scalar:064x}")
    print(f"test public key: {public_key}")
    print(f"window: [0x{window_start:x}, 0x{window_start + args.window_size:x})")
    print(f"known test shift: {shift}")
    print(f"expected shifted scalar: 0x{(scalar + shift):x}")
    print(f"recovery correct: {result.returncode == 0 and recovered == scalar}")
    if result.returncode != 0 or recovered != scalar:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
