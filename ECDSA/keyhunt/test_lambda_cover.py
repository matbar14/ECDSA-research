#!/usr/bin/env python3
"""Generate a 135-bit key and recover it through automatic shifted coverage."""

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
    parser.add_argument("--binary", default="./eng_bsgs")
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    local_start = 1 << 134
    local_size = 4096
    cover_start = local_start + (1 << 20)
    cover_size = local_size * 3 + 513
    scalar = cover_start + local_size * 3 + rng.randrange(513)
    public_key = compressed_public_key(scalar)
    command = [
        args.binary,
        "--target", public_key,
        "--start", hex(local_start),
        "--size", str(local_size),
        "--baby-rows", "64",
        "--threads", "2",
        "--cover-start", hex(cover_start),
        "--cover-size", str(cover_size),
        "--lambda-priority",
        "--lambda-sectors", "16",
        "--glv-orbit",
    ]
    result = subprocess.run(command, text=True, capture_output=True)
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    match = re.search(r"recovered private key: 0x([0-9a-fA-F]{64})", result.stdout)
    recovered = int(match.group(1), 16) if match else None
    print(f"generated public key: {public_key}")
    print(f"generated scalar: 0x{scalar:064x}")
    print(f"covered interval: [0x{cover_start:x}, 0x{cover_start + cover_size:x})")
    print(f"local BSGS window: [0x{local_start:x}, 0x{local_start + local_size:x})")
    print(f"recovery correct: {result.returncode == 0 and recovered == scalar}")
    raise SystemExit(0 if result.returncode == 0 and recovered == scalar else 1)


if __name__ == "__main__":
    main()
