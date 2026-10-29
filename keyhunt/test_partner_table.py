#!/usr/bin/env python3
"""Black-box test for --partner-table: recovery of out-of-window orbit partners.

The target scalar is built as k = start + i0*m + lambda*j0 (mod n): it lies
far outside the searched window, but its giant-step image equals
lambda*j0*G, whose K3=x^3 matches table row j0. Without --partner-table the
engine skips the out-of-interval branch (negative control); with it the
branch is exact-verified against the target and the scalar is recovered.
"""

import argparse
import random
import re
import subprocess

from cryptography.hazmat.primitives.asymmetric import ec

LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

START = 1 << 134
SIZE = 1 << 20
M = 1024


def pubkey_hex(scalar):
    numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def run(binary, target, extra):
    cmd = [binary, "--target", target, "--start", hex(START),
           "--size", str(SIZE), "--segment", str(SIZE), "--threads", "2"] + extra
    result = subprocess.run(cmd, text=True, capture_output=True)
    match = re.search(r"recovered private key: 0x([0-9a-fA-F]+)",
                      result.stdout)
    return result, (int(match.group(1), 16) if match else None)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    parser.add_argument("--seed", type=int, default=20261030)
    args = parser.parse_args()
    rng = random.Random(args.seed)

    i0 = rng.randrange(50, 400)
    j0 = rng.randrange(1, M)
    partner_scalar = (START + i0 * M + LAMBDA * j0) % N
    window_scalar = START + i0 * M + j0
    failures = []

    # 1. negative control: partner key must NOT be recovered without the flag
    result, rec = run(args.binary, pubkey_hex(partner_scalar), ["--glv-orbit"])
    ok = rec is None and result.returncode != 0
    print(f"[1] negative control (no --partner-table): "
          f"recovered={rec is not None} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        failures.append(1)

    # 2. partner recovery through the six-branch check
    result, rec = run(args.binary, pubkey_hex(partner_scalar),
                      ["--glv-orbit", "--partner-table"])
    pr = re.search(r"partner recoveries: (\d+)", result.stdout)
    ok = rec == partner_scalar
    print(f"[2] partner recovery (6-branch): correct={ok}, "
          f"partner recoveries={pr.group(1) if pr else '?'} -> "
          f"{'PASS' if ok else 'FAIL'}")
    if not ok:
        failures.append(2)

    # 3. partner recovery through the phase-signature router
    result, rec = run(args.binary, pubkey_hex(partner_scalar),
                      ["--glv-orbit", "--phase-signature", "--partner-table"])
    routes = re.search(r"signature routes: (\d+)", result.stdout)
    anom = re.search(r"signature anomalies \(6-branch fallback\): (\d+)",
                     result.stdout)
    ok = rec == partner_scalar
    print(f"[3] partner recovery (signature): correct={ok}, "
          f"routes={routes.group(1) if routes else '?'}, "
          f"anomalies={anom.group(1) if anom else '?'} -> "
          f"{'PASS' if ok else 'FAIL'}")
    if not ok:
        failures.append(3)

    # 4. regression: in-window target is still recovered with the flag on
    result, rec = run(args.binary, pubkey_hex(window_scalar),
                      ["--glv-orbit", "--phase-signature", "--partner-table"])
    ok = rec == window_scalar
    print(f"[4] in-window regression with --partner-table: correct={ok} -> "
          f"{'PASS' if ok else 'FAIL'}")
    if not ok:
        failures.append(4)

    print(f"partner-table black-box target scalar: 0x{partner_scalar:064x}")
    print(f"partner-table black-box window row j0={j0}, giant index i0={i0}")
    print("ALL PASS" if not failures else f"FAILED CASES: {failures}")
    raise SystemExit(0 if not failures else 1)


if __name__ == "__main__":
    main()
