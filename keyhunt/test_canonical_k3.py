#!/usr/bin/env python3
"""Black-box test for --table-scalars + --canonical-k3: real multi-row buckets.

The baby table is built from an explicit scalar list that contains all three
x-representatives of one orbit (j3, lambda*j3, lambda^2*j3) plus background
rows, so a giant hitting that orbit sees a genuine 3-row K3 bucket - the
full-range scenario modeled by phase_signature_economy_experiment.py.

Cases:
  [1] 6-branch baseline on the 3-row bucket (partner target out of window)
  [2] phase-signature routing: all rows collapse into one candidate pair
  [3] --canonical-k3: the 3 rows merge into one canonical row at build time
  [4] --canonical-k3 in-window regression (partner-table off)
  [5] --canonical-k3 on the continuous-window table (no injected collisions)
"""

import argparse
import random
import re
import subprocess
import tempfile
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
LAMBDA2 = 0xAC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

START = 1 << 134
SIZE = 1 << 20
J3 = 777
BACKGROUND = 6


def pubkey_hex(scalar):
    numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def run(binary, target, size, extra):
    cmd = [binary, "--target", target, "--start", hex(START),
           "--size", str(size), "--segment", str(size), "--threads", "2"] + extra
    result = subprocess.run(cmd, text=True, capture_output=True)
    match = re.search(r"recovered private key: 0x([0-9a-fA-F]+)",
                      result.stdout)
    recovered = int(match.group(1), 16) if match else None
    counter = {}
    for name in ("signature routes", "signature pair collapses",
                 "signature pair mismatches", "canonical rows merged",
                 "partner recoveries", "exact candidate checks",
                 "orbit factors scheduled", "table rows stored"):
        m = re.search(rf"{re.escape(name)}: (\d+)", result.stdout)
        if m:
            counter[name] = int(m.group(1))
    if result.returncode not in (0, 1):
        print(result.stdout)
        print(result.stderr)
    return recovered, counter, result


def write_rows_file(directory):
    rows = [J3, LAMBDA * J3 % N, LAMBDA2 * J3 % N]
    rng = random.Random(20261031)
    rows += [rng.randrange(1, 1 << 100) for _ in range(BACKGROUND)]
    path = Path(directory) / "table_scalars.txt"
    lines = ["# orbit representatives of j3 plus background rows"]
    lines += [hex(r) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return path, len(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    args = parser.parse_args()

    failures = []
    with tempfile.TemporaryDirectory(
            dir="/home/z/my-project") as tmpdir:
        rows_path, raw_rows = write_rows_file(tmpdir)
        merged_rows = raw_rows - 2  # 3 representatives of one orbit -> 1 row

        def table_file():
            return ["--table-scalars", str(rows_path)]

        # 1. baseline: 6-branch check on the multi-row bucket, partner target
        i0 = 213
        partner = (START + i0 * raw_rows + LAMBDA * J3) % N
        rec, ctr, _ = run(args.binary, pubkey_hex(partner), SIZE,
                          ["--glv-orbit", "--partner-table"] + table_file())
        checks1 = ctr.get("exact candidate checks", -1)
        factors1 = ctr.get("orbit factors scheduled", -1)
        ok = (rec == partner and ctr.get("partner recoveries", 0) >= 1
              and 1 <= checks1 <= 6 and 1 <= factors1 <= 6)
        print(f"[1] 6-branch baseline, 3-row bucket: correct={rec == partner}, "
              f"partner={ctr.get('partner recoveries', '?')}, "
              f"branches={factors1}, exact={checks1} -> "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append(1)

        # 2a. signature routing, partner target, NO --partner-table: the
        #     routed pair leaves the window, every later row of the bucket
        #     routes to the SAME pair and collapses (negative control).
        rec, ctr, res = run(args.binary, pubkey_hex(partner), SIZE,
                            ["--glv-orbit", "--phase-signature"] + table_file())
        ok = (rec is None and res.returncode == 1
              and ctr.get("signature routes", 0) == 3
              and ctr.get("signature pair collapses", 0) == 2
              and ctr.get("signature pair mismatches", 0) == 0
              and ctr.get("exact candidate checks", -1) == 0)
        print(f"[2a] signature collapse (pair out of window): "
              f"recovered={rec is not None}, "
              f"routes={ctr.get('signature routes', '?')}, "
              f"collapses={ctr.get('signature pair collapses', '?')}, "
              f"mismatches={ctr.get('signature pair mismatches', '?')}, "
              f"exact={ctr.get('exact candidate checks', '?')} -> "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append("2a")

        # 2b. signature routing with --partner-table: the first routed row
        #     already contains the giant scalar in its pair, so the target is
        #     recovered within at most two exact checks.
        rec, ctr, _ = run(args.binary, pubkey_hex(partner), SIZE,
                          ["--glv-orbit", "--partner-table",
                           "--phase-signature"] + table_file())
        ok = (rec == partner
              and ctr.get("signature routes", 0) >= 1
              and ctr.get("signature pair mismatches", 0) == 0
              and 1 <= ctr.get("exact candidate checks", 99) <= 2)
        print(f"[2b] signature partner recovery: correct={rec == partner}, "
              f"routes={ctr.get('signature routes', '?')}, "
              f"exact={ctr.get('exact candidate checks', '?')} -> "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append("2b")

        # 3. canonical storage: three rows merge into one canonical row
        partner_canon = (START + i0 * merged_rows + LAMBDA * J3) % N
        rec, ctr, _ = run(args.binary, pubkey_hex(partner_canon), SIZE,
                          ["--glv-orbit", "--partner-table",
                           "--phase-signature", "--canonical-k3"] + table_file())
        ok = (rec == partner_canon
              and ctr.get("canonical rows merged", 0) == 2
              and ctr.get("table rows stored", 0) == merged_rows
              and ctr.get("signature routes", 0) == 1
              and ctr.get("signature pair collapses", 0) == 0
              and ctr.get("signature pair mismatches", 0) == 0
              and ctr.get("partner recoveries", 0) >= 1
              and 1 <= ctr.get("exact candidate checks", 99) <= 2)
        print(f"[3] canonical storage, merged bucket: correct={rec == partner_canon}, "
              f"merged={ctr.get('canonical rows merged', '?')}, "
              f"stored={ctr.get('table rows stored', '?')}, "
              f"routes={ctr.get('signature routes', '?')}, "
              f"exact={ctr.get('exact candidate checks', '?')} -> "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append(3)

        # 4. canonical storage, in-window target, no partner-table
        window_scalar = START + i0 * merged_rows + J3
        rec, ctr, _ = run(args.binary, pubkey_hex(window_scalar), SIZE,
                          ["--glv-orbit", "--phase-signature",
                           "--canonical-k3"] + table_file())
        ok = (rec == window_scalar
              and ctr.get("canonical rows merged", 0) == 2
              and ctr.get("partner recoveries", 0) == 0
              and ctr.get("exact candidate checks", 99) == 1)
        print(f"[4] canonical in-window regression: correct={rec == window_scalar}, "
              f"exact={ctr.get('exact candidate checks', '?')} -> "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append(4)

        # 5. canonical storage on the continuous-window table: no collisions
        rng = random.Random(20260926)
        scalar = START + rng.randrange(SIZE)
        rec, ctr, _ = run(args.binary, pubkey_hex(scalar), SIZE,
                          ["--glv-orbit", "--phase-signature",
                           "--canonical-k3", "--baby-rows", "1024"])
        ok = (rec == scalar and ctr.get("canonical rows merged", -1) == 0
              and ctr.get("exact candidate checks", 99) == 1)
        print(f"[5] canonical continuous-window regression: "
              f"correct={rec == scalar}, "
              f"merged={ctr.get('canonical rows merged', '?')} -> "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append(5)

    print("ALL PASS" if not failures else f"FAILED CASES: {failures}")
    raise SystemExit(0 if not failures else 1)


if __name__ == "__main__":
    main()
