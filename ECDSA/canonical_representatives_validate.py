#!/usr/bin/env python3
"""Валидация файла представителей на движке: одна орбитальная цель
(k = START + i*stride + s*, s* = +/-lambda^a*j, j in [1, H]) на каноническом
режиме с маршрутизацией и partner-table, окно (3H)^2, stride = H строк."""
import argparse
import json
import random
import re
import subprocess
import threading
import time
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
START = 1 << 134


def pubkey_hex(scalar):
    numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default="./eng_bsgs")
    ap.add_argument("--table", required=True)
    ap.add_argument("--heads", type=int, required=True)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()

    H = args.heads
    stride = H
    size = (3 * H) ** 2
    giant_count = (size - 1) // stride + 1

    rng = random.Random(20261102)
    j = rng.randrange(1, H + 1)
    a = rng.randrange(3)
    sign = rng.choice((1, -1))
    i = rng.randrange(giant_count // 2, giant_count)
    s = sign * pow(LAMBDA, a, N) * j % N
    scalar = (START + i * stride + s) % N
    print(f"target: head={j} rotation={a} sign={sign:+d} i={i} "
          f"stride={stride} window={size:.3e}")

    cmd = [args.binary, "--target", pubkey_hex(scalar), "--start", hex(START),
           "--size", str(size), "--segment", str(size),
           "--threads", str(args.threads),
           "--table-scalars", args.table,
           "--glv-orbit", "--phase-signature", "--partner-table"]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, text=True, capture_output=True)
    wall = time.perf_counter() - t0

    match = re.search(r"recovered private key: 0x([0-9a-fA-F]+)", proc.stdout)
    recovered = int(match.group(1), 16) if match else None
    counters = {}
    for name in ("table rows stored", "canonical rows merged", "giant steps",
                 "invariant hits", "exact candidate checks",
                 "EC scalar multiplications", "EC point additions",
                 "total EC operations", "signature routes",
                 "partner recoveries"):
        m = re.search(rf"{re.escape(name)}: (\d+)", proc.stdout)
        if m:
            counters[name] = int(m.group(1))
    ok = recovered == scalar
    if not ok:
        print("== stdout tail ==")
        print("\n".join(proc.stdout.splitlines()[-15:]))
        print("== stderr ==")
        print(proc.stderr)
        print(f"== returncode: {proc.returncode} ==")
    for line in proc.stdout.splitlines():
        if any(k in line for k in ("table rows stored", "EC scalar",
                                   "EC point", "total EC", "giant steps",
                                   "table build seconds", "total seconds")):
            print("  " + line)
    print(f"recovered: {'YES' if recovered is not None else 'NO'}, "
          f"exact match: {ok}, wall {wall:.2f}s")
    result = {"heads": H, "head": j, "rotation": a, "sign": sign,
              "stride index": i, "correct": ok, "wall seconds": round(wall, 2),
              "counters": counters}
    out = Path("../representative_file_validation.json")
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(f"[saved] {out}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
