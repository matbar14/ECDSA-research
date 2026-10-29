#!/usr/bin/env python3
"""Black-box верификация файла представителей полного диапазона.

Проверяет, что persistent-файл канонических представителей
(representatives_canonical_H*.txt, см. make_representatives_file.py)
работает через движок end-to-end:

  [1] каноничность строк: для выборки строк файла канонический
      представитель головы j совпадает с сохранённой строкой;
  [2] partner-recovery: цель k = start + i*H + s*,
      s* = +/-lambda^a * c (c - строка файла), восстанавливается
      --table-scalars FILE --glv-orbit --phase-signature --partner-table,
      счётчик partner recoveries растёт, восстановленный скаляр точен;
  [3] core-tiling recovery: цель k = start + i*H + j0 (ядро [1, H])
      восстанавливается той же композицией;
  [4] negative control: та же partner-цель БЕЗ --partner-table не
      восстанавливается (интервальный пре-фильтр отбрасывает ветви вне
      локального окна) - partner-покрытие файла действительно требует
      partner-expansion.

Харнесс передаёт движку только public key/start/size - private key
используется лишь для сверки внутри тестового процесса.
"""

import argparse
import random
import re
import subprocess
import time
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
START = 1 << 134


def canonical_representative(j):
    best = None
    v = j % N
    for _ in range(3):
        candidate = min(v, N - v) if v else 0
        if best is None or candidate < best:
            best = candidate
        v = v * LAMBDA % N
    return best


def pubkey_hex(scalar):
    numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def run_engine(binary, target, size, extra):
    cmd = [binary, "--target", target, "--start", hex(START),
           "--size", str(size), "--segment", str(size), "--threads", "2"] + extra
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, text=True, capture_output=True)
    wall = time.perf_counter() - t0
    match = re.search(r"recovered private key: 0x([0-9a-fA-F]+)", proc.stdout)
    recovered = int(match.group(1), 16) if match else None
    counters = {}
    for name in ("table rows stored", "partner recoveries",
                 "signature routes", "exact candidate checks"):
        m = re.search(rf"{re.escape(name)}: (\d+)", proc.stdout)
        if m:
            counters[name] = int(m.group(1))
    return {"recovered": recovered, "counters": counters,
            "wall": wall, "stdout": proc.stdout, "stderr": proc.stderr}


def load_rows(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            rows.append(int(line, 16))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    parser.add_argument("--file", default="../representatives_canonical_H50000.txt")
    parser.add_argument("--partner-trials", type=int, default=4)
    parser.add_argument("--core-trials", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20261030)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    rows = load_rows(args.file)
    H = len(rows)
    stride = H                      # движок стридит по числу строк файла
    size = stride * 4096            # окно проверки: 4096 гигант-позиций
    print(f"file: {args.file}")
    print(f"rows: {H}, check window: {size} ({4096} giant positions)")

    failures = []

    # [1] каноничность выборки строк
    sample = rng.sample(range(H), min(200, H))
    bad = sum(1 for idx in sample
              if canonical_representative(idx + 1) != rows[idx])
    print(f"[1] canonical spot-check: {len(sample) - bad}/{len(sample)} rows "
          f"match c(j)")
    if bad:
        failures.append(f"canonicality: {bad} rows mismatch")

    # [2] partner-recovery: s* = +/-lambda^a * c
    base_extra = ["--glv-orbit", "--phase-signature", "--partner-table",
                  "--table-scalars", args.file]
    for _ in range(args.partner_trials):
        idx = rng.randrange(H)
        c = rows[idx]
        a = rng.randrange(3)
        sign = rng.choice((1, -1))
        i = rng.randrange(2048, 4096)
        s = sign * pow(LAMBDA, a, N) * c % N
        scalar = (START + i * stride + s) % N
        run = run_engine(args.binary, pubkey_hex(scalar), size, base_extra)
        ok = run["recovered"] == scalar
        partners = run["counters"].get("partner recoveries", 0)
        print(f"[2] partner r={a} sign={sign:+d} row_idx={idx} -> "
              f"{'PASS' if ok else 'FAIL'}, partner recoveries={partners}, "
              f"wall={run['wall']:.1f}s")
        if not ok:
            failures.append(f"partner idx={idx} a={a} sign={sign}")

    # [3] core-tiling recovery: s* = j0 из ядра [1, H]
    for _ in range(args.core_trials):
        j0 = rng.randrange(1, H)
        i = rng.randrange(2048, 4096)
        scalar = START + i * stride + j0
        run = run_engine(args.binary, pubkey_hex(scalar), size, base_extra)
        ok = run["recovered"] == scalar
        print(f"[3] core j0={j0} -> {'PASS' if ok else 'FAIL'}, "
              f"wall={run['wall']:.1f}s")
        if not ok:
            failures.append(f"core j0={j0}")

    # [4] negative control: partner-цель без --partner-table
    idx = rng.randrange(H)
    c = rows[idx]
    i = rng.randrange(2048, 4096)
    s = pow(LAMBDA, 1, N) * c % N
    scalar = (START + i * stride + s) % N
    run = run_engine(args.binary, pubkey_hex(scalar), size,
                     ["--glv-orbit", "--phase-signature",
                      "--table-scalars", args.file])
    if run["recovered"] is None:
        print(f"[4] negative control (no --partner-table): PASS "
              f"(не восстановлено, как и ожидалось)")
    else:
        print(f"[4] negative control: FAIL (неожиданное восстановление)")
        failures.append("negative control recovered unexpectedly")

    if failures:
        print(f"\nFAILURES: {len(failures)}")
        for f in failures:
            print(" -", f)
        raise SystemExit(1)
    print("\nALL PASS")


if __name__ == "__main__":
    main()
