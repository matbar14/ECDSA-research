#!/usr/bin/env python3
"""CBT1 H=4M, file-mode: сессия на бинарной таблице против текстового пути.

  [R] recovery в file-mode: контрольная цель build-прогона
      (start + 2*stride + 777) восстанавливается через --table-binary;
      stride == 4M, merged rows 0, insert collisions 0, exact checks 1;
  [C] стоимость старта сессии: "table build seconds" в file-mode - это
      mmap + валидация + checksum, БЕЗ 4M EC-умножений (build mults ~1
      против 4,000,000 в RAM-режиме);
  [S] парный скан RAM vs mmap: 20M гигант-шагов на дальнюю цель вне
      окна, steps/s по тикеру (прогресс каждые 5 c) + steady-RSS.

Артефакт пишется атомарно (resume-безопасность среды).
"""

import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ENGINE = ROOT / "keyhunt" / "eng_bsgs"
TABLE = ROOT / "tables" / "table_cbt1_H4000000.cbt1"
CANON = ROOT / "tables" / "representatives_canonical_H4000000.txt"
ARTIFACT = ROOT / "results" / "cbt1_4m_filemode.json"

START = 0x4000000000000000000000000000000000
STRIDE = 4_000_000
TARGET = START + 2 * STRIDE + 777
FAR_TARGET = START + 1_000_000_000_000_000
SCAN_SIZE = 80_000_000_000_000         # 20M гигант-шагов при стриде 4M


def pubkey_hex(scalar: int) -> str:
    from cryptography.hazmat.primitives.asymmetric import ec
    numbers = ec.derive_private_key(scalar, ec.SECP256K1()) \
        .public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def run(cmd, timeout=300):
    t0 = time.monotonic()
    proc = subprocess.run(cmd, text=True, capture_output=True, timeout=timeout)
    return proc, time.monotonic() - t0


def grab(out, pattern):
    m = re.search(pattern, out)
    return m.group(1) if m else None


def ticker_samples(out):
    samples = []
    for m in re.finditer(
            r"progress: elapsed ([0-9.]+)s \| giant steps \d+ \| "
            r"steps/s ([0-9.]+) \| rss kB (\d+)", out):
        samples.append({"elapsed": float(m.group(1)),
                        "steps/s": float(m.group(2)),
                        "rss kB": int(m.group(3))})
    return samples


def main():
    if ARTIFACT.exists():
        print(f"[skip] artifact exists: {ARTIFACT}")
        return
    if not TABLE.exists():
        raise SystemExit(f"[E] missing table: {TABLE}")

    base = [str(ENGINE), "--start", hex(START), "--threads", "2"]

    # ---- [R]+[C] file-mode recovery and session start cost ----
    cmd_r = base + [
        "--target", pubkey_hex(TARGET),
        "--size", "500000000000", "--segment", "500000000000",
        "--progress-interval", "0",
        "--table-binary", str(TABLE), "--glv-orbit", "--partner-table",
    ]
    proc_r, wall_r = run(cmd_r)
    if proc_r.returncode != 0:
        print(proc_r.stdout[-1500:], proc_r.stderr[-800:], sep="\n")
        raise SystemExit(f"[E] file-mode run exit {proc_r.returncode}")
    rec = grab(proc_r.stdout, r"recovered private key: 0x([0-9a-f]+)")
    rec_ok = rec is not None and int(rec, 16) == TARGET
    stride_r = grab(proc_r.stdout, r"giant stride \(post-merge rows\): (\d+)")
    mmap_backend = "table backend: mmap" in proc_r.stdout
    cnt = {name: grab(proc_r.stdout, rf"{name}: (\d+)")
           for name in ("merged rows", "insert collisions",
                        "exact candidate checks", "invariant hits",
                        "EC scalar multiplications")}
    build_secs_r = grab(proc_r.stdout, r"table build seconds: ([0-9.]+)")
    build_mults_r = grab(proc_r.stdout,
                         r"EC scalar multiplications \(build\): (\d+)")

    print(f"[R] recovery ok={rec_ok}, stride={stride_r}, backend mmap={mmap_backend}")
    print(f"[C] session start (mmap+validate): {build_secs_r} s, "
          f"build EC mults: {build_mults_r}")

    # ---- [S] paired RAM vs mmap scan, 20M giant steps ----
    rates = {}
    for tag, extra in (
        ("ram", ["--target", pubkey_hex(FAR_TARGET),
                 "--size", str(SCAN_SIZE), "--segment", str(SCAN_SIZE),
                 "--progress-interval", "5",
                 "--table-scalars", str(CANON),
                 "--glv-orbit", "--canonical-k3"]),
        ("mmap", ["--target", pubkey_hex(FAR_TARGET),
                  "--size", str(SCAN_SIZE), "--segment", str(SCAN_SIZE),
                  "--progress-interval", "5",
                  "--table-binary", str(TABLE), "--glv-orbit"]),
    ):
        cmd = base + extra
        proc, wall = run(cmd, timeout=300)
        # Exit 1 is the engine's "target not found in bounded interval" -
        # the expected outcome for a far target outside the scan window.
        if proc.returncode not in (0, 1):
            print(proc.stdout[-1500:], proc.stderr[-500:], sep="\n")
            raise SystemExit(f"[E] {tag} scan exit {proc.returncode}")
        samples = ticker_samples(proc.stdout)
        build_secs = grab(proc.stdout, r"table build seconds: ([0-9.]+)")
        rss_build = grab(proc.stdout, r"RSS after table build: (\d+) kB")
        hwm_build = grab(proc.stdout, r"RSS after table build: \d+ kB \(peak (\d+)")
        rates[tag] = {
            "wall seconds": round(wall, 3),
            "table build seconds": float(build_secs) if build_secs else None,
            "rss after build kB": int(rss_build) if rss_build else None,
            "hwm after build kB": int(hwm_build) if hwm_build else None,
            "ticker samples": samples,
        }
        meds = sorted(s["steps/s"] for s in samples)
        med = meds[len(meds) // 2] if meds else 0
        rates[tag]["steps/s median"] = med
        print(f"[S] {tag}: build {rates[tag]['table build seconds']} s, "
              f"RSS after build {rates[tag]['rss after build kB']} kB "
              f"(hwm {rates[tag]['hwm after build kB']} kB), "
              f"samples {len(samples)}, median {med:.0f} steps/s")

    ram_med = rates["ram"]["steps/s median"]
    mmap_med = rates["mmap"]["steps/s median"]
    ratio = mmap_med / max(ram_med, 1)
    print(f"[S] mmap/ram median ratio: {ratio:.3f}")

    ok = (rec_ok and stride_r == str(STRIDE) and mmap_backend
          and cnt["merged rows"] == "0" and cnt["insert collisions"] == "0"
          and cnt["exact candidate checks"] == "1"
          and build_mults_r is not None and int(build_mults_r) < 100)
    print(f"\nCBT1 H=4M file-mode: {'ALL PASS' if ok else 'FAILED'}")

    artifact = {
        "table": TABLE.name,
        "table bytes": TABLE.stat().st_size,
        "recovery": {
            "target hex": hex(TARGET),
            "recovered hex": rec,
            "correct": rec_ok,
            "stride": stride_r,
            "mmap backend reported": mmap_backend,
            "counters": cnt,
            "table build seconds (mmap+validate)": build_secs_r,
            "build EC scalar mults": build_mults_r,
            "wall seconds": round(wall_r, 3),
            "pass": ok,
        },
        "paired scan (20M giants, far target)": {
            "scan size": SCAN_SIZE,
            "far target hex": hex(FAR_TARGET),
            "ram": rates["ram"],
            "mmap": rates["mmap"],
            "mmap/ram steps/s ratio": round(ratio, 4),
        },
    }
    tmp = ARTIFACT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(artifact, indent=1))
    tmp.rename(ARTIFACT)
    print(f"[artifact] {ARTIFACT}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
