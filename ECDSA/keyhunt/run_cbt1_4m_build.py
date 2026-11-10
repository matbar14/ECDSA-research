#!/usr/bin/env python3
"""CBT1 H=4M: build + dump из expanded-файла (merge-режим), фазовый RSS.

Запускает eng_bsgs в merge-режиме (--table-scalars expanded + --canonical-k3)
с --dump-table-binary: тот же путь, что реальная сессия текстового формата,
плюс сериализация поискового состояния в CBT1. Внешний сэмплер 0.25 c
дополняет фазовый таймлайн движка (VmRSS/VmHWM на границах load/pre-pass/
build) непрерывной кривой - в том числе на фазе dump, которую таймлайн
движка не покрывает. Артефакт пишется атомарно (resume-безопасность).
"""

import hashlib
import json
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent      # ECDSA-research/
ENGINE = ROOT / "keyhunt" / "eng_bsgs"
EXPANDED = ROOT / "tables" / "representatives_expanded_H4000000.txt"
DUMP = ROOT / "tables" / "table_cbt1_H4000000.cbt1"
ARTIFACT = ROOT / "results" / "cbt1_4m_build.json"

START = 0x4000000000000000000000000000000000        # 2^134, как в h_curve
STRIDE = 4_000_000
TARGET = START + 2 * STRIDE + 777                   # контрольная цель [1]
SIZE = 500_000_000_000                              # гигант-шагов: 125k

def compressed_pubkey_hex(scalar: int) -> str:
    from cryptography.hazmat.primitives.asymmetric import ec
    numbers = ec.derive_private_key(scalar, ec.SECP256K1()) \
        .public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


CMD = [
    str(ENGINE),
    "--target", compressed_pubkey_hex(TARGET),
    "--start", hex(START), "--size", str(SIZE),
    "--segment", str(SIZE), "--threads", "2", "--progress-interval", "0",
    "--glv-orbit", "--canonical-k3", "--partner-table",
    "--table-scalars", str(EXPANDED),
    "--dump-table-binary", str(DUMP),
]


def sample_rss(proc, stop_evt, samples):
    tick = time.monotonic()
    while not stop_evt.is_set():
        try:
            with open(f"/proc/{proc.pid}/status") as f:
                rss = hwm = None
                for line in f:
                    if line.startswith("VmRSS:"):
                        rss = int(line.split()[1])
                    elif line.startswith("VmHWM:"):
                        hwm = int(line.split()[1])
            if rss is not None:
                now = time.monotonic() - tick
                samples.append([round(now, 3), rss, hwm])
        except (FileNotFoundError, ProcessLookupError):
            return
        stop_evt.wait(0.25)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def grab(out, pattern):
    m = re.search(pattern, out)
    return m.group(1) if m else None


def main():
    if ARTIFACT.exists():
        print(f"[skip] artifact exists: {ARTIFACT}")
        return
    if DUMP.exists():
        raise SystemExit(f"[E] dump file already exists, refuse to rerun: {DUMP}")
    if not EXPANDED.exists():
        raise SystemExit(f"[E] missing input: {EXPANDED}")

    print(f"[run] {' '.join(CMD)}")
    t0 = time.monotonic()
    proc = subprocess.Popen(CMD, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    samples = []
    stop_evt = threading.Event()
    sampler = threading.Thread(target=sample_rss,
                               args=(proc, stop_evt, samples), daemon=True)
    sampler.start()
    try:
        out, err = proc.communicate(timeout=540)
    except subprocess.TimeoutExpired:
        proc.kill()
        out, err = proc.communicate()
        stop_evt.set()
        raise SystemExit("[E] engine timeout 540 s")
    wall = time.monotonic() - t0
    stop_evt.set()
    sampler.join(timeout=1)

    if proc.returncode != 0:
        print(out[-2000:] if out else "", err[-1000:] if err else "",
              sep="\n")
        raise SystemExit(f"[E] engine exit {proc.returncode}")

    rec = grab(out, r"recovered private key: 0x([0-9a-f]+)")
    if rec is None or int(rec, 16) != TARGET:
        raise SystemExit(f"[E] control target not recovered (rec={rec})")

    dump_bytes = DUMP.stat().st_size
    digest = sha256_file(DUMP)

    # Downsample: keep every 4th sample (1 s resolution) for the artifact.
    kept = samples[::4] if len(samples) > 4000 else samples
    peak_seen = max((s[2] for s in samples), default=None)

    artifact = {
        "mode": "expanded->merge build + dump-table-binary (CBT1)",
        "heads": STRIDE,
        "input": {
            "file": EXPANDED.name,
            "bytes": EXPANDED.stat().st_size,
            "rows": 12_000_000,
            "sha256": json.loads((ROOT / "tables" / "representatives_files.json")
                                 .read_text())["expanded_files"][-1]["sha256"],
        },
        "cmd": CMD,
        "returncode": proc.returncode,
        "wall seconds": round(wall, 3),
        "recovered scalar hex": rec,
        "dump file": {
            "path": DUMP.name,
            "bytes": dump_bytes,
            "expected bytes": 256 + 65 * STRIDE,
            "sha256": digest,
        },
        "engine": {
            "table rows stored": grab(out, r"table rows stored: (\d+)"),
            "build merged rows": grab(out, r"build merged rows (\d+)"),
            "giant stride": grab(out, r"giant stride \(post-merge rows\): (\d+)"),
            "table build seconds": grab(out, r"table build seconds: ([0-9.]+)"),
            "rss after load kB": grab(out, r"RSS after table load: (\d+)"),
            "hwm after load kB": grab(out, r"RSS after table load: \d+ kB \(peak (\d+)"),
            "rss after prepass kB": grab(out, r"RSS after canonical pre-pass: (\d+)"),
            "hwm after prepass kB": grab(out, r"pre-pass: \d+ kB \(peak (\d+)"),
            "rss after build kB": grab(out, r"RSS after table build: (\d+)"),
            "hwm after build kB": grab(out, r"table build: \d+ kB \(peak (\d+)"),
            "ec scalar mults": grab(out, r"EC scalar multiplications: (\d+)"),
            "giant steps": grab(out, r"giant steps: (\d+)"),
        },
        "sampler": {
            "interval seconds": 0.25,
            "samples": len(samples),
            "peak VmHWM kB seen": peak_seen,
            "series (t s, VmRSS kB, VmHWM kB), 1 s resolution": kept,
        },
    }
    tmp = ARTIFACT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(artifact, indent=1))
    tmp.rename(ARTIFACT)

    e = artifact["engine"]
    print(f"[ok] wall {wall:.1f}s, build {e['table build seconds']}s, "
          f"rows {e['table rows stored']} (merged in {e['build merged rows']}), "
          f"stride {e['giant stride']}")
    print(f"[rss] load {e['rss after load kB']} kB (hwm {e['hwm after load kB']}), "
          f"pre-pass {e['rss after prepass kB']} kB (hwm {e['hwm after prepass kB']}), "
          f"build {e['rss after build kB']} kB (hwm {e['hwm after build kB']}), "
          f"sampler peak {peak_seen} kB")
    print(f"[dump] {dump_bytes} B (expected {256 + 65 * STRIDE}), "
          f"sha256 {digest[:16]}...")
    print(f"[artifact] {ARTIFACT}")


if __name__ == "__main__":
    main()
