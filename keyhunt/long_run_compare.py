#!/usr/bin/env python3
"""Длинный прогон (а): canonical поверх representatives_canonical_H400000.txt.

Вопрос прогона: удерживает ли canonical 400k производительность
(~330k шагов/с), hit-экономику и корректность на многочанковом
покрытии с kill/resume, и сколько стоит checkpoint/resume
(пересборка таблицы на каждой сессии) в % от общего времени.

Дизайн:
- покрытие: 100 внешних chunk по 5e11 = 5e13 кандидатов от START = 2^134;
  chunk = --size = --segment (гигант-сетка чанко-локальна, 5e11 кратно
  stride 4e5, сетка выровнена между chunk);
- цель: одна, глубоко в покрытии (chunk 85, гигант-индекс 625000,
  c = 12345 - канонический представитель головы 12345); позиция
  зеркалится фолдом при необходимости; seed 20261104 фиксирует дизайн;
- сессии: движок убивается SIGKILL по --session-cap секунд, --resume
  продолжает с checkpoint (пересборка таблицы = цена checkpoint);
- steps/s каждые 60s: --progress-interval 60 движка, строки
  "progress: ..." собираются по сессиям;
- JSON атомарен (tmp+rename) после каждой сессии, дозапускResume-safe.

Запуск:  python3 long_run_compare.py [--session-cap 120] [--budget 540]
Выход:   ../long_run_canonical.json
"""

import argparse
import json
import os
import re
import signal
import subprocess
import threading
import time
from pathlib import Path

START = 1 << 134
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
CHUNK = 500_000_000_000
CHUNKS = 100
STRIDE = 400_000                 # строки в canonical-файле H400000
TARGET_CHUNK = 85
TARGET_GIANT_INDEX = 625_000
TARGET_RESIDUAL = 12_345         # c: канонический представитель головы 12345
CANON_FILE = Path("../representatives_canonical_H400000.txt")
SEED = 20261104

COUNTER_KEYS = (
    "giant steps", "giant steps scheduled", "giant steps planned",
    "invariant hits", "search hits", "build hits",
    "exact candidate checks", "orbit factors scheduled",
    "EC scalar multiplications", "EC scalar multiplications (build)",
    "EC point additions", "total EC operations",
    "canonical rows merged", "orbit duplicate rows skipped",
    "partner recoveries", "raw baby rows built", "unique baby rows built",
    "duplicate baby keys",
)
TIMING_KEYS = ("table build seconds", "giant search seconds",
               "giant steps/second", "total seconds",
               "time to solution seconds",
               "cumulative elapsed seconds (checkpoint)")

PROGRESS_RE = re.compile(
    r"progress: elapsed ([0-9.]+)s \| giant steps (\d+) \| "
    r"steps/s ([0-9.]+) \| rss kB (-?\d+)")


def pubkey_hex(scalar):
    from cryptography.hazmat.primitives.asymmetric import ec
    numbers = ec.derive_private_key(scalar % N, ec.SECP256K1()) \
        .public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def expected_scalar():
    offset = (TARGET_CHUNK * CHUNK + TARGET_GIANT_INDEX * STRIDE
              + TARGET_RESIDUAL)
    return START + offset


def run_session(binary, cmd_extra, cap, workdir):
    """Одна сессия: запуск до cap секунд, затем SIGKILL. stdout стримится."""
    proc = subprocess.Popen([str(binary)] + cmd_extra, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            cwd=workdir)
    lines = []
    err_chunks = []

    def drain_err():
        for line in proc.stderr:
            err_chunks.append(line)

    def kill():
        try:
            proc.send_signal(signal.SIGKILL)
        except ProcessLookupError:
            pass

    err_thread = threading.Thread(target=drain_err, daemon=True)
    err_thread.start()
    timer = threading.Timer(cap, kill)
    timer.start()
    killed = False
    try:
        for line in proc.stdout:
            lines.append(line.rstrip("\n"))
        proc.wait()
        killed = proc.returncode == -signal.SIGKILL
    finally:
        timer.cancel()
    err_thread.join(timeout=2)
    err = "".join(err_chunks)
    out = "\n".join(lines)

    counters = {}
    for name in COUNTER_KEYS:
        m = re.search(rf"{re.escape(name)}: (\d+)", out)
        if m:
            counters[name] = int(m.group(1))
    timing = {}
    for name in TIMING_KEYS:
        m = re.search(rf"{re.escape(name)}: ([0-9.eE+-]+)", out)
        if m:
            timing[name] = float(m.group(1))
    progress = [{"elapsed s": float(e), "giant steps": int(g),
                 "steps/s": float(s), "rss kb": int(r)}
                for e, g, s, r in PROGRESS_RE.findall(out)]
    stride_line = re.search(r"giant stride \(post-merge rows\): (\d+)", out)
    recovered = re.search(r"recovered private key: 0x([0-9a-fA-F]+)", out)
    return {
        "cmd": " ".join(cmd_extra),
        "killed": killed,
        "returncode": proc.returncode,
        "recovered": int(recovered.group(1), 16) if recovered else None,
        "stride line": int(stride_line.group(1)) if stride_line else None,
        "counters": counters,
        "timing": timing,
        "progress samples": progress,
        "stderr tail": err.strip().splitlines()[-3:],
        "peak rss kb": max((p["rss kb"] for p in progress), default=0),
    }


def save_json(path, results):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(results, indent=2) + "\n")
    os.replace(tmp, path)


def verdict(results, expect_scalar):
    """Контрольные проверки длинного прогона (после финальной сессии)."""
    bad = []
    sessions = results["sessions"]
    last = sessions[-1]
    if last["recovered"] is None or sessions[-1]["killed"]:
        return {"verdict": "INCOMPLETE", "reasons": ["run not finished"]}
    if last["recovered"] % N != expect_scalar % N:
        bad.append("recovered scalar != seeded target")
    strides = {s["stride line"] for s in sessions if s["stride line"]}
    if strides != {STRIDE}:
        bad.append(f"giant stride lines {strides} != {{{STRIDE}}}")
    if last["counters"].get("build hits", 0) != 0:
        bad.append("build hits != 0 on a canonical file")
    if last["counters"].get("search hits", 0) < 1:
        bad.append("search hits < 1 at TTS")
    scans = [s for s in sessions if s["timing"].get("giant search seconds", 0) > 60]
    for i, s in enumerate(scans):
        if len(s["progress samples"]) < 1:
            bad.append(f"session {i}: no 60s progress sample despite "
                       "search > 60s")
    rates = [p["steps/s"] for s in sessions for p in s["progress samples"]]
    if rates and max(rates) / max(min(rates), 1) > 1.5:
        bad.append(f"steps/s unstable: min {min(rates):.0f}, "
                   f"max {max(rates):.0f}")
    killed = sum(1 for s in sessions if s["killed"])
    if killed < 1:
        bad.append("kill/resume cycle never exercised")
    builds = [s["timing"].get("table build seconds", 0.0) for s in sessions]
    # Killed sessions lose the engine's tail timing (printed only at normal
    # exit) - use the harness-measured wall seconds as the denominator.
    totals = [s.get("wall seconds",
                    s["timing"].get("total seconds", 0.0)) for s in sessions]
    rebuild = sum(builds[1:])
    total = sum(totals)
    results["checkpoint cost"] = {
        "sessions": len(sessions),
        "killed": killed,
        "rebuild seconds (sessions 2+)": round(rebuild, 3),
        "total seconds (all sessions)": round(total, 3),
        "checkpoint cost pct": round(100 * rebuild / total, 2) if total else None,
    }
    results["steps/s samples"] = rates
    return {"verdict": "PASS" if not bad else "FAIL", "reasons": bad}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--out", default="../long_run_canonical.json")
    parser.add_argument("--session-cap", type=float, default=120.0,
                        help="SIGKILL after this many seconds per session")
    parser.add_argument("--max-sessions", type=int, default=8)
    parser.add_argument("--budget", type=float, default=480.0,
                        help="wall budget for this invocation, seconds")
    args = parser.parse_args()

    out = Path(args.out).resolve()
    workdir = Path(__file__).resolve().parent
    binary = (workdir / args.binary).resolve()
    if not (workdir / CANON_FILE).resolve().exists():
        raise SystemExit(f"missing {CANON_FILE}")

    expect_scalar = expected_scalar()
    target_pub = pubkey_hex(expect_scalar)
    cover_size = CHUNKS * CHUNK

    base_cmd = ["--target", target_pub, "--start", hex(START),
                "--size", str(CHUNK), "--segment", str(CHUNK),
                "--threads", str(args.threads),
                "--table-scalars", str((workdir / CANON_FILE).resolve()),
                "--glv-orbit", "--phase-signature", "--partner-table",
                "--cover-start", hex(START), "--cover-size", hex(cover_size),
                "--checkpoint-file", str(workdir / "cp_long_run.txt"),
                "--checkpoint-steps", "1",
                "--progress-interval", "60"]

    results = None
    if out.exists():
        try:
            loaded = json.loads(out.read_text())
            if loaded.get("design", {}).get("seed") == SEED:
                results = loaded
        except (json.JSONDecodeError, OSError):
            results = None
    fresh = results is None
    if fresh:
        results = {
            "design": {
                "seed": SEED,
                "question": "canonical 400k long cover: steps/s stability, "
                            "hit economy, checkpoint cost pct",
                "start": hex(START),
                "chunk": CHUNK, "chunks": CHUNKS,
                "cover size": cover_size,
                "stride": STRIDE,
                "target chunk": TARGET_CHUNK,
                "target giant index": TARGET_GIANT_INDEX,
                "target residual c": TARGET_RESIDUAL,
                "target scalar": hex(expect_scalar),
                "target pubkey": target_pub,
                "session cap s": args.session_cap,
            },
            "sessions": [],
        }
    results.setdefault("sessions", [])

    invoked = time.monotonic()
    while True:
        sessions = results["sessions"]
        if sessions and sessions[-1]["killed"] is False:
            break  # движок завершился сам - прогон исполнен
        if sessions and time.monotonic() - invoked > args.budget:
            results["note"] = "budget reached, continue on next invocation"
            break
        if len(sessions) >= args.max_sessions:
            results["note"] = "max sessions reached"
            break
        cmd_extra = list(base_cmd)
        if sessions:
            cmd_extra.append("--resume")
        s_start = time.monotonic()
        session = run_session(binary, cmd_extra, args.session_cap, workdir)
        session["wall seconds"] = round(time.monotonic() - s_start, 3)
        sessions.append(session)
        save_json(out, results)
        tag = "KILL" if session["killed"] else (
            "FOUND" if session["recovered"] is not None else "EXIT")
        print(f"[session {len(sessions)}] {tag}: "
              f"build {session['timing'].get('table build seconds', 0):.1f}s, "
              f"search {session['timing'].get('giant search seconds', 0):.1f}s, "
              f"samples {len(session['progress samples'])}, "
              f"search hits {session['counters'].get('search hits')}, "
              f"build hits {session['counters'].get('build hits')}",
              flush=True)

    sessions = results["sessions"]
    if sessions and sessions[-1]["killed"] is False:
        results["verdict block"] = verdict(results, expect_scalar)
        save_json(out, results)
        v = results["verdict block"]
        print(f"verdict: {v['verdict']} {v.get('reasons', [])}")
        if "checkpoint cost" in results:
            cc = results["checkpoint cost"]
            print(f"checkpoint cost: {cc['rebuild seconds (sessions 2+)']}s "
                  f"rebuild / {cc['total seconds (all sessions)']}s total = "
                  f"{cc['checkpoint cost pct']}%")
            print(f"steps/s samples: min {min(results['steps/s samples']):.0f}, "
                  f"max {max(results['steps/s samples']):.0f}, "
                  f"n={len(results['steps/s samples'])}")
    else:
        save_json(out, results)
        print(f"sessions so far: {len(sessions)}, "
              f"continue with the same command (resume-safe)")


if __name__ == "__main__":
    main()
