#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Релизный long-cover прогон file-mode (CBT1, H = 4M) поверх --table-binary.

Гейт пройден: эквивалентность трёх источников (in-memory / file / resume)
бит-в-бит на 100k и 4M (equivalence_4m.json) - здесь тот же класс плана
прогоняется как релизный сценарий: kill/resume-цепочки, fingerprint-reject
чужой таблицей, полный miss-проход, метрики скорости и RSS.

Геометрия (H = 4M, стрид = row_count = 4M, тайлинг-ядро j в [1, H]
непрерывно => любой офсет chunk-локальной сетки достижим):

  START = 2^134; CHUNK = G*stride = 1e6 * 4e6 = 4e12; 12 chunk => 4.8e13.

Цели (offset от START; позиция гиганта m = chunk*G + g, d = offset - m*stride):

  A  chunk 0,  g=54321,   j=314159    - ранняя находка (1 сессия);
  B  chunk 3,  g=7777,    j=2345678   - kill-цепочка (cap 4 c);
  E  chunk 7,  partner: i=7G+313159, d = -lambda*c (c=271829) - орбита c,
     прямого d в [1,H] нет => находка только через partner-ветку
     (partner recoveries >= 1), kill-цепочка (cap 6 c);
  D  chunk 11, g=G-1,     j=H-1       - последняя точка покрытия,
     длинная kill-цепочка (cap 5 c);
  R  chunk 5,  g=999999,  j=1234567   - fingerprint: сессия kill (cap 4),
     затем --resume против ДРУГОЙ валидной CBT1 (100k) => exit 2
     "different plan", чекпоинт байт-в-бит не тронут, затем resume
     правильной таблицей до находки;
  M  вне покрытия (offset = cover + 12345) - негативный контроль: rc=1,
     chunks attempted == 12, гигант-шаги == (12+1)*G == 13,000,000 точно
     (cover-цикл + fallback исходного интервала), 0 хитов.

Протокол сессий: сценарий - фазовая машина (main -> [reject] -> main);
cap > 0 убивает сессию SIGKILL через cap секунд, следующая --resume.
Состояние атомарно пишется в JSON после каждой сессии - прогон
возобновляем между вызовами. Проверки:

  [1] восстановленный k == засеянному (mod n) во всех found-сценариях;
  [2] финальный checkpoint: chunks completed == ожидаемому, found == 1;
  [3] непрерывность цепочки: "chunks completed before this session"
      каждой resume-сессии == chunks completed чекпоинта предыдущей;
  [4] found == 0 во всех убитых сессиях; счётчики гигант-шагов в
      чекпоинтах монотонны; финальный cp counter giant steps == stdout;
  [5] R: reject-сессия rc == 2, "different plan" в stderr, digest
      чекпоинта до/после совпадает, затем обычная находка;
  [6] M: rc == 1, 0 recovery, 0 хитов, 13M гигант-шагов точно;
  [7] E: partner recoveries >= 1 в финальной сессии.
"""

import argparse
import hashlib
import json
import os
import re
import struct
import subprocess
import sys
import threading
import time
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
START = 1 << 134

H = 4_000_000
STRIDE = H
G = 1_000_000
CHUNK = G * STRIDE
COVER_CHUNKS = 12
COVER_SIZE = COVER_CHUNKS * CHUNK

TABLE = Path("../tables/table_cbt1_H4000000.cbt1")
CANON_100K = Path("../tables/representatives_canonical_H100000.txt")
WRONG_TABLE = Path("../tables/table_cbt1_H100000.cbt1")
HARD_SESSION_TIMEOUT = 300

COUNT_KEYS = (
    "giant steps", "giant steps scheduled", "giant steps planned",
    "invariant hits", "search hits", "orbit factors scheduled",
    "exact candidate checks", "partner recoveries", "merged rows",
    "insert collisions", "EC scalar multiplications",
    "total EC operations", "Lambda calculations",
)
TIMING_KEYS = (
    "table build seconds", "giant search seconds", "giant steps/second",
    "total seconds", "time to solution seconds",
    "cumulative elapsed seconds (checkpoint)",
)
TICKER_RE = re.compile(
    r"progress: elapsed ([0-9.]+)s \| giant steps (\d+) "
    r"\| steps/s ([0-9.]+) \| rss kB (\d+)")

failures = []


def pubkey_hex(scalar):
    numbers = ec.derive_private_key(
        scalar % N, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def build_targets():
    t = {
        "A": {"offset": 54321 * STRIDE + 314_159},
        "B": {"offset": 3 * CHUNK + 7777 * STRIDE + 2_345_678},
        "D": {"offset": COVER_SIZE - 1},
        "R": {"offset": 5 * CHUNK + 999_999 * STRIDE + 1_234_567},
        "M": {"offset": COVER_SIZE + 12_345},
    }
    i_giant, c = 7 * G + 313_159, 271_829
    s = (-LAMBDA * c) % N
    seed = START + i_giant * STRIDE + s
    assert seed < N, "partner seed wrapped mod n - pick another c"
    t["E"] = {"offset": i_giant * STRIDE + s,
              "partner": {"i giant": i_giant, "c": c}}
    for name, spec in t.items():
        spec["scalar"] = (START + spec["offset"]) % N
        assert spec["scalar"] == START + spec["offset"], name
        spec["pubkey"] = pubkey_hex(spec["scalar"])
        spec["chunk index"] = min(spec["offset"] // CHUNK, COVER_CHUNKS - 1)
    return t


SCENARIOS = {
    # name: (target, expect_found, expect_chunks, cap, checkpoint, reject)
    "A_early": ("A", True, 1, 0, True, False),
    "B_mid_kill": ("B", True, 4, 4, True, False),
    "E_partner": ("E", True, 8, 6, True, False),
    "D_last": ("D", True, 12, 5, True, False),
    "R_fingerprint_reject": ("R", True, 6, 4, True, True),
    "M_miss_full": ("M", False, COVER_CHUNKS, 0, False, False),
}


def save_json(out, results):
    tmp = out.with_suffix(".tmp")
    tmp.write_text(json.dumps(results, indent=1) + "\n")
    os.replace(tmp, out)


def peak_rss_kb(pid, sink, stop_event):
    status = Path(f"/proc/{pid}/status")
    while not stop_event.is_set():
        try:
            for line in status.read_text().splitlines():
                if line.startswith("VmHWM:"):
                    sink[0] = max(sink[0], int(line.split()[1]))
                    break
        except (FileNotFoundError, ProcessLookupError):
            return
        stop_event.wait(0.002)


def grab(pattern, text, cast=float):
    m = re.search(pattern, text)
    return cast(m.group(1)) if m else None


def run_session(binary, extra, cap, workdir, role, table_label):
    """Одна сессия движка; cap > 0 - принудительный SIGKILL через cap секунд."""
    cmd = [str(binary)] + extra
    sink = [0]
    stop_event = threading.Event()
    proc = subprocess.Popen(cmd, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, cwd=str(workdir))
    monitor = threading.Thread(target=peak_rss_kb,
                               args=(proc.pid, sink, stop_event), daemon=True)
    monitor.start()
    killed = False
    t0 = time.perf_counter()
    try:
        stdout, stderr = proc.communicate(
            timeout=cap if cap > 0 else HARD_SESSION_TIMEOUT)
    except subprocess.TimeoutExpired:
        killed = True
        proc.kill()
        stdout, stderr = proc.communicate()
    stop_event.set()
    monitor.join(timeout=0.2)

    session = {
        "role": role,
        "table": table_label,
        "killed": killed,
        "wall seconds": round(time.perf_counter() - t0, 3),
        "returncode": proc.returncode if not killed else None,
        "recovered": grab(r"recovered private key: 0x([0-9a-fA-F]{64})",
                          stdout, str),
        "already recovered": "already recovered" in stdout,
        "chunks attempted": grab(r"automatic chunks attempted: (\d+)",
                                 stdout, int),
        "resumed chunks before": grab(
            r"chunks completed before this session: (\d+)", stdout, int),
        "counters": {k: v for k in COUNT_KEYS
                     for v in [grab(rf"{re.escape(k)}: (\d+)", stdout, int)]
                     if v is not None},
        "timing": {k: v for k in TIMING_KEYS
                   for v in [grab(rf"{re.escape(k)}: ([0-9.eE+-]+)", stdout)]
                   if v is not None},
        "rss after table load kb": grab(
            r"RSS after table load: (\d+) kB", stdout, int),
        "ticker": [{"elapsed": float(e), "giant steps": int(s),
                    "steps/s": float(r), "rss kB": int(m)}
                   for e, s, r, m in TICKER_RE.findall(stdout)],
        "peak rss kb": sink[0],
    }
    if not killed and proc.returncode not in (0, 1):
        session["stderr tail"] = stderr[-500:]
    session["stderr tail"] = session.get("stderr tail") or stderr[-240:]
    return session


def parse_cp(path):
    if not path.exists():
        return None
    fields = {}
    for line in path.read_text().splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            fields[key.strip()] = value.strip()
    return fields


def cp_digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() \
        if path.exists() else None


def ensure_wrong_table(workdir, binary):
    """Валидная, но ДРУГАЯ CBT1 (100k) для fingerprint-reject сценария."""
    if WRONG_TABLE.exists():
        return False
    seed = START + 7 * 100_000 + 12_345
    chunk100k = 300_000 * 100_000
    cmd = [str(binary), "--target", pubkey_hex(seed),
           "--start", hex(START), "--size", str(10 * chunk100k),
           "--segment", str(10 * chunk100k), "--threads", "2",
           "--progress-interval", "0",
           "--table-scalars", str((workdir / CANON_100K).resolve()),
           "--glv-orbit", "--canonical-k3", "--partner-table",
           "--dump-table-binary", str((workdir / WRONG_TABLE).resolve())]
    proc = subprocess.run(cmd, text=True, capture_output=True,
                          cwd=str(workdir), timeout=120)
    assert proc.returncode == 0, proc.stderr[-400:]
    data = (workdir / WRONG_TABLE).read_bytes()[:256]
    assert data[:4] == b"CBT1"
    assert struct.unpack_from(">Q", data, 32)[0] == 100_000
    return True


def next_main_session(state, cap, cp_path, workdir, binary, args, base_extra):
    """Одна main-сессия; --resume только если чекпоинт уже на диске."""
    resume = cp_path is not None and cp_path.exists()
    extra = list(base_extra)
    if resume:
        extra += ["--resume"]
    return run_session(binary, extra, cap, workdir, "main", "correct")


def run_scenario(name, results, args, targets, workdir, binary):
    out = Path(args.out).resolve()
    tname, expect_found, expect_chunks, cap, use_cp, reject = SCENARIOS[name]
    spec = targets[tname]
    cp_path = workdir / f"cp_lc_{name}.txt" if use_cp else None

    base_extra = [
        "--target", spec["pubkey"],
        "--start", hex(START), "--size", str(CHUNK),
        "--segment", str(CHUNK), "--threads", str(args.threads),
        "--progress-interval", "5",
        "--cover-start", hex(START), "--cover-size", hex(COVER_SIZE),
        "--table-binary", str((workdir / TABLE).resolve()),
        "--glv-orbit", "--partner-table",
    ]
    if cp_path:
        base_extra += ["--checkpoint-file", str(cp_path),
                       "--checkpoint-steps", "1"]

    state = results["scenarios"].setdefault(name, {
        "target": tname, "expect found": expect_found,
        "expect chunks": expect_chunks, "cap seconds": cap,
        "phase": "main", "sessions": []})
    invoked = time.monotonic()
    while not state.get("done"):
        if time.monotonic() - invoked > args.budget:
            state["note"] = "budget reached, continue on next invocation"
            break
        if len(state["sessions"]) >= args.max_sessions:
            state["note"] = "max sessions reached"
            break
        phase = state["phase"]
        if phase == "main":
            session = next_main_session(
                state, cap, cp_path, workdir, binary, args, base_extra)
            state["sessions"].append(session)
            if not session["killed"]:
                state["phase"] = "reject" if reject else "done"
                if not reject:
                    state["done"] = True
        elif phase == "reject":
            session = {
                "role": "reject", "table": "wrong(100k)",
                "cp digest before": cp_digest(cp_path),
            }
            extra = list(base_extra) + [
                "--table-binary", str((workdir / WRONG_TABLE).resolve()),
                "--resume",
            ]
            # --table-binary дублируется: поздний аргумент перекрывает ранний
            session.update(run_session(binary, extra, 0, workdir,
                                       "reject", "wrong(100k)"))
            session["cp digest after"] = cp_digest(cp_path)
            state["sessions"].append(session)
            state["phase"] = "after_reject"
        elif phase == "after_reject":
            session = next_main_session(
                state, 0, cp_path, workdir, binary, args, base_extra)
            state["sessions"].append(session)
            if not session["killed"]:
                state["done"] = True
        if cp_path:
            cp = parse_cp(cp_path)
            if cp:
                state["sessions"][-1]["checkpoint after"] = {
                    "chunks completed": int(cp.get("chunks completed", -1)),
                    "found": cp.get("found"),
                    "giant steps": int(
                        cp.get("counter giant steps", -1)),
                    "exact checks": int(
                        cp.get("counter exact candidate checks", -1)),
                }
        save_json(out, results)


def verdict(name, results, targets):
    tname, expect_found, expect_chunks, cap, use_cp, reject = SCENARIOS[name]
    spec = targets[tname]
    state = results["scenarios"].get(name)
    if not state or not state["sessions"]:
        return None
    sessions = state["sessions"]
    bad = []
    main_sessions = [s for s in sessions if s["role"] == "main"]
    last = sessions[-1]
    final_cp = last.get("checkpoint after")

    if state.get("phase") not in ("done",) and not state.get("done"):
        bad.append(f"scenario incomplete (phase {state.get('phase')})")
    if last["killed"]:
        bad.append("last session was killed (scenario unfinished)")
        return {"verdict": "FAIL", "reasons": bad}

    if expect_found:
        recovered = last["recovered"]
        if recovered is None:
            bad.append("target not recovered")
        elif (int(recovered, 16) - spec["scalar"]) % N != 0:
            bad.append("recovered scalar mismatch")
        if last["returncode"] != 0:
            bad.append(f"rc={last['returncode']}")
        if not final_cp:
            bad.append("no checkpoint after final session")
        else:
            if final_cp["chunks completed"] != expect_chunks:
                bad.append(f"chunks completed "
                           f"{final_cp['chunks completed']} != {expect_chunks}")
            if final_cp["found"] != "1":
                bad.append("checkpoint found flag != 1")
            gs_stdout = last["counters"].get("giant steps")
            if gs_stdout is not None:
                if final_cp["giant steps"] != gs_stdout:
                    bad.append("final cp giant steps != stdout giant steps")
            elif not last.get("already recovered"):
                bad.append("final session lacks stdout counters")
            if last.get("already recovered") and not any(
                    s["counters"].get("giant steps") is not None
                    for s in main_sessions[:-1]):
                bad.append("vacuous scenario: only already-recovered resume")
        # chain continuity over main sessions
        for i in range(len(main_sessions) - 1):
            prev, nxt = main_sessions[i], main_sessions[i + 1]
            prev_cp = prev.get("checkpoint after")
            if not prev_cp or nxt.get("resumed chunks before") is None:
                continue
            if nxt["resumed chunks before"] != prev_cp["chunks completed"]:
                bad.append(
                    f"main session {i + 1} resumed at "
                    f"{nxt['resumed chunks before']}, checkpoint had "
                    f"{prev_cp['chunks completed']}")
            if prev_cp["found"] != "0":
                bad.append(f"main session {i}: found flag set before final")
        gs = [s["checkpoint after"]["giant steps"]
              for s in main_sessions if s.get("checkpoint after")
              and s["checkpoint after"]["giant steps"] >= 0]
        if gs != sorted(gs):
            bad.append("checkpoint giant steps counters not monotone")
        if cap > 0 and not any(s["killed"] for s in main_sessions):
            bad.append("kill/resume cycle never exercised")
        if name == "E_partner":
            pr = last["counters"].get("partner recoveries") or 0
            if pr < 1:
                bad.append(f"partner recoveries {pr} < 1")
    else:
        if last["recovered"] is not None:
            bad.append("unexpected recovery in negative control")
        if last["returncode"] != 1:
            bad.append(f"rc={last['returncode']} != 1")
        if last.get("chunks attempted") != expect_chunks:
            bad.append(f"chunks attempted {last.get('chunks attempted')} "
                       f"!= {expect_chunks}")
        steps = last["counters"].get("giant steps")
        if steps != (COVER_CHUNKS + 1) * G:
            bad.append(f"giant steps {steps} != {(COVER_CHUNKS + 1) * G}")
        for zero_key in ("search hits", "invariant hits",
                         "exact candidate checks", "partner recoveries"):
            if last["counters"].get(zero_key, 0) != 0:
                bad.append(f"{zero_key} != 0 in negative control")

    if reject:
        rej = [s for s in sessions if s["role"] == "reject"]
        if not rej:
            bad.append("reject session missing")
        else:
            r = rej[0]
            if r["returncode"] != 2:
                bad.append(f"reject rc={r['returncode']} != 2")
            if "different plan" not in r.get("stderr tail", ""):
                bad.append("reject stderr lacks 'different plan'")
            if r.get("cp digest before") != r.get("cp digest after"):
                bad.append("checkpoint modified by rejected resume")
            if r.get("recovered") is not None:
                bad.append("reject session recovered a key")

    return {"verdict": "PASS" if not bad else "FAIL",
            "reasons": bad,
            "sessions": len(sessions),
            "killed sessions": sum(1 for s in sessions if s["killed"]),
            "chunks completed": (final_cp or {}).get("chunks completed"),
            "total seconds": last["timing"].get("total seconds"),
            "table load seconds (last)": last["timing"].get(
                "table build seconds"),
            "giant steps/second": last["timing"].get("giant steps/second"),
            "tts seconds": last["timing"].get("time to solution seconds"),
            "peak rss kb": max(s["peak rss kb"] for s in sessions)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--out", default="../results/long_cover_cbt1_4m.json")
    parser.add_argument("--scenarios", default="",
                        help="comma list; default = all incomplete")
    parser.add_argument("--session-cap", type=int, default=None,
                        help="override kill cap (seconds) for all scenarios")
    parser.add_argument("--max-sessions", type=int, default=40)
    parser.add_argument("--budget", type=int, default=420,
                        help="wall budget per invocation, seconds")
    args = parser.parse_args()

    out = Path(args.out).resolve()
    workdir = Path(__file__).resolve().parent
    binary = (workdir / args.binary).resolve()

    if not (workdir / TABLE).resolve().exists():
        raise SystemExit(f"missing {TABLE}")
    if not (workdir / CANON_100K).resolve().exists():
        raise SystemExit(f"missing {CANON_100K}")

    results = {
        "heads": H, "stride": STRIDE, "giants per chunk": G,
        "chunk": CHUNK, "cover chunks": COVER_CHUNKS,
        "cover size": COVER_SIZE, "start": hex(START),
        "threads": args.threads, "mode": "file (CBT1 mmap)",
        "table": {"path": str(TABLE), "bytes": (workdir / TABLE).stat().st_size,
                  "sha256": sha256_of((workdir / TABLE).resolve())},
        "cache sidecar present at start":
            Path(str((workdir / TABLE).resolve()) + ".cache").exists(),
        "scenarios": {}}
    if out.exists():
        old = json.loads(out.read_text())
        if old.get("table", {}).get("sha256") == results["table"]["sha256"]:
            results["scenarios"] = old.get("scenarios", {})
        else:
            print("[setup] table changed - starting fresh results")

    if args.session_cap is not None:
        for key in SCENARIOS:
            t, f, c, _, cp, r = SCENARIOS[key]
            SCENARIOS[key] = (t, f, c, args.session_cap, cp, r)

    dumped = ensure_wrong_table(workdir, binary)
    if dumped:
        print(f"[setup] dumped {WRONG_TABLE} "
              f"({(workdir / WRONG_TABLE).stat().st_size} bytes)")
    results["wrong table"] = {
        "path": str(WRONG_TABLE), "bytes": (workdir / WRONG_TABLE).stat().st_size,
        "sha256": sha256_of((workdir / WRONG_TABLE).resolve())}

    targets = build_targets()
    results["targets"] = {k: {kk: (hex(vv) if kk == "scalar" and
                                   isinstance(vv, int) else vv)
                              for kk, vv in spec.items()}
                          for k, spec in targets.items()}

    wanted = ([s.strip() for s in args.scenarios.split(",") if s.strip()]
              or [n for n in SCENARIOS
                  if not results["scenarios"].get(n, {}).get("verdict")])
    for name in wanted:
        if name not in SCENARIOS:
            raise SystemExit(f"unknown scenario {name}")
        state = results["scenarios"].get(name, {})
        print(f"[scenario {name}] phase={state.get('phase', 'new')} "
              f"sessions so far: {len(state.get('sessions', []))}")
        run_scenario(name, results, args, targets, workdir, binary)
        v = verdict(name, results, targets)
        if v:
            results["scenarios"][name]["verdict"] = v
            print(f"[scenario {name}] {v['verdict']} "
                  f"({v['sessions']} sessions, "
                  f"{v['killed sessions']} killed) "
                  f"{'; '.join(v['reasons']) if v['reasons'] else ''}")
        save_json(out, results)

    done = [n for n in SCENARIOS
            if results["scenarios"].get(n, {}).get("verdict", {})
            .get("verdict") == "PASS"]
    failed = [n for n in SCENARIOS
              if results["scenarios"].get(n, {}).get("verdict", {})
              .get("verdict") == "FAIL"]
    print(f"\ncompleted scenarios: {len(done)}/{len(SCENARIOS)} PASS, "
          f"{len(failed)} FAIL")
    for name in SCENARIOS:
        v = results["scenarios"].get(name, {}).get("verdict")
        if v:
            print(f"  {name:22s} {v['verdict']:4s} sessions={v['sessions']} "
                  f"killed={v['killed sessions']} "
                  f"rss_peak={v['peak rss kb']}kB "
                  f"steps/s={v['giant steps/second']}")
    if failed:
        failures.extend(failed)
    print(f"[artifact] {out}")
    sys.exit(0 if not failures else 1)


if __name__ == "__main__":
    main()
