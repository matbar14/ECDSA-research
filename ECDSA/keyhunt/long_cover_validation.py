#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Длинная валидация --cover поверх файла представителей (PROGRESS_MEMO.md,
«следующий шаг» после checkpoint/resume): автоматическое покрытие
--cover-start/--cover-size с контрольными целями, kill/resume-циклами и
сравнением обоих режимов таблицы на длинном масштабе.

Конфигурация (H = 400000, mirror H-кривой):

  cover      16 внешних chunk по 5e11 = 8e12 кандидатов от START = 2^134;
  canonical  --table-scalars representatives_canonical_H400000.txt
             (H = 400k строк, стрид H) + --glv-orbit --phase-signature
             --partner-table - производственная конфигурация полного
             покрытия: строки файла суть 1..H, поэтому тайлинг
             i*stride + [1, H] непрерывен и ЛЮБОЙ офсет достижим;
  expanded   --table-scalars representatives_expanded_H400000.txt
             (3H = 1.2M строк: j, lambda*j, lambda^2*j) + --canonical-k3
             (merge в 400k) + те же флаги. ВАЖНО (выяснено валидацией):
             pre-pass канонизации стоит ДО m = table_scalars.size(),
             поэтому стрид = ЧИСЛО СТРОК ПОСЛЕ merge = 400k - тот же
             класс покрытия/скорости, что canonical (гигант-шаги и окно
             совпадают; разница только в сборке: вход 1.2M + merge);
  expanded_plain  тот же файл БЕЗ --canonical-k3: стрид = 1.2M строк
             входа, ядро покрытия [1, H] после каждой гигант-позиции,
             офсеты из (H, 3H) в chunk-локальной сетке недостижимы.

Отсюда дизайн целей (offset от START). ВАЖНО: гигант-сетка чанко-локальна
(позиции chunk_start + g*stride, а chunk = 5e11 не кратен strides 1.2e6),
поэтому остатки целей считаются от начала СВОЕГО chunk:

  A  249999612345  chunk 0, 5e11 mod обоих стридов сетка = глобальная:
                   1.2e6*208333 + 12345, остаток 12345 при ОБОИХ
                   стридах (1.2e6 = 3*4e5 кратно) - находка обоих режимов;
  B  500001207777  chunk 1, delta = 5e11 + 1.2e6 + 7777: остаток 7777
                   в chunk-локальной сетке ОБОИХ стридов - оба режима;
  C  7501000123456                            остаток canonical/merge
                   123456 - находка обоих полных режимов; в plain
                   expanded остаток 523456 > H - ПРОПУСК (негативный
                   контроль структуры покрытия);
  D  8e12 - 1      последняя точка покрытия, остаток canonical/merge
                   399999 <= H - находка; plain expanded 799999 > H -
                   пропуск;
  E  partner-цель  START + 11e6*H - lambda*c (c = 313159): canonical и
                   merge через --partner-table, plain expanded - пропуск.

Протокол сессий: каждый сценарий - цикл запусков движка; сессии с
cap > 0 принудительно убиваются SIGKILL через cap секунд (среда убивает
долгие процессы - checkpoint ровно для этого), каждая следующая сессия
--resume. Состояние атомарно пишется в JSON после каждой сессии, поэтому
валидация возобновляема между вызовами. Проверки на сценарий:

  [1] восстановленный ключ == засеянный (mod n);
  [2] финальный checkpoint: chunks completed == ожидаемому числу chunk -
      каждый chunk учтён ровно один раз (нет дыр и двойного счёта);
  [3] цепочка непрерывности: "chunks completed before this session"
      каждой resume-сессии == checkpoint предыдущей сессии;
  [4] found-флаг в checkpoint: 0 во всех убитых сессиях, 1 в финальной;
  [5] у сценариев с cap > 0 есть хотя бы одна убитая сессия;
  [6] кумулятивные счётчики (giant steps по checkpoint) монотонны;
  [7] негативные контроли expanded: rc=1, восстановлений нет, cover
      пройден целиком (attempted == 16).
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
import threading
import time
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
START = 1 << 134

H = 400_000
CHUNK = 500_000_000_000
COVER_CHUNKS = 16
CANON_FILE = Path("../tables/representatives_canonical_H400000.txt")
EXPANDED_FILE = Path("../tables/representatives_expanded_H400000.txt")

COUNT_KEYS = (
    "giant steps", "invariant hits", "exact candidate checks",
    "EC scalar multiplications", "EC point additions",
    "EC scalar multiplications (build)",
    "EC scalar multiplications (giant search, cumulative)",
    "EC scalar multiplications (exact candidate checks)",
    "total EC operations", "partner recoveries",
    "signature pair collapses", "signature pair mismatches",
    "signature anomalies (6-branch fallback)",
)


def pubkey_hex(scalar):
    numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def generate_expanded(path):
    """3H строк: j, lambda*j, lambda^2*j для голов j = 1..H.

    Детерминирована головами; уникальность всех 3H скаляров проверяется.
    """
    rows = []
    seen = set()
    for j in range(1, H + 1):
        for v in (j, j * LAMBDA % N, j * LAMBDA % N * LAMBDA % N):
            if v in seen:
                raise SystemExit(f"duplicate expanded scalar at head {j}")
            seen.add(v)
            rows.append(v)
    lines = [f"# expanded representatives, heads 1..{H}, 3 rows per head",
             f"# rows: {len(rows)} (stride = rows = {3 * H})"]
    lines += [hex(v) for v in rows]
    path.write_text("\n".join(lines) + "\n")
    return len(rows)


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def build_targets():
    t = {
        "A": {"offset": 1_200_000 * 208_333 + 12_345},
        "B": {"offset": CHUNK + 1_200_000 + 7_777},
        "C": {"offset": 7_500_000_000_000 + 1_000_000_000 + 123_456},
        "D": {"offset": COVER_CHUNKS * CHUNK - 1},
    }
    i_giant, c = 11_000_000, 313_159
    s = (-pow(LAMBDA, 1, N) * c) % N
    t["E"] = {"offset": i_giant * H + s, "partner": (i_giant, c, 1, -1)}
    for name, spec in t.items():
        spec["scalar"] = (START + spec["offset"]) % N
        spec["pubkey"] = pubkey_hex(spec["scalar"])
        spec["chunk index"] = min(spec["offset"] // CHUNK, COVER_CHUNKS - 1)
    return t


SCENARIOS = {
    "canonical_A": ("canonical", "A", True, 1, 0),
    "canonical_B": ("canonical", "B", True, 2, 0),
    "canonical_E": ("canonical", "E", True, 9, 0),
    "canonical_C": ("canonical", "C", True, 16, 12),
    "canonical_D": ("canonical", "D", True, 16, 12),
    "expanded_A": ("expanded", "A", True, 1, 0),
    "expanded_B": ("expanded", "B", True, 2, 0),
    "expanded_C": ("expanded", "C", True, 16, 0),
    "expanded_D": ("expanded", "D", True, 16, 0),
    "expanded_E": ("expanded", "E", True, 9, 0),
    "expanded_plain_C": ("expanded_plain", "C", False, 16, 0),
    "expanded_plain_D": ("expanded_plain", "D", False, 16, 0),
    "expanded_plain_E": ("expanded_plain", "E", False, 16, 0),
}


def save_json(out, results):
    tmp = out.with_suffix(".tmp")
    tmp.write_text(json.dumps(results, indent=2) + "\n")
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


def run_session(binary, extra, cap, workdir):
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
        stdout, stderr = proc.communicate(timeout=cap or None)
    except subprocess.TimeoutExpired:
        killed = True
        proc.kill()
        stdout, stderr = proc.communicate()
    stop_event.set()
    monitor.join(timeout=0.2)

    m = re.search(r"recovered private key: 0x([0-9a-fA-F]{64})", stdout)
    recovered = int(m.group(1), 16) if m else None
    counters = {k: int(v) for k in COUNT_KEYS
                for v in [next(iter(re.findall(
                    rf"{re.escape(k)}: (\d+)", stdout)), -1)] if int(v) >= 0}
    timing = {k: float(v) for k in ("table build seconds",
                                    "giant search seconds",
                                    "giant steps/second", "total seconds",
                                    "time to solution seconds")
              for v in [next(iter(re.findall(
                  rf"{re.escape(k)}: ([0-9.eE+-]+)", stdout)), "")] if v}
    chunks_att = re.search(r"automatic chunks attempted: (\d+)", stdout)
    resumed_before = re.search(
        r"chunks completed before this session: (\d+)", stdout)
    session = {
        "killed": killed,
        "wall seconds": round(time.perf_counter() - t0, 3),
        "returncode": proc.returncode if not killed else None,
        "recovered": hex(recovered) if recovered is not None else None,
        "chunks attempted": int(chunks_att.group(1)) if chunks_att else None,
        "resumed chunks before": int(resumed_before.group(1))
        if resumed_before else None,
        "counters": counters,
        "timing": timing,
        "peak rss kb": sink[0],
    }
    if not killed and proc.returncode not in (0, 1):
        session["stderr tail"] = stderr[-400:]
    return session


def parse_checkpoint(path):
    if not path.exists():
        return None
    fields = {}
    for line in path.read_text().splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            fields[key.strip()] = value.strip()
    return fields


def run_scenario(name, results, args, targets, workdir, binary):
    out = Path(args.out).resolve()
    mode, tname, expect_found, expect_chunks, cap = SCENARIOS[name]
    spec = targets[tname]
    cp_path = workdir / f"cp_{name}.txt" if expect_found else None

    if mode == "canonical":
        extra = ["--table-scalars", str(CANON_FILE),
                 "--glv-orbit", "--phase-signature", "--partner-table"]
    elif mode == "expanded":
        extra = ["--table-scalars", str(EXPANDED_FILE),
                 "--glv-orbit", "--phase-signature", "--partner-table",
                 "--canonical-k3"]
    else:  # expanded_plain: стрид 1.2M, дыры покрытия (H, 3H)
        extra = ["--table-scalars", str(EXPANDED_FILE),
                 "--glv-orbit", "--phase-signature", "--partner-table"]
    extra += ["--target", spec["pubkey"], "--start", hex(START),
              "--size", str(CHUNK), "--segment", str(CHUNK),
              "--threads", str(args.threads),
              "--cover-start", hex(START),
              "--cover-size", hex(COVER_CHUNKS * CHUNK)]
    if cp_path:
        extra += ["--checkpoint-file", str(cp_path), "--checkpoint-steps", "1"]

    state = results["scenarios"].setdefault(
        name, {"mode": mode, "target": tname, "expect found": expect_found,
               "expect chunks": expect_chunks, "sessions": []})
    invoked = time.monotonic()
    while True:
        sessions = state["sessions"]
        if sessions and sessions[-1]["killed"] is False:
            break  # движок завершился сам - сценарий исполнен
        if sessions and time.monotonic() - invoked > args.budget:
            state["note"] = "budget reached, continue on next invocation"
            break
        if len(sessions) >= args.max_sessions:
            state["note"] = "max sessions reached"
            break
        cmd_extra = list(extra)
        if sessions:
            cmd_extra.append("--resume")
        session = run_session(binary, cmd_extra, cap, workdir)
        sessions.append(session)
        if cp_path:
            cp = parse_checkpoint(cp_path)
            if cp:
                session["checkpoint after"] = {
                    "chunks completed": int(cp.get("chunks completed", -1)),
                    "found": cp.get("found"),
                    "giant steps": int(cp.get("counter giant steps", -1)),
                    "exact checks": int(
                        cp.get("counter exact candidate checks", -1)),
                }
        save_json(out, results)


def verdict(name, results, targets):
    mode, tname, expect_found, expect_chunks, cap = SCENARIOS[name]
    spec = targets[tname]
    state = results["scenarios"].get(name)
    if not state or not state["sessions"]:
        return None
    sessions = state["sessions"]
    bad = []
    last = sessions[-1]
    final_cp = last.get("checkpoint after")

    if last["killed"]:
        bad.append("last session was killed (scenario unfinished)")
        return {"verdict": "FAIL", "reasons": bad}

    if expect_found:
        recovered = last["recovered"]
        if recovered is None:
            bad.append("target not recovered")
        else:
            k = spec["scalar"]
            if (int(recovered, 16) - k) % N != 0:
                bad.append("recovered scalar mismatch")
        if last["returncode"] != 0:
            bad.append(f"rc={last['returncode']}")
        if not final_cp:
            bad.append("no checkpoint after final session")
        else:
            done = final_cp["chunks completed"]
            if done != expect_chunks:
                bad.append(f"chunks completed {done} != {expect_chunks}")
            if final_cp["found"] != "1":
                bad.append("checkpoint found flag != 1")
        chain = [s.get("checkpoint after") for s in sessions[:-1]]
        for i, cp in enumerate(chain):
            if not cp:
                continue
            if cp["found"] != "0":
                bad.append(f"session {i}: found flag set before final")
            nxt = sessions[i + 1].get("resumed chunks before")
            if nxt is not None and nxt != cp["chunks completed"]:
                bad.append(
                    f"session {i + 1} resumed at {nxt}, "
                    f"checkpoint had {cp['chunks completed']}")
        gs = [cp["giant steps"] for cp in chain
              if cp and cp["giant steps"] >= 0]
        if gs != sorted(gs):
            bad.append("giant steps counters not monotone")
        if cap > 0 and not any(s["killed"] for s in sessions):
            bad.append("kill/resume cycle never exercised")
    else:
        if last["recovered"] is not None:
            bad.append("unexpected recovery in negative control")
        if last["returncode"] != 1:
            bad.append(f"rc={last['returncode']} != 1")
        att = last.get("chunks attempted")
        if att != expect_chunks:
            bad.append(f"chunks attempted {att} != {expect_chunks}")

    return {"verdict": "PASS" if not bad else "FAIL",
            "reasons": bad,
            "sessions": len(sessions),
            "killed sessions": sum(1 for s in sessions if s["killed"]),
            "chunks completed": (final_cp or {}).get("chunks completed"),
            "total seconds": last["timing"].get("total seconds"),
            "build seconds": last["timing"].get("table build seconds"),
            "giant steps/second": last["timing"].get("giant steps/second"),
            "tts seconds": last["timing"].get("time to solution seconds"),
            "checks per hit": last["counters"].get("exact candidate checks"),
            "peak rss kb": max(s["peak rss kb"] for s in sessions)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--out", default="../results/long_cover_validation.json")
    parser.add_argument("--scenarios", default="",
                        help="comma list; default = all incomplete")
    parser.add_argument("--session-cap", type=int, default=None,
                        help="override kill cap (seconds) for all scenarios")
    parser.add_argument("--max-sessions", type=int, default=24)
    parser.add_argument("--budget", type=int, default=480,
                        help="wall budget per invocation, seconds")
    args = parser.parse_args()

    out = Path(args.out).resolve()
    workdir = Path(__file__).resolve().parent
    binary = (workdir / args.binary).resolve()

    if not CANON_FILE.resolve().exists():
        raise SystemExit(f"missing {CANON_FILE}")
    expanded_meta = {}
    if not EXPANDED_FILE.exists():
        print(f"[setup] generating {EXPANDED_FILE} ...")
        rows = generate_expanded(EXPANDED_FILE)
        print(f"[setup] {rows} rows written")
    expanded_meta = {"rows": 3 * H, "sha256": sha256_of(EXPANDED_FILE),
                     "stride": 3 * H}

    results = {"heads": H, "chunk": CHUNK, "cover chunks": COVER_CHUNKS,
               "cover size": COVER_CHUNKS * CHUNK, "start": hex(START),
               "threads": args.threads,
               "canonical file": {"path": str(CANON_FILE),
                                  "sha256": sha256_of(CANON_FILE)},
               "expanded file": expanded_meta, "scenarios": {}}
    if out.exists():
        old = json.loads(out.read_text())
        if old.get("canonical file", {}).get("sha256") == \
                results["canonical file"]["sha256"] and \
                old.get("expanded file", {}).get("sha256") == \
                expanded_meta["sha256"]:
            results["scenarios"] = old.get("scenarios", {})
        else:
            print("[setup] input files changed - starting fresh results")

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
        if args.session_cap is not None:
            SCENARIOS[name] = SCENARIOS[name][:4] + (args.session_cap,)
        print(f"[scenario {name}] sessions so far: "
              f"{len(results['scenarios'].get(name, {}).get('sessions', []))}")
        run_scenario(name, results, args, targets, workdir, binary)
        v = verdict(name, results, targets)
        if v:
            results["scenarios"][name]["verdict"] = v
            print(f"[scenario {name}] {v['verdict']} "
                  f"({v['sessions']} sessions, "
                  f"{v['killed sessions']} killed) "
                  f"{'; '.join(v['reasons']) if v['reasons'] else ''}")
        save_json(out, results)

    done = [n for n in SCENARIOS if results["scenarios"].get(n, {}).get("verdict")]
    print(f"\ncompleted scenarios: {len(done)}/{len(SCENARIOS)}")
    for name in SCENARIOS:
        v = results["scenarios"].get(name, {}).get("verdict")
        if v:
            print(f"  {name:14s} {v['verdict']}  sessions={v['sessions']} "
                  f"killed={v['killed sessions']}")


if __name__ == "__main__":
    main()
