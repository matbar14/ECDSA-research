#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Эквивалентность трёх источников таблицы на --cover-плане (100k / 4M).

План: COVER_CHUNKS внешних chunk по CHUNK = G*stride от START = 2^134,
контрольные цели в chunk 2 и chunk 3 (j в [1, H] - ядро тайлинга) и
последняя точка покрытия. Один и тот же план прогоняется тремя
источниками таблицы:

  in-memory  --table-scalars canonical_H + --canonical-k3 (RAM-таблица);
  file       --table-binary CBT1 (mmap-таблица, тот же контент);
  resume     file + --checkpoint-file, первая сессия убивается SIGKILL
             на границе chunk (poll checkpoint), вторая --resume до конца.

Сверка бит-в-бит по каждому target:
  - восстановленный k == засеянному скаляру во всех источниках;
  - вектор СТРОГИХ счётчиков (invariant hits, search hits, orbit factors
    scheduled, exact checks, partner recoveries, merged rows, insert
    collisions, giant steps planned/scheduled) идентичен во всех трёх
    источниках; EC-суммы не сравниваются: build-EC в RAM-режиме существует,
    в file-mode равен 0 по конструкции;
  - «giant steps» на found-целях зависят от планировщика потоков (other-
    thread slack при остановке скана), поэтому сравниваются только на
    miss-плане (цель вне покрытия): полный проход детерминирован,
    (COVER_CHUNKS+1)*G шагов бит-в-бит во всех источниках.
Плюс цепочка чекпоинта у resume: found=0 в убитой сессии,
"chunks completed before" == чанкам убитой сессии, found=1 в финале.

Exit 0 iff все источники всех целей бит-в-бит эквивалентны.
"""

import argparse
import json
import re
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
START = 0x4000000000000000000000000000000000
COVER_CHUNKS = 4
G = 300_000                      # гигант-шагов на chunk (скан ~1 c)

failures = []


def pubkey_hex(scalar: int) -> str:
    from cryptography.hazmat.primitives.asymmetric import ec
    numbers = ec.derive_private_key(scalar % N, ec.SECP256K1()) \
        .public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def counter(out, name):
    m = re.search(rf"{re.escape(name)}: (-?\d+)", out)
    return int(m.group(1)) if m else None


def parse_cp(path: Path):
    if not path.exists():
        return None
    fields = {}
    for line in path.read_text().splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            fields[k.strip()] = v.strip()
    return fields


def run_capture(args_list, timeout=240):
    proc = subprocess.run(args_list, text=True, capture_output=True,
                          timeout=timeout)
    return proc


def kill_after_chunks(args_list, cp_path: Path, threshold: int, timeout=240):
    """Сессия с checkpoint; SIGKILL как только chunks >= threshold."""
    proc = subprocess.Popen(args_list, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        cp = parse_cp(cp_path)
        if cp and int(cp.get("chunks completed", -1)) >= threshold:
            proc.kill()
            out, err = proc.communicate()
            return proc.returncode, out, err, True
        if proc.poll() is not None:
            out, err = proc.communicate()
            return proc.returncode, out, err, False
        time.sleep(0.02)
    proc.kill()
    out, err = proc.communicate()
    return proc.returncode, out, err, True


STRICT_COUNTERS = ("invariant hits", "search hits",
                   "orbit factors scheduled", "exact candidate checks",
                   "partner recoveries", "merged rows",
                   "insert collisions", "giant steps planned")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--heads", type=int, required=True)
    ap.add_argument("--canonical", required=True,
                    help="representatives_canonical_H*.txt")
    ap.add_argument("--table-binary", required=True,
                    help="CBT1 (создаётся, если нет, из canonical)")
    ap.add_argument("--out", required=True, help="артефакт JSON")
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()

    H = args.heads
    stride = H
    chunk = G * stride
    binary = Path("./eng_bsgs").resolve()
    canonical = Path(args.canonical).resolve()
    table = Path(args.table_binary).resolve()
    tmp = table.parent

    # Dump the CBT1 once if absent (in-memory build from the same file).
    if not table.exists():
        print(f"[dump] building {table.name} from {canonical.name}...")
        seed_target = START + 7 * stride + 12345
        dump_cmd = [str(binary), "--target", pubkey_hex(seed_target),
                    "--start", hex(START), "--size", hex(10 * chunk),
                    "--segment", hex(10 * chunk), "--threads", "2",
                    "--progress-interval", "0",
                    "--table-scalars", str(canonical),
                    "--glv-orbit", "--canonical-k3", "--partner-table",
                    "--dump-table-binary", str(table)]
        proc = run_capture(dump_cmd)
        assert proc.returncode == 0 and table.exists(), proc.stderr[-400:]
    else:
        print(f"[dump] using existing {table.name}")

    # Targets: chunk 2, chunk 3, last point of coverage (all found) and
    # a miss plan outside the cover (full deterministic pass). j in [1, H]
    # (core tiling, stride == H -> contiguous coverage).
    j2 = H // 2 + 12345
    j3 = H // 3 + 77777
    targets = {
        "chunk2": {"offset": 2 * chunk + j2, "kill after": 1,
                   "found": True},
        "chunk3": {"offset": 3 * chunk + j3, "kill after": 2,
                   "found": True},
        "last": {"offset": COVER_CHUNKS * chunk - 1, "kill after": 2,
                 "found": True},
        "miss": {"offset": COVER_CHUNKS * chunk + 12345, "kill after": 1,
                 "found": False},
    }
    for spec in targets.values():
        spec["scalar"] = (START + spec["offset"]) % N
        spec["pubkey"] = pubkey_hex(spec["scalar"])

    results = {"heads": H, "stride": stride, "chunk": chunk,
               "cover chunks": COVER_CHUNKS, "G": G,
               "canonical": canonical.name, "table": table.name,
               "targets": {}, "sources": {}}
    out_path = Path(args.out)

    base_geom = ["--start", hex(START), "--size", str(chunk),
                 "--segment", str(chunk), "--threads", str(args.threads),
                 "--progress-interval", "0",
                 "--cover-start", hex(START),
                 "--cover-size", hex(COVER_CHUNKS * chunk)]
    source_flags = {
        "in-memory": ["--table-scalars", str(canonical),
                      "--glv-orbit", "--canonical-k3", "--partner-table"],
        "file": ["--table-binary", str(table),
                 "--glv-orbit", "--partner-table"],
    }

    for tname, spec in targets.items():
        print(f"[target {tname}] offset {spec['offset']}")
        runs = {}
        for source in ("in-memory", "file"):
            cmd = ([str(binary), "--target", spec["pubkey"]] + base_geom
                   + source_flags[source])
            proc = run_capture(cmd)
            runs[source] = proc
        # resume source: killed first session + resume to completion
        cp = tmp / f"cp_eq_{tname}.txt"
        if cp.exists():
            cp.unlink()
        cmd0 = ([str(binary), "--target", spec["pubkey"]] + base_geom
                + source_flags["file"]
                + ["--checkpoint-file", str(cp), "--checkpoint-steps", "1"])
        rc0, out0, err0, killed = kill_after_chunks(
            cmd0, cp, spec["kill after"])
        cp_kill = parse_cp(cp)
        cmd1 = [str(binary), "--target", spec["pubkey"]] + base_geom \
            + source_flags["file"] \
            + ["--checkpoint-file", str(cp), "--checkpoint-steps", "1",
               "--resume"]
        proc_r = run_capture(cmd1)
        runs["resume"] = proc_r
        runs["resume-meta"] = {
            "killed session rc": rc0, "killed": killed,
            "cp kill": cp_kill,
        }

        # ---- bit-for-bit comparison ----
        rec = {}
        cnt = {}
        for source in ("in-memory", "file", "resume"):
            out = runs[source].stdout
            m = re.search(r"recovered private key: 0x([0-9a-f]+)", out)
            rec[source] = int(m.group(1), 16) if m else None
            cnt[source] = {name: counter(out, name) for name in
                           STRICT_COUNTERS + ("giant steps",)}
        rec_ok = (all(v == spec["scalar"] for v in rec.values())
                  if spec["found"] else
                  all(v is None for v in rec.values()))
        strict_equal = all(
            cnt[s][k] == cnt["in-memory"][k]
            for s in ("file", "resume") for k in STRICT_COUNTERS)
        zeros_ok = all(cnt[s]["merged rows"] == 0
                       and cnt[s]["insert collisions"] == 0
                       for s in ("in-memory", "file", "resume"))
        # giant steps: bit-for-bit only on the deterministic miss plan
        plan_steps = (COVER_CHUNKS + 1) * G
        if spec["found"]:
            steps_ok = all(0 < cnt[s]["giant steps"] <= plan_steps
                           for s in ("in-memory", "file", "resume"))
        else:
            steps_ok = all(cnt[s]["giant steps"] == plan_steps
                           for s in ("in-memory", "file", "resume"))
        hits_zero_ok = (True if spec["found"] else
                        all(cnt[s]["search hits"] == 0
                            and cnt[s]["invariant hits"] == 0
                            for s in ("in-memory", "file", "resume")))
        rcs_ok = all(runs[s].returncode == (0 if spec["found"] else 1)
                     for s in ("in-memory", "file", "resume"))
        # checkpoint chain
        cp_final = parse_cp(cp)
        chain_ok = (cp_kill is not None and cp_final is not None
                    and cp_kill.get("found") == "0"
                    and int(cp_kill.get("chunks completed", -1))
                    >= spec["kill after"]
                    and counter(runs["resume"].stdout,
                                "chunks completed before this session")
                    == int(cp_kill.get("chunks completed", -1))
                    and (cp_final.get("found") == ("1" if spec["found"]
                                                   else "0")))
        ok = (rec_ok and strict_equal and zeros_ok and steps_ok
              and hits_zero_ok and rcs_ok and chain_ok)
        status = "PASS" if ok else "FAIL"
        print(f"  k bit-for-bit all sources: {rec_ok}; strict counters "
              f"equal: {strict_equal}; steps ok: {steps_ok}; "
              f"cp chain: {chain_ok} -> {status}")
        if not ok:
            print("   counters:", json.dumps(cnt, indent=2))
            for s in ("in-memory", "file", "resume"):
                print(f"   --- {s} tail:", runs[s].stdout[-300:])
        if not ok:
            failures.append(tname)
        results["targets"][tname] = {
            "offset": spec["offset"], "scalar hex": hex(spec["scalar"]),
            "expect found": spec["found"],
            "recovered": {s: (hex(v) if v else None)
                          for s, v in rec.items()},
            "counters": cnt, "k bit-for-bit": rec_ok,
            "strict counters equal": strict_equal,
            "steps ok": steps_ok, "rcs ok": rcs_ok,
            "cp kill": cp_kill, "cp chain ok": chain_ok,
            "pass": ok,
        }
        tmp_json = out_path.with_suffix(".tmp")
        tmp_json.write_text(json.dumps(results, indent=1))
        tmp_json.rename(out_path)

    print(f"\nEquivalence H={H}: "
          f"{'ALL PASS' if not failures else 'FAILED: ' + str(failures)}")
    print(f"[artifact] {out_path}")
    sys.exit(0 if not failures else 1)


if __name__ == "__main__":
    main()
