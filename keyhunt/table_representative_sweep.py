#!/usr/bin/env python3
"""Движковый замер таблицы представителей полного диапазона.

PROGRESS_MEMO.md («Каноническая K3-таблица», следующий шаг): разрежение
baby-таблицы по орбитам на движке -- `--table-scalars` принимает произвольные
скаляры, поэтому таблицу можно собрать из канонических представителей орбит,
разбросанных по всему [1, n), и измерить пропускную способность giant-фазы
при уменьшенном числе строк.

Три варианта таблицы на одних головах орбит j in [0, H):

  continuous  непрерывное окно [0, 3H) из --baby-rows (производственный
              базовый режим, stride = 3H);
  expanded    явный список (j, lambda*j, lambda^2*j) для каждой головы -
              3H строк, "наивное" хранение трёх x-представителей орбиты,
              K3-bucket содержит 3 строки (stride = 3H);
  canonical   канонический представитель орбиты
              c = min_r min(lambda^r*j mod n, n - lambda^r*j mod n) -
              H строк, одна строка на орбиту (stride = H).

Замеры.

  throughput  цель заведомо вне покрытия -> полный giant-проход по окну
              S = (3H)^2; сравниваются строки таблицы, время сборки,
              giant steps/second и пиковый RSS при одном и том же окне.

  economy     12 орбитальных целей на таблицу: k = start + i*stride + s*,
              s* = +/-lambda^a*j (a in {0,1,2}); шесть-ветвевой базовый
              режим против (d01,d12)-маршрутизации; на expanded-таблице
              дополнительно --canonical-k3 (входные 3H строк сливаются в H).

Каждое восстановление сверяется с точным скаляром цели, поэтому PASS/FAIL
в выводе означает фактическую корректность композиции
"каноническое хранение + partner-expansion", а не только счётчики.
"""

import argparse
import json
import os
import random
import re
import subprocess
import tempfile
import threading
import time
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
LAMBDA2 = 0xAC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE
START = 1 << 134

COUNTERS = (
    "giant steps",
    "invariant hits",
    "orbit factors scheduled",
    "exact candidate checks",
    "EC scalar multiplications",
    "EC point additions",
    "total EC operations",
    "signature routes",
    "signature anomalies (6-branch fallback)",
    "signature pair collapses",
    "signature pair mismatches",
    "orbit duplicate rows skipped",
    "canonical rows merged",
    "partner recoveries",
    "raw baby rows built",
    "unique baby rows built",
)


def pubkey_hex(scalar):
    numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def canonical_representative(j):
    """c = min_r min(lambda^r*j mod n, n - lambda^r*j mod n) - зеркало
    canonical_orbit_scalar из eng_bsgs.cpp."""
    best = None
    v = j % N
    for _ in range(3):
        candidate = min(v, N - v) if v else 0
        if best is None or candidate < best:
            best = candidate
        v = v * LAMBDA % N
    return best


def orbit_expansion(j):
    """Три x-представителя орбиты: скаляры j, lambda*j, lambda^2*j."""
    return [j % N, j * LAMBDA % N, j * LAMBDA2 % N]


def write_scalars(path, scalars, comment):
    lines = [f"# {comment}"] + [hex(s) for s in scalars]
    Path(path).write_text("\n".join(lines) + "\n")


def peak_rss_kb(pid, sink, stop_event):
    """Опрос VmHWM из /proc - пик resident-памяти процесса движка."""
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


def run_engine(binary, target, size, extra):
    cmd = [binary, "--target", target, "--start", hex(START),
           "--size", str(size), "--segment", str(size),
           "--threads", str(args.threads)] + extra
    sink = [0]
    stop_event = threading.Event()
    t0 = time.perf_counter()
    proc = subprocess.Popen(cmd, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    monitor = threading.Thread(target=peak_rss_kb,
                               args=(proc.pid, sink, stop_event),
                               daemon=True)
    monitor.start()
    stdout, stderr = proc.communicate()
    stop_event.set()
    monitor.join(timeout=0.1)
    wall = time.perf_counter() - t0

    match = re.search(r"recovered private key: 0x([0-9a-fA-F]+)", stdout)
    recovered = int(match.group(1), 16) if match else None
    counters = {}
    for name in COUNTERS:
        m = re.search(rf"{re.escape(name)}: (\d+)", stdout)
        if m:
            counters[name] = int(m.group(1))
    timing = {}
    for name in ("table build seconds", "giant search seconds",
                 "giant steps/second", "total seconds"):
        m = re.search(rf"{re.escape(name)}: ([0-9.eE+-]+)", stdout)
        if m:
            timing[name] = float(m.group(1))
    stored = re.search(r"table rows stored: (\d+)", stdout)
    if stored:
        counters["table rows stored"] = int(stored.group(1))
    if proc.returncode not in (0, 1):
        print(stdout)
        print(stderr)
    return {
        "returncode": proc.returncode,
        "recovered": recovered,
        "counters": counters,
        "timing": timing,
        "wall seconds": round(wall, 4),
        "peak rss kb": sink[0],
        "stdout": stdout,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    parser.add_argument("--heads", type=int, default=100_000,
                        help="число голов орбит H")
    parser.add_argument("--trials", type=int, default=12,
                        help="орбитальных целей на таблицу")
    parser.add_argument("--threads", type=int,
                        default=min(8, os.cpu_count() or 1))
    parser.add_argument("--output", default="../table_representative_sweep.json")
    global args
    args = parser.parse_args()

    rng = random.Random(20261101)
    # j = 0 даёт нулевые скаляры орбиты, движок их отвергает - головы с 1.
    heads = list(range(1, args.heads + 1))
    stride_c = args.heads                    # canonical: одна строка на орбиту
    stride_e = 3 * args.heads                # expanded: три строки на орбиту
    size = stride_e * stride_e               # одно и то же окно для всех

    failures = []
    results = {"heads": args.heads, "window size": size,
               "start": hex(START), "trials": args.trials,
               "throughput": {}, "economy": {}}

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        expanded_path = tmp / "expanded_scalars.txt"
        canonical_path = tmp / "canonical_scalars.txt"
        expanded = []
        canonical = []
        for j in heads:
            expanded.extend(orbit_expansion(j))
            canonical.append(canonical_representative(j))
        write_scalars(expanded_path, expanded,
                      f"3 x-representatives per orbit, {len(expanded)} rows")
        write_scalars(canonical_path, canonical,
                      f"one canonical representative per orbit, {len(canonical)} rows")
        results["table files"] = {
            "expanded rows": len(expanded),
            "canonical rows": len(canonical),
        }
        print(f"[setup] heads={args.heads}, window size={size}, "
              f"expanded rows={len(expanded)}, canonical rows={len(canonical)}")

        # ---------- throughput: цель вне покрытия, полный проход ----------
        outside = START + size + rng.randrange(1, 1 << 32)
        target_out = pubkey_hex(outside)
        throughput = {}
        plans = [
            ("continuous 3H rows", ["--baby-rows", str(stride_e),
                                    "--glv-orbit", "--phase-signature"]),
            ("expanded 3H rows", ["--table-scalars", str(expanded_path),
                                  "--glv-orbit", "--phase-signature"]),
            ("canonical H rows", ["--table-scalars", str(canonical_path),
                                  "--glv-orbit", "--phase-signature"]),
            ("expanded->canonical merge",
             ["--table-scalars", str(expanded_path), "--glv-orbit",
              "--phase-signature", "--canonical-k3"]),
        ]
        for name, extra in plans:
            run = run_engine(args.binary, target_out, size, extra)
            rows_stored = (run["counters"].get("table rows stored")
                           or run["counters"].get("raw baby rows built"))
            throughput[name] = {
                "table rows stored": rows_stored,
                "unique baby keys": run["counters"].get("unique baby rows built"),
                "canonical rows merged": run["counters"].get("canonical rows merged", 0),
                "giant steps": run["counters"].get("giant steps"),
                "EC scalar multiplications": run["counters"].get("EC scalar multiplications"),
                "EC point additions": run["counters"].get("EC point additions"),
                "total EC operations": run["counters"].get("total EC operations"),
                "build seconds": run["timing"].get("table build seconds"),
                "search seconds": run["timing"].get("giant search seconds"),
                "giant steps/second": run["timing"].get("giant steps/second"),
                "wall seconds": run["wall seconds"],
                "peak rss kb": run["peak rss kb"],
                "recovered unexpectedly": run["recovered"] is not None,
            }
            steps_s = throughput[name]["giant steps/second"]
            steps_s_text = f"{steps_s:.0f}" if steps_s is not None else "n/a"
            print(f"[throughput] {name}: rows="
                  f"{throughput[name]['table rows stored']}, "
                  f"steps={throughput[name]['giant steps']}, "
                  f"steps/s={steps_s_text}, "
                  f"EC-ops={throughput[name]['total EC operations']}, "
                  f"rss={run['peak rss kb']} kB")
            if run["recovered"] is not None:
                failures.append(f"throughput {name}: unexpected recovery")
        results["throughput"] = throughput

        # ---------- economy: орбитальные цели ----------
        # Цели строятся под stride конкретной таблицы:
        # k = start + i*stride + s*, s* = +/-(lambda^a * j).
        def partner_targets(stride, giant_count, count):
            targets = []
            for _ in range(count):
                j = rng.randrange(args.heads)
                a = rng.randrange(3)
                sign = rng.choice((1, -1))
                i = rng.randrange(giant_count // 2, giant_count)
                s = sign * pow(LAMBDA, a, N) * j % N
                # скаляр цели обязан лежать в [1, n) для derive_private_key;
                # сравнение восстановленного значения идёт по модулю n
                scalar = (START + i * stride + s) % N
                targets.append({
                    "scalar": scalar, "head": j, "rotation": a, "sign": sign,
                    "stride index": i, "pubkey": pubkey_hex(scalar),
                })
            return targets

        giant_count_c = (size - 1) // stride_c + 1
        giant_count_e = (size - 1) // stride_e + 1
        targets_c = partner_targets(stride_c, giant_count_c, args.trials)
        targets_e = partner_targets(stride_e, giant_count_e, args.trials)

        economy_plans = [
            ("canonical six-branch", targets_c, stride_c,
             ["--table-scalars", str(canonical_path), "--glv-orbit",
              "--partner-table"]),
            ("canonical routed", targets_c, stride_c,
             ["--table-scalars", str(canonical_path), "--glv-orbit",
              "--phase-signature", "--partner-table"]),
            ("expanded routed", targets_e, stride_e,
             ["--table-scalars", str(expanded_path), "--glv-orbit",
              "--phase-signature", "--partner-table"]),
            ("expanded routed + canonical merge", targets_e, stride_e,
             ["--table-scalars", str(expanded_path), "--glv-orbit",
              "--phase-signature", "--partner-table", "--canonical-k3"]),
        ]
        economy = {}
        for name, targets, stride, base_extra in economy_plans:
            runs = []
            for t in targets:
                run = run_engine(args.binary, t["pubkey"], size, base_extra)
                ok = run["recovered"] == t["scalar"]
                hit = run["recovered"] is not None
                entry = {
                    "correct": ok,
                    "head": t["head"], "rotation": t["rotation"],
                    "sign": t["sign"], "stride index": t["stride index"],
                    "stride": stride,
                    **{k: v for k, v in run["counters"].items()
                       if k in ("exact candidate checks",
                                "orbit factors scheduled",
                                "EC scalar multiplications",
                                "EC point additions",
                                "total EC operations",
                                "signature routes",
                                "signature pair collapses",
                                "signature pair mismatches",
                                "signature anomalies (6-branch fallback)",
                                "partner recoveries",
                                "table rows stored",
                                "canonical rows merged")},
                    "search seconds": run["timing"].get("giant search seconds"),
                }
                runs.append(entry)
                status = "PASS" if ok else "FAIL"
                print(f"[economy] {name}: head={t['head']} r={t['rotation']} "
                      f"sign={t['sign']:+d} -> {status}, "
                      f"exact={entry.get('exact candidate checks')}, "
                      f"EC-ops={entry.get('total EC operations')}, "
                      f"routes={entry.get('signature routes')}")
                if not ok:
                    failures.append(
                        f"{name}: head={t['head']} r={t['rotation']} "
                        f"sign={t['sign']} stride={stride}")
            exact = [r["exact candidate checks"] for r in runs
                     if r.get("exact candidate checks") is not None]
            ec_ops = [r["total EC operations"] for r in runs
                      if r.get("total EC operations") is not None]
            economy[name] = {
                "stride": stride,
                "table rows stored": runs[0].get("table rows stored"),
                "all correct": all(r["correct"] for r in runs),
                "avg exact checks per hit": round(sum(exact) / len(exact), 3)
                if exact else None,
                "avg total EC ops per run": round(sum(ec_ops) / len(ec_ops), 1)
                if ec_ops else None,
                "partner recoveries": sum(r.get("partner recoveries", 0)
                                          for r in runs),
                "signature pair collapses": sum(
                    r.get("signature pair collapses", 0) for r in runs),
                "signature pair mismatches": sum(
                    r.get("signature pair mismatches", 0) for r in runs),
                "signature anomalies": sum(
                    r.get("signature anomalies (6-branch fallback)", 0)
                    for r in runs),
                "runs": runs,
            }
            print(f"[economy] {name}: all_correct="
                  f"{economy[name]['all correct']}, "
                  f"avg exact/hit={economy[name]['avg exact checks per hit']}")
        results["economy"] = economy

    results["failures"] = failures
    out = Path(args.output)
    out.write_text(json.dumps(results, indent=2) + "\n")

    print("\n=== summary ===")
    for name, tp in throughput.items():
        steps_s = tp["giant steps/second"]
        steps_s_text = f"{steps_s:.0f}" if steps_s is not None else "n/a"
        print(f"throughput {name}: rows={tp['table rows stored']}, "
              f"steps/s={steps_s_text}, "
              f"EC-ops={tp['total EC operations']}, "
              f"rss={tp['peak rss kb']} kB")
    for name, ec_data in economy.items():
        print(f"economy {name}: all_correct={ec_data['all correct']}, "
              f"avg exact/hit={ec_data['avg exact checks per hit']}")
    if failures:
        print(f"\nFAILURES: {len(failures)}")
        for f in failures:
            print(" -", f)
        raise SystemExit(1)
    print("\nALL PASS")


if __name__ == "__main__":
    main()
