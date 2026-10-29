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
    "merged rows",
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


def partner_targets(rng, stride, giant_count, count, heads):
    """Детерминированный список орбитальных целей для данного seed-потока."""
    targets = []
    for _ in range(count):
        j = rng.randrange(heads)
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


def save_json(out, results):
    """Атомарная запись частичных результатов (resume между вызовами)."""
    tmp = out.with_suffix(".tmp")
    tmp.write_text(json.dumps(results, indent=2) + "\n")
    os.replace(tmp, out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    parser.add_argument("--heads", type=int, default=100_000,
                        help="число голов орбит H")
    parser.add_argument("--trials", type=int, default=12,
                        help="орбитальных целей на таблицу")
    parser.add_argument("--threads", type=int,
                        default=min(8, os.cpu_count() or 1))
    parser.add_argument("--phases", default="tp,ec",
                        help="запускаемые фазы: tp (throughput), ec (economy)")
    parser.add_argument("--output", default="../table_representative_sweep.json")
    global args
    args = parser.parse_args()
    run_tp = "tp" in args.phases.split(",")
    run_ec = "ec" in args.phases.split(",")

    heads = list(range(1, args.heads + 1))
    stride_c = args.heads                    # canonical: одна строка на орбиту
    stride_e = 3 * args.heads                # expanded: три строки на орбиту
    size = stride_e * stride_e               # одно и то же окно для всех

    # ---- resume: дозапуск в существующий JSON, цели детерминированы ----
    out = Path(args.output)
    results = None
    if out.exists():
        try:
            loaded = json.loads(out.read_text())
            if loaded.get("heads") == args.heads:
                results = loaded
        except (json.JSONDecodeError, OSError):
            results = None
    fresh = results is None
    if fresh:
        results = {"heads": args.heads, "window size": size,
                   "start": hex(START), "trials": args.trials,
                   "throughput": {}, "economy": {}, "failures": []}
    results.setdefault("throughput", {})
    results.setdefault("economy", {})
    results.setdefault("failures", [])
    results["trials"] = args.trials
    results["window size"] = size
    results["start"] = hex(START)

    # Детерминированные цели: отдельные seed-потоки, независимые от того,
    # какие фазы/планы уже выполнены (устойчивый resume).
    rng_out = random.Random(20261101)
    outside = START + size + rng_out.randrange(1, 1 << 32)
    giant_count_c = (size - 1) // stride_c + 1
    giant_count_e = (size - 1) // stride_e + 1
    targets_c = partner_targets(random.Random(20261102), stride_c,
                                giant_count_c, args.trials, args.heads)
    targets_e = partner_targets(random.Random(20261103), stride_e,
                                giant_count_e, args.trials, args.heads)

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
        if fresh:
            print(f"[setup] heads={args.heads}, window size={size}, "
                  f"expanded rows={len(expanded)}, canonical rows={len(canonical)}")
        else:
            done_tp = len(results["throughput"])
            done_ec = sum(1 for v in results["economy"].values()
                          if v.get("all correct"))
            print(f"[resume] heads={args.heads}: throughput планов готово "
                  f"{done_tp}/4, economy планов готово {done_ec}/4")

        # ---------- throughput: цель вне покрытия, полный проход ----------
        throughput_plans = [
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
        if run_tp:
            for name, extra in throughput_plans:
                if results["throughput"].get(name, {}).get(
                        "table rows stored") is not None:
                    print(f"[throughput] {name}: уже выполнен, пропуск")
                    continue
                run = run_engine(args.binary, pubkey_hex(outside), size, extra)
                rows_stored = (run["counters"].get("table rows stored")
                               or run["counters"].get("raw baby rows built"))
                tp_entry = {
                    "table rows stored": rows_stored,
                    "unique baby keys": run["counters"].get("unique baby rows built"),
                    "merged rows": run["counters"].get("merged rows", 0),
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
                results["throughput"][name] = tp_entry
                save_json(out, results)
                steps_s = tp_entry["giant steps/second"]
                steps_s_text = f"{steps_s:.0f}" if steps_s is not None else "n/a"
                print(f"[throughput] {name}: rows={rows_stored}, "
                      f"steps={tp_entry['giant steps']}, "
                      f"steps/s={steps_s_text}, "
                      f"EC-ops={tp_entry['total EC operations']}, "
                      f"rss={run['peak rss kb']} kB")
                if run["recovered"] is not None:
                    results["failures"].append(
                        f"throughput {name}: unexpected recovery")
                    save_json(out, results)

        # ---------- economy: орбитальные цели ----------
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
        if run_ec:
            for name, targets, stride, base_extra in economy_plans:
                prior = results["economy"].get(name)
                if prior and prior.get("all correct") \
                        and len(prior.get("runs", [])) == args.trials:
                    print(f"[economy] {name}: уже выполнен, пропуск")
                    continue
                runs = []
                for t in targets:
                    run = run_engine(args.binary, t["pubkey"], size, base_extra)
                    ok = run["recovered"] == t["scalar"]
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
                                    "merged rows")},
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
                        results["failures"].append(
                            f"{name}: head={t['head']} r={t['rotation']} "
                            f"sign={t['sign']} stride={stride}")
                exact = [r["exact candidate checks"] for r in runs
                         if r.get("exact candidate checks") is not None]
                ec_ops = [r["total EC operations"] for r in runs
                          if r.get("total EC operations") is not None]
                results["economy"][name] = {
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
                save_json(out, results)
                print(f"[economy] {name}: all_correct="
                      f"{results['economy'][name]['all correct']}, "
                      f"avg exact/hit="
                      f"{results['economy'][name]['avg exact checks per hit']}")

    failures = results["failures"]
    save_json(out, results)

    print("\n=== summary ===")
    for name, tp in results["throughput"].items():
        steps_s = tp["giant steps/second"]
        steps_s_text = f"{steps_s:.0f}" if steps_s is not None else "n/a"
        print(f"throughput {name}: rows={tp['table rows stored']}, "
              f"steps/s={steps_s_text}, "
              f"EC-ops={tp['total EC operations']}, "
              f"rss={tp['peak rss kb']} kB")
    for name, ec_data in results["economy"].items():
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
