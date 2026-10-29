#!/usr/bin/env python3
"""Рандомизированная проверка 400k-аномалии сборки: direct canonical vs merge.

PROGRESS_MEMO.md, «H-кривая», вывод 2: при H = 400k merge-сборка (4.41s)
обогнала прямой canonical-вход (5.36s), хотя merge делает строго больше
работы до EC (читает 1.2M строк, канонизирует и сортирует их; EC-работа
строк одинакова). Гипотеза memo - подогрев файлового кэша предыдущим
планом; в одном порядке планов гипотезу отделить от скрытой стоимости
direct canonical невозможно.

Этот скрипт снимает обе сборки back-to-back с перемешанным порядком
(seed 20261105), 3 повтора, против цели вне покрытия. Условия совпадают
с h_curve_400k.json: START = 2^134, окно S = (3H)^2, 2 потока.

Обе таблицы строят одинаковые 400k EC-строк, поэтому сравниваются
"table build seconds" (и, для контроля, giant steps/second).

Вердикт:
  mean(canonical) ~ mean(merge)  -> аномалия была кэш/порядок планов,
                                    точка 400k недостоверна (cold-cache)
  mean(canonical) > mean(merge) на ~20%  -> у direct canonical скрытая
                                    стоимость на 400k: фиксовать режим
                                    до длинного прогона.

Дополнительно проверяется статический предусловие "неинициализированный
_O": curve.Init() (eng_bsgs.cpp main, строка ~1019) вызывается до парсинга
--table-scalars и canonical pre-pass, поэтому ни canonical, ни merge не
платят за порядок инициализации (см. вывод скрипта).

Каждый прогон записывается атомарно в JSON сразу после запуска - скрипт
resume-безопасен в убивающей среде.
"""

import argparse
import json
import os
import random
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "keyhunt"))
from table_representative_sweep import (  # noqa: E402
    START, canonical_representative, orbit_expansion, pubkey_hex,
)

OUT_DEFAULT = Path(__file__).resolve().parent / "h_curve_400k_anomaly.json"
SEED = 20261105


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


def run_engine(binary, target, size, extra, threads):
    cmd = [binary, "--target", target, "--start", hex(START),
           "--size", str(size), "--segment", str(size),
           "--threads", str(threads)] + extra
    sink = [0]
    stop_event = threading.Event()
    t0 = time.perf_counter()
    proc = subprocess.Popen(cmd, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    monitor = threading.Thread(target=peak_rss_kb,
                               args=(proc.pid, sink, stop_event), daemon=True)
    monitor.start()
    stdout, stderr = proc.communicate()
    stop_event.set()
    monitor.join(timeout=0.1)
    wall = time.perf_counter() - t0
    timing = {}
    for name in ("table build seconds", "giant search seconds",
                 "giant steps/second", "total seconds"):
        m = re.search(rf"{re.escape(name)}: ([0-9.eE+-]+)", stdout)
        if m:
            timing[name] = float(m.group(1))
    counters = {}
    for name in ("giant steps", "table rows stored", "canonical rows merged",
                 "raw baby rows built", "unique baby rows built",
                 "EC scalar multiplications", "EC point additions",
                 "total EC operations"):
        m = re.search(rf"{re.escape(name)}: (\d+)", stdout)
        if m:
            counters[name] = int(m.group(1))
    if proc.returncode not in (0, 1):
        print(stdout, file=sys.stderr)
        print(stderr, file=sys.stderr)
    return {
        "returncode": proc.returncode,
        "timing": timing, "counters": counters,
        "wall seconds": round(wall, 4), "peak rss kb": sink[0],
        "stdout": stdout,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="keyhunt/eng_bsgs")
    parser.add_argument("--heads", type=int, default=400_000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--pairs", type=int, default=6,
                        help="пар C/M при --design paired")
    parser.add_argument("--design", choices=("shuffle", "paired"),
                        default="paired",
                        help="shuffle: 2*repeats перемешанных прогонов; "
                             "paired: чередующиеся пары M,C / C,M "
                             "(парная разность гасит дрейф среды)")
    parser.add_argument("--window-divisor", type=int, default=8,
                        help="окно делится на K, чтобы пара прогонов "
                             "укладывалась в короткий интервал дрейфа; "
                             "build от окна не зависит")
    parser.add_argument("--threads", type=int, default=2,
                        help="h_curve снимался при 2 потоках")
    parser.add_argument("--output", default=str(OUT_DEFAULT))
    args = parser.parse_args()

    heads = list(range(1, args.heads + 1))
    stride_c = args.heads
    stride_e = 3 * args.heads
    size = (stride_e * stride_e) // args.window_divisor
    # точный (3H)^2 делится на 8 нацело, но подстраховка на целые H
    size -= size % stride_e

    out = Path(args.output)
    base = None
    if out.exists():
        try:
            loaded = json.loads(out.read_text())
            if loaded.get("heads") == args.heads:
                base = loaded
        except (json.JSONDecodeError, OSError):
            base = None
    fresh = base is None
    if fresh:
        if args.design == "paired":
            plan = []
            for i in range(args.pairs):
                first = "merge" if i % 2 == 1 else "canonical"
                plan += [first, "merge" if first == "canonical" else "canonical"]
            note = ("строго чередующиеся пары, ориентация C,M / M,C "
                    "переключается по чётности i; парная разность "
                    "build(canonical) - build(merge) гасит дрейф среды")
        else:
            rng_seq = random.Random(SEED)
            plan = ["canonical", "merge"] * args.repeats
            rng_seq.shuffle(plan)
            note = (f"2*{args.repeats} перемешанных прогонов "
                    f"(seed {SEED}); при дрейфе среды смещает "
                    f"режимы в разные полосы скорости")
        base = {
            "heads": args.heads, "seed": SEED, "design": args.design,
            "order": plan, "window divisor": args.window_divisor,
            "window size": size, "start": hex(START),
            "threads": args.threads, "note": note,
            "runs": [], "verdict": None,
        }
        base["outside target"] = pubkey_hex(
            START + size + random.Random(SEED + 1).randrange(1, 1 << 32))
    rng_plan = base["order"]

    with __import__("tempfile").TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        expanded_path = tmp / "expanded_scale.txt"
        canonical_path = tmp / "canonical_scale.txt"
        expanded, canonical = [], []
        for j in heads:
            expanded.extend(orbit_expansion(j))
            canonical.append(canonical_representative(j))
        for path, rows, comment in (
                (expanded_path, expanded, "3 x-представителя на орбиту"),
                (canonical_path, canonical, "canonical представитель на орбиту")):
            lines = [f"# {comment}, {len(rows)} строк"] + [hex(s) for s in rows]
            path.write_text("\n".join(lines) + "\n")

        plans = {
            "canonical": ["--table-scalars", str(canonical_path),
                          "--glv-orbit", "--phase-signature"],
            "merge": ["--table-scalars", str(expanded_path), "--glv-orbit",
                      "--phase-signature", "--canonical-k3"],
        }
        if fresh:
            print(f"[setup] heads={args.heads}, window={size}, "
                  f"threads={args.threads}")
            print(f"[setup] порядок прогонов: {' -> '.join(rng_plan)}")

        resume_from = len(base["runs"])
        for run_index in range(resume_from, len(rng_plan)):
            name = rng_plan[run_index]
            t0 = time.perf_counter()
            run = run_engine(args.binary, base["outside target"], size,
                             plans[name], args.threads)
            entry = {
                "run index": run_index, "mode": name,
                **run["timing"], "giant steps": run["counters"].get("giant steps"),
                "table rows stored": run["counters"].get("table rows stored"),
                "canonical rows merged": run["counters"].get("canonical rows merged"),
                "raw baby rows built": run["counters"].get("raw baby rows built"),
                "unique baby rows built": run["counters"].get("unique baby rows built"),
                "total EC operations": run["counters"].get("total EC operations"),
                "wall seconds": run["wall seconds"],
                "peak rss kb": run["peak rss kb"],
                "probe overhead seconds": round(time.perf_counter() - t0
                                                - run["wall seconds"], 3),
            }
            base["runs"].append(entry)
            _atomic_write(out, base)
            print(f"[run {run_index}] {name}: build={entry['table build seconds']}s "
                  f"search={entry['giant search seconds']}s "
                  f"steps/s={entry['giant steps/second']:,.0f} "
                  f"rss={entry['peak rss kb'] // 1024} MB "
                  f"merged={entry['canonical rows merged']}")

    # ---- вердикт по завершённым прогонам ----
    runs = base["runs"]
    by_mode = {}
    for r in runs:
        by_mode.setdefault(r["mode"], []).append(r["table build seconds"])
    canon = by_mode.get("canonical", [])
    merge = by_mode.get("merge", [])
    verdict = {
        "completed": [{"mode": mo, "n": len(v),
                       "mean build s": round(sum(v) / len(v), 3),
                       "min build s": round(min(v), 3),
                       "max build s": round(max(v), 3),
                       "median build s": round(sorted(v)[len(v) // 2], 3),
                       "runs": [round(x, 3) for x in v]}
                      for mo, v in by_mode.items()],
        "order": base.get("order", []),
        "per-run": [{"mode": r["mode"],
                     "build": r["table build seconds"],
                     "steps/s": r["giant steps/second"]} for r in runs],
    }
    if base.get("design") == "paired":
        # парная разность: canonical - merge внутри каждой пары (i=2k, 2k+1)
        pairs = []
        for k in range(len(runs) // 2):
            a, b = runs[2 * k], runs[2 * k + 1]
            ca = a if a["mode"] == "canonical" else b
            me = b if a["mode"] == "canonical" else a
            if ca["mode"] != "canonical" or me["mode"] != "merge":
                continue
            pairs.append({
                "pair": k,
                "canonical build": ca["table build seconds"],
                "merge build": me["table build seconds"],
                "diff canonical-merge": round(
                    ca["table build seconds"] - me["table build seconds"], 3),
            })
        verdict["pairs"] = pairs
        if pairs:
            diffs = sorted(p["diff canonical-merge"] for p in pairs)
            med_diff = diffs[len(diffs) // 2]
            med_merge = sorted(p["merge build"] for p in pairs)[
                len(pairs) // 2]
            verdict["median diff s"] = med_diff
            verdict["median diff / median merge"] = round(
                med_diff / med_merge, 3)
            if med_diff > 0.10 * med_merge:
                verdict["state"] = ("NOT cache: прямая canonical стабильно "
                                    "медленнее merge >10% - скрытая "
                                    "стоимость реальна; фикс режима до "
                                    "длинного прогона")
            elif med_diff < -0.05 * med_merge:
                verdict["state"] = (f"cache/order: canonical БЫСТРЕЕ merge "
                                    f"на {-med_diff:.2f}s - соответствует "
                                    f"структуре работ (merge несёт лишний "
                                    f"pre-EC проход); записанная точка "
                                    f"5.36/4.41 недостоверна")
            else:
                verdict["state"] = ("равно в пределах шума (<=10%): "
                                    "аномалия 5.36/4.41 не является "
                                    "свойством режима - кэш/порядок "
                                    "планов; canonical остаётся "
                                    "оптимизированным режимом")
        else:
            verdict["state"] = "incomplete"
    else:
        if canon and merge:
            c, m = sorted(canon), sorted(merge)
            verdict["canonical over merge (median)"] = round(
                c[len(c) // 2] / m[len(m) // 2], 3)
            verdict["state"] = ("shuffle: см. paired-верdict "
                                "(h_curve_400k_anomaly_paired.json) - "
                                "дрейф среды подтверждён")
        else:
            verdict["state"] = "incomplete"
    base["verdict"] = verdict
    _atomic_write(out, base)
    print("\n=== anomaly verdict ===")
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


def _atomic_write(out: Path, data):
    tmp = out.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    os.replace(tmp, out)


if __name__ == "__main__":
    main()
