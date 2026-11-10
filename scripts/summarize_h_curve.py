#!/usr/bin/env python3
"""Сводная H-кривая: база H=100k (table_representative_sweep.json) плюс
точки 50k/200k/400k (h_curve_H*.json). Печатает таблицу и пишет
h_curve_summary.json с планами throughput (EC-ops, RSS, сборка)."""
import json
import sys
from pathlib import Path

ROOT = Path("/home/z/my-project/ECDSA")
POINTS = [50_000, 100_000, 200_000, 400_000]
FILES = {50_000: ROOT / "h_curve_H50000.json",
         100_000: ROOT / "table_representative_sweep.json",
         200_000: ROOT / "h_curve_H200000.json",
         400_000: ROOT / "h_curve_H400000.json"}

PLANS = ("continuous 3H rows", "expanded 3H rows", "canonical H rows",
         "expanded->canonical merge")


def main():
    summary = {"points": {}, "failures": []}
    print(f"{'H':>7} {'plan':<26} {'rows':>8} {'build s':>8} "
          f"{'search s':>9} {'steps/s':>9} {'EC mult':>8} {'EC add':>9} "
          f"{'total EC':>9} {'RSS MB':>7} {'wall s':>7}")
    for heads in POINTS:
        path = FILES[heads]
        if not path.exists():
            print(f"  [skip] {path.name} отсутствует")
            continue
        data = json.loads(path.read_text())
        if data.get("failures"):
            summary["failures"].append({heads: data["failures"]})
        summary["points"][str(heads)] = {}
        for plan in PLANS:
            run = data["throughput"].get(plan)
            if run is None:
                continue
            rss_mb = round(run["peak rss kb"] / 1024, 1)
            row = {
                "rows": run["table rows stored"],
                "build seconds": run["build seconds"],
                "search seconds": run["search seconds"],
                "giant steps/second": run["giant steps/second"],
                "EC scalar multiplications": run["EC scalar multiplications"],
                "EC point additions": run["EC point additions"],
                "total EC operations": run["total EC operations"],
                "peak rss mb": rss_mb,
                "wall seconds": run["wall seconds"],
                "canonical rows merged": run.get("canonical rows merged", 0),
            }
            summary["points"][str(heads)][plan] = row
            print(f"{heads:>7} {plan:<26} {run['table rows stored']:>8} "
                  f"{run['build seconds']:>8.2f} {run['search seconds']:>9.2f} "
                  f"{run['giant steps/second']:>9.0f} "
                  f"{run['EC scalar multiplications']:>8} "
                  f"{run['EC point additions']:>9} "
                  f"{run['total EC operations']:>9} {rss_mb:>7} "
                  f"{run['wall seconds']:>7.2f}")
        # экономика: агрегированные поля плана
        econ = data.get("economy", {})
        for name, plan_data in econ.items():
            if not isinstance(plan_data, dict) or "runs" not in plan_data:
                continue
            n = len(plan_data["runs"])
            ok = sum(1 for r in plan_data["runs"] if r.get("correct"))
            entry = {
                "targets": n, "correct": ok,
                "avg exact checks per hit":
                    plan_data.get("avg exact checks per hit"),
                "avg total EC ops per run":
                    plan_data.get("avg total EC ops per run"),
                "partner recoveries": plan_data.get("partner recoveries"),
                "signature pair mismatches":
                    plan_data.get("signature pair mismatches"),
            }
            summary["points"][str(heads)][f"economy {name}"] = entry
            print(f"        economy {name}: {ok}/{n} correct, "
                  f"exact/hit {entry['avg exact checks per hit']}, "
                  f"EC-ops/run {entry['avg total EC ops per run']}")
    if summary["failures"]:
        print("!! ЕСТЬ FAILURES:", summary["failures"], file=sys.stderr)
    out = ROOT / "h_curve_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(f"[saved] {out}")


if __name__ == "__main__":
    main()
