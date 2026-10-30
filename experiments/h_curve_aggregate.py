#!/usr/bin/env python3
"""Агрегация H-кривой (50k/100k/200k/400k) и график память/сборка/пропускная
способность. Вход: h_curve_<H>.json из table_representative_sweep.py.
Выход: h_curve_summary.json + h_curve.png."""

import json
from pathlib import Path

import matplotlib.font_manager as fm
fm.fontManager.addfont('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')

import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parent
HEADS = [50_000, 100_000, 200_000, 400_000]
SOURCE = {50_000: "results/h_curve_50k.json", 100_000: "results/h_curve_100k.json",
          200_000: "results/h_curve_200k.json", 400_000: "results/h_curve_400k.json"}
VARIANTS = [
    ("continuous 3H rows", "continuous (3H строк)", "#8a8a8a", "o"),
    ("expanded 3H rows", "expanded (3H строк)", "#4C72B0", "s"),
    ("canonical H rows", "canonical (H строк)", "#55A868", "^"),
    ("expanded->canonical merge", "expanded->canonical merge", "#C44E52", "D"),
]
ECONOMY = [
    ("canonical six-branch", "six-branch (без маршрутизации)", "#8a8a8a", "o"),
    ("canonical routed", "canonical + сигнатурная маршрутизация", "#55A868", "^"),
]

summary = {"points": [], "economy": {}, "notes": {}}
data = {}
for H in HEADS:
    path = ROOT / SOURCE[H]
    r = json.loads(path.read_text())
    data[H] = r
    tp = r["throughput"]
    point = {"heads": H, "window size": r["window size"], "throughput": {}}
    for key, *_ in VARIANTS:
        v = tp[key]
        point["throughput"][key] = {
            "table rows stored": v["table rows stored"],
            "giant steps": v["giant steps"],
            "giant steps/second": v["giant steps/second"],
            "EC scalar multiplications": v["EC scalar multiplications"],
            "EC point additions": v["EC point additions"],
            "total EC operations": v["total EC operations"],
            "build seconds": v["build seconds"],
            "search seconds": v["search seconds"],
            "wall seconds": v["wall seconds"],
            "peak rss mb": round(v["peak rss kb"] / 1024, 1),
            "table bytes per row": round(
                v["peak rss kb"] * 1024 / v["table rows stored"], 1),
        }
    summary["points"].append(point)

# экономика: свертка по H (маршрутизация не должна зависеть от H)
for H in HEADS:
    for key, _label, _c, _m in ECONOMY:
        e = data[H]["economy"][key]
        acc = summary["economy"].setdefault(key, {"avg exact per hit": [],
                                                  "all correct": []})
        acc["avg exact per hit"].append(e["avg exact checks per hit"])
        acc["all correct"].append(e["all correct"])
for key, acc in summary["economy"].items():
    vals = acc["avg exact per hit"]
    acc["avg exact per hit by H"] = {str(h): v
                                     for h, v in zip(HEADS, vals)}
    acc["mean over H"] = round(sum(vals) / len(vals), 3)
    acc["all correct everywhere"] = all(acc["all correct"])
    acc.pop("avg exact per hit")
    acc.pop("all correct")

can = {p["heads"]: p["throughput"]["canonical H rows"] for p in summary["points"]}
exp = {p["heads"]: p["throughput"]["expanded 3H rows"] for p in summary["points"]}
summary["notes"] = {
    "rss compression canonical/expanded": {
        h: round(exp[h]["peak rss mb"] / can[h]["peak rss mb"], 2) for h in HEADS},
    "build speedup canonical/expanded": {
        h: round(exp[h]["build seconds"] / can[h]["build seconds"], 2) for h in HEADS},
    "giant steps canonical / expanded": {
        h: round(can[h]["giant steps"] / exp[h]["giant steps"], 1) for h in HEADS},
    "steps-per-second range over all variants and H": [
        min(p["throughput"][k]["giant steps/second"]
            for p in summary["points"] for k, *_ in VARIANTS),
        max(p["throughput"][k]["giant steps/second"]
            for p in summary["points"] for k, *_ in VARIANTS)],
}
(ROOT / "results/h_curve_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary["notes"], indent=2, ensure_ascii=False))

# ---------------- график 2x2 ----------------
fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0), constrained_layout=True)
Hk = [h / 1000 for h in HEADS]

ax = axes[0][0]
for key, label, color, marker in VARIANTS:
    y = [p["throughput"][key]["peak rss mb"] for p in summary["points"]]
    ax.plot(Hk, y, marker=marker, color=color, label=label, lw=1.8, ms=6)
ax.set_title("Пиковая память процесса (RSS)")
ax.set_ylabel("МБ")
ax.legend(fontsize=8.5, loc="upper left")
ax.grid(alpha=0.3)

ax = axes[0][1]
for key, label, color, marker in VARIANTS:
    y = [p["throughput"][key]["build seconds"] for p in summary["points"]]
    ax.plot(Hk, y, marker=marker, color=color, label=label, lw=1.8, ms=6)
ax.set_title("Сборка таблицы")
ax.set_ylabel("секунды")
ax.grid(alpha=0.3)

ax = axes[1][0]
for key, label, color, marker in VARIANTS:
    y = [p["throughput"][key]["giant steps/second"] for p in summary["points"]]
    ax.plot(Hk, y, marker=marker, color=color, label=label, lw=1.8, ms=6)
ax.set_title("Пропускная способность giant-фазы")
ax.set_ylabel("giant-шагов/сек")
ax.set_ylim(250_000, 400_000)
ax.grid(alpha=0.3)

ax = axes[1][1]
for key, label, color, marker in ECONOMY:
    y = [summary["economy"][key]["avg exact per hit by H"][str(h)]
         for h in HEADS]
    ax.plot(Hk, y, marker=marker, color=color, label=label, lw=1.8, ms=6)
ax.set_title("Экономика точных проверок (12 целей на точку)")
ax.set_ylabel("точных EC-проверок на hit")
ax.set_ylim(0, 3.4)
ax.legend(fontsize=8.5, loc="center right")
ax.grid(alpha=0.3)

for ax in (axes[0][0], axes[0][1], axes[1][0], axes[1][1]):
    ax.set_xlabel("H, голов орбит (тысяч)")
    ax.set_xticks(Hk)

fig.suptitle("H-кривая таблицы представителей полного диапазона "
             "(окно S = (3H)$^2$, start = 2$^{134}$, 2 потока)",
             fontsize=13)
fig.savefig(ROOT / "figures/h_curve.png", dpi=150)
print("saved h_curve.png")
