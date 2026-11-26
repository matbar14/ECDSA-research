#!/usr/bin/env python3
"""Render the delta-sweep results as a self-contained HTML report.

Reads results/angle_table_delta_sweep.json and writes
figures/angle_delta_sweep.html: a P(exact) heatmap over (log2 W, delta),
the wall curve P(exact single) vs delta for all three arms with theory
overlays, and the memory/verification summary.
"""

import json
import math
import os
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results", "angle_table_delta_sweep.json")
OUT = os.path.join(HERE, "..", "figures", "angle_delta_sweep.html")

ARM_LABEL = {
    "pair": "pair (проверенный 2-угловой fingerprint)",
    "orbit": "orbit (GLV-канонизация, λ_eff = 3W/2^B)",
    "sig": "sig (пара B−2 бита + 2 бита фазовой подписи)",
}
ARM_SHIFT = {"pair": 0.0, "orbit": 0.0, "sig": 2.0}


def color(p):
    """Blue (fail) -> yellow (wall) -> green (perfect)."""
    if p >= 0.999:
        return "var(--success)"
    t = min(max(p, 0.0), 1.0)
    if t < 0.5:
        a, b, u = (30, 60, 120), (180, 150, 40), t / 0.5
    else:
        a, b, u = (180, 150, 40), (40, 140, 80), (t - 0.5) / 0.5
    r = int(a[0] + (b[0] - a[0]) * u)
    g = int(a[1] + (b[1] - a[1]) * u)
    bl = int(a[2] + (b[2] - a[2]) * u)
    return f"rgb({r},{g},{bl})"


def svg_polyline(points, x0, y0, sx, sy):
    pts = " ".join(f"{x0 + x * sx:.1f},{y0 - y * sy:.1f}" for x, y in points)
    return f'<polyline points="{pts}" fill="none" stroke-width="2"/>'


def main():
    with open(RESULTS) as handle:
        data = json.load(handle)
    cells = data["cells"]
    log2s = sorted({c["log2_w"] for c in cells})
    deltas = sorted({c["delta"] for c in cells})

    # --- heatmap: pair arm, P(exact) over (L, delta) ---
    grid = {(c["log2_w"], c["delta"]): c for c in cells if c["arm"] == "pair"}
    heat_rows = []
    for delta in reversed(deltas):
        tds = [f'<th class="rl">δ={delta:+d}</th>']
        for lg in log2s:
            p = grid[(lg, delta)]["p_exact_single"]
            tds.append(f'<td style="background:{color(p)}">{p:.2f}</td>')
        heat_rows.append("<tr>" + "".join(tds) + "</tr>")
    heat_head = "<tr><th></th>" + "".join(
        f"<th>{lg}</th>" for lg in log2s) + "</tr>"
    heat = ("<table class='heat'>" + heat_head + "".join(heat_rows) + "</table>")

    # --- wall curves: mean over L per (arm, delta) + theory ---
    curves = {}
    for arm in ("pair", "orbit", "sig"):
        sel = [c for c in cells if c["arm"] == arm]
        curves[arm] = [
            (d,
             statistics.mean(c["p_exact_single"] for c in sel if c["delta"] == d))
            for d in deltas
        ]
    W_, H_ = 640, 360
    x0, y0 = 56, 20
    sx = (W_ - x0 - 16) / (deltas[-1] - deltas[0])
    sy = (H_ - y0 - 44) / 1.0
    frame = (f'<svg viewBox="0 0 {W_} {H_}" class="chart">'
             f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0+sy}" class="ax"/>'
             f'<line x1="{x0}" y1="{y0+sy}" x2="{x0+sx*(deltas[-1]-deltas[0])}" '
             f'y2="{y0+sy}" class="ax"/>')
    for gy in (0.25, 0.5, 0.75, 1.0):
        yy = y0 + sy * (1 - gy)
        frame += (f'<line x1="{x0}" y1="{yy}" x2="{W_-16}" y2="{yy}" class="grid"/>'
                  f'<text x="{x0-6}" y="{yy+4}" class="lbl" text-anchor="end">'
                  f'{gy:.2f}</text>')
    for d in deltas:
        xx = x0 + sx * (d - deltas[0])
        frame += f'<text x="{xx}" y="{H_-24}" class="lbl" text-anchor="middle">{d:+d}</text>'
    frame += (f'<text x="{W_/2}" y="{H_-6}" class="lbl" text-anchor="middle">'
              f'δ = B − ⌈log₂W⌉</text>'
              f'<text x="14" y="{y0+sy/2}" class="lbl" transform="rotate(-90 14 {y0+sy/2})" '
              f'text-anchor="middle">P(exact), 200 запросов</text>')
    palette = {"pair": "var(--chart-1)", "orbit": "var(--chart-2)",
               "sig": "var(--chart-3)"}
    tpalette = {"pair": "var(--chart-4)", "orbit": "var(--chart-5)",
                "sig": "var(--chart-6)"}
    for arm in ("pair", "orbit", "sig"):
        frame += (f'<g stroke="{palette[arm]}">'
                  + svg_polyline(curves[arm], x0, y0 + sy,
                                 sx, sy) + "</g>")
        shift = ARM_SHIFT[arm]
        lam = 3.0 if arm == "orbit" else 1.0
        # pair: exp(-2^-d); orbit: exp(-3*2^-d) (wall right by log2 3);
        # sig: exp(-2^-(d-2)) = pair two deltas right (no signature entropy).
        theory = [(d, math.exp(-lam * 2.0 ** (-(d - shift))))
                  for d in deltas]
        frame += (f'<g stroke="{tpalette[arm]}" stroke-dasharray="5 4">'
                  + svg_polyline(theory, x0, y0 + sy, sx, sy) + "</g>")
    frame += "</svg>"
    legend = "".join(
        f"<span><i style='background:{palette[a]}'></i>{ARM_LABEL[a]} "
        f"<i class='dash' style='background:{tpalette[a]}'></i>теория</span>"
        for a in ("pair", "orbit", "sig"))

    # --- bucket law check table ---
    bucket_rows = []
    for delta in deltas:
        row = [f"<tr><th>δ={delta:+d}</th>"]
        for arm in ("pair", "orbit", "sig"):
            sel = [c for c in cells if c["arm"] == arm and c["delta"] == delta]
            row.append(f"<td>{statistics.mean(c['mean_bucket_query'] for c in sel):.3f}</td>")
        row.append(f"<td>{1 + 2.0**-delta:.3f}</td>")
        row.append(f"<td>{1 + 3 * 2.0**-delta:.3f}</td>")
        row.append(f"<td>{1 + 4 * 2.0**-delta:.3f}</td>")
        bucket_rows.append("".join(row) + "</tr>")

    # --- memory / verification table at the largest window ---
    big = max(log2s)
    mem_rows = []
    for delta in (2, 4, 6):
        c = next(c for c in cells
                 if c["log2_w"] == big and c["arm"] == "pair" and c["delta"] == delta)
        o = next(c for c in cells
                 if c["log2_w"] == big and c["arm"] == "orbit" and c["delta"] == delta)
        mem_rows.append(
            f"<tr><th>{big} (2^{big})</th><th>δ={delta:+d}</th>"
            f"<td>{c['packed_bytes']/1e9:.2f} ГБ</td>"
            f"<td>{o['packed_bytes']/1e9:.2f} ГБ</td>"
            f"<td>{o['packed_bytes_per_covered_scalar']/c['packed_bytes_per_covered_scalar']:.2f}×</td>"
            f"<td>{c['verify_s_per_query']*1e3:.2f} мс</td>"
            f"<td>{o['verify_s_per_query']*1e3:.2f} мс</td></tr>")

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
body {{ font-family: var(--font-sans); color: var(--foreground);
       max-width: 1100px; margin: 0 auto; padding: 8px 4px; }}
h2 {{ font-size: 17px; margin: 22px 0 8px; }}
p.note {{ color: var(--muted-foreground); font-size: 13px; margin: 4px 0 10px; }}
table {{ border-collapse: collapse; font-size: 12px; }}
th, td {{ border: 1px solid var(--border); padding: 3px 7px; text-align: center; }}
th.rl {{ text-align: right; white-space: nowrap; }}
td {{ min-width: 44px; color: var(--foreground); }}
table.heat td {{ font-size: 10px; }}
.axis {{ stroke: var(--muted-foreground); stroke-width: 1; }}
.grid {{ stroke: var(--border); stroke-width: 1; }}
.lbl {{ fill: var(--muted-foreground); font-size: 11px; }}
.legend {{ display: flex; gap: 18px; font-size: 13px; margin: 8px 0; flex-wrap: wrap; }}
.legend i {{ display: inline-block; width: 16px; height: 3px; vertical-align: middle; margin-right: 5px; }}
.legend i.dash {{ background: repeating-linear-gradient(90deg, currentColor 0 5px, transparent 5px 9px); }}
.kv {{ color: var(--muted-foreground); font-size: 13px; }}
code {{ background: var(--code-background); color: var(--code-foreground);
        padding: 1px 5px; border-radius: 4px; font-size: 12px; }}
</style></head><body>

<h2>P(exact) по сетке (log₂W, δ), плечо pair — изолинии вертикальны</h2>
<p class="note">Столбцы: log₂W от {log2s[0]} до {log2s[-1]} (W = 2^L, окно [2^134, 2^134+W)).
Строки: δ = B − ⌈log₂W⌉. В каждой ячейке 200 случайных запросов.
P(exact) не зависит от W при фиксированном δ — качество определяется только
занятостью бакета λ = 2^{{−δ}}.</p>
{heat}

<h2>Стена качества: P(exact) по δ для трёх плеч</h2>
<div class="legend">{legend}</div>
{frame}
<p class="note">Сплошные линии — среднее по всей оси log₂W (разброс по L виден на тепловой карте),
пунктир — теория P = exp(−λ_eff·2^{{−δ}}). pair: λ_eff = 2^{{−δ}} (идеальный Пуассон,
poisson-отношение занятости 0.98–1.01 на всех ячейках). orbit: λ_eff = 3·2^{{−δ}} —
стена сдвинута вправо на log₂3 ≈ 1.58 бита, зато каждая запись покрывает до 6 скаляров орбиты.
sig: совпадает с pair при δ−2 — два бита (d01,d12) детерминированы уже сохранёнными
старшими битами и добавляют ноль энтропии; гипотеза «сдвиг стены влево на ~1.5 бита»
в табличной модели опровергнута, стена уходит на 2 бита <b>вправо</b>.</p>

<h2>Закон бакетов (среднее по L): измерение vs предсказание</h2>
<table>
<tr><th></th><th>pair</th><th>orbit</th><th>sig</th>
<th>pair = 1+2^{{−δ}}</th><th>orbit = 1+3·2^{{−δ}}</th><th>sig = 1+4·2^{{−δ}}</th></tr>
{''.join(bucket_rows)}
</table>
<p class="note">sig = pair при δ−2: средний бакет 1+4·2^{{−δ}} = (1+2^{{−(δ−2)}}).</p>

<h2>Память и отложенная верификация при W = 2^{big}</h2>
<table>
<tr><th>W</th><th>δ</th><th>pair, упаковано</th><th>orbit, упаковано</th>
<th>orbit / pair на покрываемый скаляр</th><th>pair, вериф./запрос</th><th>orbit, вериф./запрос</th></tr>
{''.join(mem_rows)}
</table>
<p class="note">Упаковка: B бит fingerprint + L бит offset (+2 бита индекса орбиты для orbit),
без выравнивания. Верификация — стоимость точной EC-перепроверки кандидата
({data['verify_cost_s_per_candidate']*1e3:.3f} мс/кандидат, измерено) × средний бакет запроса.
Хвост занятости (p99/max по занятым ключам) и все сырые ячейки: results/angle_table_delta_sweep.json.</p>

</body></html>"""
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as handle:
        handle.write(html)
    print(f"written {OUT}")


if __name__ == "__main__":
    main()
