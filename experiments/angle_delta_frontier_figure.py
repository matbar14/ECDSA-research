#!/usr/bin/env python3
"""Render the cost-vs-memory frontier report for the fingerprint table.

Two explicitly separated operating modes (the memo question: "is 100/100
the right target?"):

  filter mode   - the table pre-filters candidates for a downstream search
                  (BSGS/kangaroo).  Per-query success e^{-2^-delta}, misses
                  fall back; the useful frontier is delta = 0..2 where the
                  per-query verify cost stays ~1 ms and the fallback share
                  22-63% is absorbed by the search that runs anyway.

  standalone    - the table must answer every one of Q queries.  The wall
                  family: delta_req(Q, P_target) = log2(Q / (-ln P_target)),
                  so delta grows with log2(Q) and the 100/100 criterion
                  (P -> 1) never closes: delta_req(P=1) = infinity.

All curves are analytic (validated by the delta sweep); the only measured
constants are the verify costs (0.776 ms reference path, 0.341 ms with the
precomputed [2^j]G + coincurve reconstruction, both measured and checked
for exact-invariant equivalence).
"""

import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "angle_delta_frontier.html")

V_REF_MS = 0.776    # measured: reference deferred verify per candidate
V_FAST_MS = 0.341   # measured: precomputed [2^j]G + coincurve + 2 lambdas
LAMBDA = {"pair": 1.0, "orbit": 3.0}


def delta_req(q, p_target):
    return math.log2(q / (-math.log(p_target)))


def main():
    # load measured anchors from the delta sweep
    with open(os.path.join(HERE, "..", "results",
                           "angle_table_delta_sweep.json")) as handle:
        sweep = json.load(handle)
    anchors = {}
    for cell in sweep["cells"]:
        if cell["arm"] == "pair" and cell["log2_w"] == 27:
            anchors[cell["delta"]] = cell["p_exact_single"]

    # --- panel A: wall family delta_req(Q, P_target) ---
    qs = [1, 3, 10, 30, 100, 300, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7]
    targets = [("P=0.50", 0.50, "var(--chart-1)"),
               ("P=0.95", 0.95, "var(--chart-2)"),
               ("P=0.99", 0.99, "var(--chart-3)"),
               ("P=1-1e-3", 0.999, "var(--chart-4)"),
               ("P=1-1e-6", 0.999999, "var(--chart-5)")]
    W_, H_ = 660, 380
    x0, y0 = 60, 22
    lq_min, lq_max = 0.0, 7.0
    d_min, d_max = -2.0, 26.0
    sx = (W_ - x0 - 16) / (lq_max - lq_min)
    sy = (H_ - y0 - 48) / (d_max - d_min)
    svg = (f'<svg viewBox="0 0 {W_} {H_}" class="chart">'
           f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0+sy*(d_max-d_min)}" class="ax"/>'
           f'<line x1="{x0}" y1="{y0+sy*(d_max-d_min)}" '
           f'x2="{x0+sx*(lq_max-lq_min)}" y2="{y0+sy*(d_max-d_min)}" class="ax"/>')
    for dy in range(0, 27, 4):
        yy = y0 + sy * (d_max - dy)
        svg += (f'<line x1="{x0}" y1="{yy}" x2="{W_-16}" y2="{yy}" class="grid"/>'
                f'<text x="{x0-6}" y="{yy+4}" class="lbl" text-anchor="end">{dy}</text>')
    for lg in range(0, 8):
        xx = x0 + sx * lg
        svg += (f'<text x="{xx}" y="{H_-26}" class="lbl" text-anchor="middle">'
                f'10^{lg}</text>')
    for label, p, col in targets:
        pts = []
        for q in qs:
            lgq = math.log10(q)
            if lgq < lq_min or lgq > lq_max:
                continue
            d = delta_req(q, p)
            if d < d_min or d > d_max:
                continue
            pts.append((x0 + sx * (lgq - lq_min), y0 + sy * (d_max - d)))
        poly = " ".join(f"{a:.1f},{b:.1f}" for a, b in pts)
        svg += f'<polyline points="{poly}" fill="none" stroke="{col}" stroke-width="2"/>'
        lx, ly = pts[-1]
        svg += (f'<text x="{lx+6}" y="{ly+4}" class="lbl" fill="{col}">{label}</text>')
    # measured per-query anchors (Q=1 cells): P at delta from the sweep
    for d in (0, 2, 4, 6):
        xx = x0 + sx * (0 - lq_min)
        yy = y0 + sy * (d_max - d)
        svg += f'<circle cx="{xx}" cy="{yy}" r="4" fill="var(--foreground)"/>'
        svg += (f'<text x="{xx+8}" y="{yy-6}" class="lbl">δ={d}: '
                f'P={anchors[d]:.2f}</text>')
    svg += (f'<text x="{W_/2}" y="{H_-6}" class="lbl" text-anchor="middle">'
            f'Q запросов (log₁₀)</text>'
            f'<text x="14" y="{y0+sy*(d_max-d_min)/2}" class="lbl" '
            f'transform="rotate(-90 14 {y0+sy*(d_max-d_min)/2})" '
            f'text-anchor="middle">требуемый δ</text></svg>')
    svg += ('<p class="note">Точки — измерение (Q=1, W=2^27); кривые — '
            'δ_req(Q,P) = log₂(Q/(−ln P)). «Стена» — семейство, а не число: '
            'P=1 недостижима (δ→∞), 100/100 на Q запросов живёт на кривой '
            'P=1−1/Q.</p>')

    # --- panel B: two modes at fixed W ---
    rows = []
    for lg in (27, 30, 40):
        w = 1 << lg
        # filter mode, delta = 2
        d = 2
        mem = w * (2 * lg + d) / 8
        mem_h = (f"{mem/1e12:.2f} ТБ" if mem >= 1e12 else f"{mem/1e9:.2f} ГБ")
        rows.append(
            f"<tr><th rowspan='2'>2^{lg}</th><td>filter, δ={d}</td>"
            f"<td>{mem_h}</td><td>{(1+2.0**-d)*V_FAST_MS:.2f} мс</td>"
            f"<td>{math.exp(-2.0**-d):.2f} на запрос</td>"
            f"<td>{(1-math.exp(-2.0**-d))*100:.0f}% fallback</td>"
            f"<td>filter: fallback съедает поиск, но он и так работал бы</td></tr>")
        # standalone: Q=1e4 at P=0.95
        q, pt = 1e4, 0.95
        d = delta_req(q, pt)
        mem = w * (2 * lg + d) / 8
        mem_h = (f"{mem/1e12:.2f} ТБ" if mem >= 1e12 else f"{mem/1e9:.2f} ГБ")
        total = q * (1 + 2.0**-d) * V_FAST_MS / 1000
        rows.append(
            f"<tr><td>standalone, Q=10⁴, P=0.95</td>"
            f"<td>{mem_h}</td><td>{(1+2.0**-d)*V_FAST_MS:.2f} мс</td>"
            f"<td>{pt:.2f} на все 10⁴</td><td>—</td>"
            f"<td>total verify ≈ {total:.1f} с</td></tr>")

    # --- panel C: memory per record law ---
    mem_rows = []
    for lg in (27, 30, 34, 40, 67):
        w = 1 << lg
        for d, tag in ((2, "filter"), (6, "P≈0.98/query"),
                       (math.ceil(delta_req(1e4, 0.95)), "standalone 10⁴@0.95")):
            mem = w * (2 * lg + d) / 8
            human = (f"{mem/1e21:.2f} ЗБ" if mem >= 1e21 else
                     f"{mem/1e18:.2f} ЭБ" if mem >= 1e18 else
                     f"{mem/1e12:.2f} ТБ" if mem >= 1e12 else (
                         f"{mem/1e9:.2f} ГБ" if mem >= 1e9 else f"{mem/1e6:.0f} МБ"))
            mem_rows.append(
                f"<tr><td>2^{lg}</td><td>δ={d}</td><td>{2*lg+d} бит/запись</td>"
                f"<td>{human}</td><td>{tag}</td></tr>")

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
body {{ font-family: var(--font-sans); color: var(--foreground);
       max-width: 1000px; margin: 0 auto; padding: 8px 4px; }}
h2 {{ font-size: 17px; margin: 22px 0 8px; }}
p.note {{ color: var(--muted-foreground); font-size: 13px; margin: 4px 0 10px; }}
table {{ border-collapse: collapse; font-size: 13px; }}
th, td {{ border: 1px solid var(--border); padding: 4px 9px; text-align: left; }}
th {{ background: var(--secondary); }}
.axis {{ stroke: var(--muted-foreground); }}
.grid {{ stroke: var(--border); }}
.lbl {{ fill: var(--muted-foreground); font-size: 11px; }}
.chart {{ width: 100%; height: auto; }}
code {{ background: var(--code-background); color: var(--code-foreground);
        padding: 1px 5px; border-radius: 4px; font-size: 12px; }}
</style></head><body>

<h2>Семейство стен: δ_req(Q, P_target)</h2>
{svg}

<h2>Два режима таблицы (главный практический вывод)</h2>
<table>
<tr><th>W</th><th>режим</th><th>память</th><th>verify/запрос</th>
<th>P успеха</th><th>fallback</th><th>комментарий</th></tr>
{''.join(rows)}
</table>
<p class="note">verify/запрос = (1+2^{{−δ}})·{V_FAST_MS:.3f} мс (fast-verify,
измерено; reference-путь {V_REF_MS:.3f} мс — 2.3× медленнее). Filter-режим при
δ=2 доминирует по памяти: стена по критерию 100/100 отодвинула бы бюджет на
~log₂Q бит, память ×2^{{log₂Q·..}} растёт вместе с Q, а filter-режим платит
22% fallback вместо этого.</p>

<h2>Память на запись: M = W·(2L+δ)/8 по закону таблицы</h2>
<table>
<tr><th>W</th><th>δ</th><th>плотность</th><th>упаковано</th><th>режим</th></tr>
{''.join(mem_rows)}
</table>
<p class="note">W≈2^30–2^32 — одна машина (десятки ГБ), W≈2^35 — кластер,
W≥2^40 — распределённое хранилище, W=2^67 (BSGS-таргет) — ~1 ЭБ: структура
масштабируется, железо нет. Измеренная точка: W=2^27, δ=+2 → 0.94 ГБ
(results/angle_table_delta_sweep.json), совпадает с законом M = W·(2L+δ)/8.</p>

<h2>Ускорение верификации (пункт «precomputed [2^j]G»)</h2>
<p class="note">Reference-путь 0.776 мс/кандидат (cryptography derive 706 мкс +
инварианты). Fast-путь: precomputed [2^j]G (j=0..39) + [2^134]G + сложения
coincurve по битам offset + 2 инварианта = 0.341 мс (2.3×); EC-часть отдельно
706 → 161 мкс (4.4×). Эквивалентность fast-пути reference-верификации
проверена на 200 случайных кандидатах по точным инвариантам. Pippenger-батчинг
неприменим (база фиксированная, это не multi-scalar) — выигрыш только общий
на все кандидаты precomputed-таблицы. Ускорение верификации не двигает δ-стену
(вероятности те же), но дешевит fallback в filter-режиме.</p>

</body></html>"""
    with open(OUT, "w") as handle:
        handle.write(html)
    print(f"written {OUT}")


if __name__ == "__main__":
    main()
