#!/usr/bin/env python3
"""Гровер-механика на движке и приговорка гибриду «классика сужает - квант добивает».

Оракул-уровень Grover на Statevector (qip_engine):
  оракул  - диагональный: помеченные амплитуды домножаются на -1
            (подмножество кандидатов C: интервал или одиночные метки
            в разных ветвях Lambda-орбиты);
  диффузия - точная D = 2|s><s| - I через psi -> 2*sum(psi)/N - psi.
Верифицируем теорию P(k) = sin^2((2k+1)*theta), sin^2(theta) = |C|/N,
пик на k* = pi/(4*theta) - 1/2, инвариантность к фрагментации меток
(K интервалов vs один интервал той же суммарной массы) и выборкой
через measure().

Экстраполяция на наши измеренные параметры (источники - не выдумка):
  results/h_curve_summary.json       - окно 1.44e12 при H=400k, гигантские
                                       шаги 3.6e6, 309918 шагов/с, wall 17.1 c
  results/phase_signature_135bit.json - маршрутизация по орбитам (3 ветви),
                                       mean_sig_ms 0.062
  results/qip_shor_ecdlp.json        - Шор: 2330 лог. кубитов цит. Roetteler
Считаем: сколько Гровер-итераций реально снимает наша классическая
фильтрация (sqrt(c) - ожидание), и где проходит граница классика/квант
W*(kappa_eff) при эффективной стоимости оракула kappa_eff секунд.

Запуск из корня репо: python3 experiments/qip_grover_hybrid.py
Выход: results/qip_grover_hybrid.json
"""
import argparse
import json
import math
import random
import time
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qip_engine import Statevector  # noqa: E402


# ---------- оракул-уровень Grover ----------

def grover_curve(n, marked, k_max, x_true):
    """Точная кривая успеха: P_marked(k) и |psi[x_true]|^2 для k=0..k_max.

    marked - np.array индексов помеченных кандидатов (включая x_true).
    Оракул: psi[marked] *= -1; диффузия: psi -> 2*sum(psi)/N - psi.
    """
    psi = np.full(n, 1.0 / math.sqrt(n), dtype=np.complex128)
    inv = 1.0 / math.sqrt(n)
    p_marked = np.empty(k_max + 1)
    p_true = np.empty(k_max + 1)
    mset = np.zeros(n, dtype=bool)
    mset[marked] = True
    for k in range(k_max + 1):
        p_marked[k] = float((np.abs(psi[mset]) ** 2).sum())  # масса на C
        p_true[k] = float(abs(psi[x_true]) ** 2)
        if k == k_max:
            break
        psi[mset] *= -1.0                       # оракул
        s = psi.sum()                           # диффузия 2|s><s| - I
        psi = 2.0 * s / n - psi
    return p_marked, p_true


def analyze_curve(p_marked, c_over_n):
    """Теория против замера: theta, k*, пик, k90."""
    theta = math.asin(math.sqrt(c_over_n))
    k_theory = math.pi / (4.0 * theta) - 0.5
    k_peak = int(np.argmax(p_marked))
    peak = float(p_marked[k_peak])
    k90 = int(np.argmax(p_marked >= 0.90)) if (p_marked >= 0.90).any() else None
    return {
        "theta": theta,
        "k theory peak": k_theory,
        "k measured peak": k_peak,
        "peak deviation": abs(k_peak - k_theory),
        "peak P marked": round(peak, 6),
        "k first >= 0.90": k90,
    }


# ---------- сценарии ----------

def run_case(n, c, mode, x_true=None, rng=None):
    """Один случай: c помеченных кандидатов, mode: interval|fragmented."""
    if x_true is None:
        x_true = rng.randrange(n)
    if mode == "interval":
        start = max(0, min(x_true - c // 2, n - c))
        marked = np.arange(start, start + c)
    else:  # fragmented: одиночные метки в разных ветвях орбиты
        step = n // 7
        marked = np.array([(x_true + i * step) % n for i in range(c)])
    k_theory = math.pi / (4.0 * math.asin(math.sqrt(c / n))) - 0.5
    k_max = int(math.ceil(k_theory * 1.4))
    t0 = time.perf_counter()
    p_marked, p_true = grover_curve(n, marked, k_max, x_true)
    dt = time.perf_counter() - t0
    res = {"n": n, "marked c": c, "mode": mode, "k max": k_max,
           "iterations seconds": round(dt, 2)}
    res.update(analyze_curve(p_marked, c / n))
    res["p marked at k90"] = (round(float(p_marked[res["k first >= 0.90"]]), 4)
                              if res["k first >= 0.90"] is not None else None)
    res["p true at peak"] = round(float(p_true[res["k measured peak"]]), 6)
    return res


def shots_check(n, c, shots, rng):
    """Выборка после k* итераций: доля попаданий в C и в x_true."""
    x_true = rng.randrange(n)
    start = (x_true - c // 2) % n
    marked = (start + np.arange(c)) % n          # интервал с заворотом
    k_theory = int(math.pi / (4.0 * math.asin(math.sqrt(c / n))) - 0.5)
    psi = np.full(n, 1.0 / math.sqrt(n), dtype=np.complex128)
    mset = np.zeros(n, dtype=bool)
    mset[marked] = True
    for _ in range(k_theory):
        psi[mset] *= -1.0
        psi = 2.0 * psi.sum() / n - psi
    w = n.bit_length() - 1
    hit, true_hits = 0, 0
    for _ in range(shots):
        sv = Statevector(w)
        sv.psi = psi.copy()
        idx, _ = sv.measure(rng)
        if mset[idx]:
            hit += 1
        if idx == x_true:
            true_hits += 1
    return {"n": n, "c": c, "k used": k_theory, "shots": shots,
            "hit C sampled": hit, "P C exact": round(float(p_marked_exact(n, mset, k_theory)), 4),
            "hit x sampled": true_hits, "P x exact": round(float(p_true_exact(n, mset, x_true, k_theory)), 4)}


def p_marked_exact(n, mset, k):
    psi = np.full(n, 1.0 / math.sqrt(n), dtype=np.complex128)
    for _ in range(k):
        psi[mset] *= -1.0
        psi = 2.0 * psi.sum() / n - psi
    return float((np.abs(psi[mset]) ** 2).sum())


def p_true_exact(n, mset, x, k):
    psi = np.full(n, 1.0 / math.sqrt(n), dtype=np.complex128)
    for _ in range(k):
        psi[mset] *= -1.0
        psi = 2.0 * psi.sum() / n - psi
    return float(abs(psi[x]) ** 2)


# ---------- реальные параметры ----------

def real_parameters():
    """Экстраполяция на измеренные константы H-кривой и маршрутизации."""
    hs = json.loads(Path("results/h_curve_summary.json").read_text())
    pt = hs["points"][3]                      # H=400k
    canon = pt["throughput"]["canonical H rows"]
    W = pt["window size"]                     # 1.44e12
    M = canon["table rows stored"]            # 400000
    gs = canon["giant steps"]                 # 3600000
    rate = canon["giant steps/second"]        # 309918
    search_s = canon["search seconds"]        # 11.6
    build_s = canon["build seconds"]          # 5.36
    wall_s = canon["wall seconds"]            # 17.1

    k_full = math.pi / 4.0 * math.sqrt(W)     # Гровер по всему окну
    rho_orbit = 3                             # маршрутизация по орбитам (замер)
    k_routed = math.pi / 4.0 * math.sqrt(W / rho_orbit)
    c_six = 2.75                              # кандидатов на попадание (замер)
    k_posthit = 0.0                           # после попадания квант не нужен

    # граница классика/квант: (W/M)/rate = (pi/4)*sqrt(W)*kappa
    # -> sqrt(W*) = (pi/4)*rate*M*kappa
    def crossover(kappa):
        return (math.pi / 4.0 * rate * M * kappa) ** 2

    # паззл 135: классика и Гровер
    W135 = 2 ** 135
    gs135 = W135 / M
    years135 = gs135 / rate / 3.156e7
    k135 = math.pi / 4.0 * math.sqrt(W135)

    return {
        "source": "results/h_curve_summary.json (H=400k, canonical H rows) + "
                  "results/phase_signature_135bit.json (orbit_mode)",
        "window W": W, "table rows M": M,
        "giant steps measured": gs, "steps per second": rate,
        "search seconds": search_s, "build seconds": build_s,
        "wall seconds": wall_s,
        "grover full window": {
            "oracle calls": round(k_full, 0),
            "ratio to classical giant steps": round(k_full / gs, 2),
            "note": "каждый вызов оракула = скалярное умножение EC в "
                    "суперпозиции (кси*10^6 логических вентилей); "
                    "классический шаг = одно сложение точек (~микросекунды)",
        },
        "hybrid signature routing (rho=3)": {
            "oracle calls": round(k_routed, 0),
            "iterations removed": round(100 * (1 - k_routed / k_full), 1),
            "reduction factor": round(math.sqrt(rho_orbit), 3),
            "classical routing cost ms": 0.062,
        },
        "hybrid post-hit (six-branch c=2.75)": {
            "oracle calls": k_posthit,
            "note": "после попадания таблицы кандидат ~2.75 скаляра - "
                    "проверяется классически, Гровер не нужен вовсе",
        },
        "crossover W*(kappa_eff) = ((pi/4)*rate*M*kappa)^2": {
            "kappa_eff 1 us": f"2^{math.log2(crossover(1e-6)):.1f}",
            "kappa_eff 1 ms": f"2^{math.log2(crossover(1e-3)):.1f}",
            "kappa_eff 1 s": f"2^{math.log2(crossover(1.0)):.1f}",
            "note": "окно, при котором Гровер догоняет наш классический "
                    "скан: левее - классика быстрее, правее - Гровер дешевле "
                    "по секундам, но оба времени при больших W астрономичны; "
                    "практичный квантовый путь - только Шор",
        },
        "puzzle 135 (W=2^135)": {
            "classical giant steps": f"{gs135:.3e}",
            "classical years (309918 шагов/с)": f"{years135:.2e}",
            "grover oracle calls": f"{k135:.3e}",
            "grover years (kappa_eff=1 s)": f"{k135 / 3.156e7:.2e}",
            "shor": "2330 лог. кубитов, ~1.26e11 Тоффоли "
                    "(Roetteler et al. 2017) - вне таблиц, вне поиска",
        },
    }


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20261110)
    ap.add_argument("--out", default="results/qip_grover_hybrid.json")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    report = {"description": "Гровер на движке: интервальный оракул, "
                             "фрагментация, гибрид с классической фильтрацией",
              "seed": args.seed, "cases": [], "shots validation": None,
              "real parameters": None}

    # w=18: полный набор сужений
    n18 = 2 ** 18
    for c, mode in ((1, "interval"), (2, "interval"), (3, "interval"),
                    (3, "fragmented"), (6, "interval"), (16, "interval")):
        r = run_case(n18, c, mode, rng=rng)
        report["cases"].append(r)
        print(f"[w=18 c={c:2d} {mode:10s}] k*={r['k theory peak']:7.1f} "
              f"замер={r['k measured peak']:4d} (откл {r['peak deviation']:.1f}) "
              f"P_peak={r['peak P marked']:.4f} k90={r['k first >= 0.90']} "
              f"за {r['iterations seconds']} с")

    # w=20: подтверждение масштабирования
    n20 = 2 ** 20
    for c in (1, 3):
        r = run_case(n20, c, "interval", rng=rng)
        r["note"] = "подтверждение 1/sqrt(c) на большем регистре"
        report["cases"].append(r)
        print(f"[w=20 c={c:2d} interval  ] k*={r['k theory peak']:7.1f} "
              f"замер={r['k measured peak']:4d} P_peak={r['peak P marked']:.4f} "
              f"за {r['iterations seconds']} с")

    # редукция итераций от фильтрации: замер против теории 1/sqrt(c)
    base = {int(math.log2(r["n"])): r["k measured peak"]
            for r in report["cases"] if r["marked c"] == 1}
    red = []
    for r in report["cases"]:
        if r["mode"] == "interval":
            w = int(math.log2(r["n"]))
            red.append({"c": r["marked c"], "w": w,
                        "reduction measured":
                            round(r["k measured peak"] / base[w], 3),
                        "reduction theory 1/sqrt(c)":
                            round(1.0 / math.sqrt(r["marked c"]), 3)})
    report["filtering reduction"] = red

    # выборочная валидация measure() против точных вероятностей
    sc = shots_check(2 ** 14, 3, 400, rng)
    report["shots validation"] = sc
    print(f"[shots w=14 c=3] попало в C: {sc['hit C sampled']}/400 "
          f"(точно {sc['P C exact']}), в x: {sc['hit x sampled']}/400 "
          f"(точно {sc['P x exact']})")

    # реальные параметры
    report["real parameters"] = real_parameters()
    rp = report["real parameters"]
    print(f"[real W=1.44e12] Гровер всё окно: {rp['grover full window']['oracle calls']:.0f} "
          f"оракулов ({rp['grover full window']['ratio to classical giant steps']}x "
          f"гигантских шагов); маршрутизация rho=3 снимает "
          f"{rp['hybrid signature routing (rho=3)']['iterations removed']}% итераций")

    out = Path(args.out)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"[saved] {out}")


if __name__ == "__main__":
    main()
