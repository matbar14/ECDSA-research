#!/usr/bin/env python3
"""Огибающая statevector-симуляции: сколько кубитов помещается в память.

Вектор состояния: 2^w комплексных амплитуд по 16 байт (complex128) =>
w кубитов стоит 16*2^w байт, каждый следующий кубит удваивает память.
Замеряем время типовых слоёв (1-кубитный слой, CNOT-слой, QFT-16 гейтами)
на доступной машине и экстраполируем на целевой бюджет памяти.

Запуск из корня репо: python3 experiments/qip_statevector.py [--max-gib 0.5]
Выход: results/qip_statevector_envelope.json
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qip_engine import Statevector  # noqa: E402


def bench_width(w, rng):
    sv = Statevector(w)
    res = {"qubits": w, "statevector MB": 16.0 * 2 ** w / 2 ** 20}

    # 1-кубитный слой: H на каждом кубите
    t0 = time.perf_counter()
    for q in range(w):
        sv.h(q)
    res["h layer ms"] = round((time.perf_counter() - t0) * 1000, 2)

    # CNOT-слой: цепочка CX
    t0 = time.perf_counter()
    for q in range(w - 1):
        sv.cx(q, q + 1)
    res["cx layer ms"] = round((time.perf_counter() - t0) * 1000, 2)

    # QFT-16 гейтами на младших 16 кубитах (если помещаются)
    if w >= 16:
        t0 = time.perf_counter()
        sv.qft_pow2(list(range(15, -1, -1)))
        res["qft16 gates ms"] = round((time.perf_counter() - t0) * 1000, 2)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-gib", type=float, default=0.5,
                    help="бюджет памяти на вектор состояния в GiB")
    ap.add_argument("--target-gib", type=float, default=32.0,
                    help="бюджет целевой машины для экстраполяции")
    args = ap.parse_args()

    q_here = int(np.floor(np.log2(args.max_gib * 2 ** 30 / 16)))
    q_target = int(np.floor(np.log2(args.target_gib * 2 ** 30 / 16)))
    widths = sorted({q for q in (20, 22, 24, 25, q_here) if q <= q_here})

    print(f"бюджет здесь {args.max_gib} GiB -> максимум "
          f"{q_here} кубитов; целевой {args.target_gib} GiB -> "
          f"{q_target} кубитов")
    runs = []
    for w in widths:
        run = bench_width(w, None)
        runs.append(run)
        print(f"[w={w:2d}] {run['statevector MB']:8.1f} MB | "
              f"H-слой {run['h layer ms']:9.2f} ms | "
              f"CX-слой {run['cx layer ms']:9.2f} ms | "
              f"QFT16 {run.get('qft16 gates ms', float('nan')):10.2f} ms")

    # линейная в 2^w экстраполяция QFT16 на целевой максимум
    pts = [(r["qubits"], r["qft16 gates ms"]) for r in runs
           if "qft16 gates ms" in r]
    if len(pts) >= 2:
        (w1, t1), (w2, t2) = pts[0], pts[-1]
        scale = 2.0 ** ((q_target - w2) / (w2 - w1))
        qft_target = t2 * scale
    else:
        qft_target = None

    envelope = {
        "statevector bytes per qubit": 16,
        "qubits max here (GiB=%.2f)" % args.max_gib: q_here,
        "qubits max target (GiB=%.0f)" % args.target_gib: q_target,
        "bench": runs,
        "qft16 extrapolation": {
            "target qubits": q_target,
            "seconds": round(qft_target / 1000, 1) if qft_target else None,
            "note": "линейная в 2^w экстраполяция по замерам здесь; "
                    "на 8-ядерной машине ожидание то же (однопоточный "
                    "NumPy) +- кэш-эффекты",
        },
    }
    out = Path("results/qip_statevector_envelope.json")
    out.write_text(json.dumps(envelope, indent=2, ensure_ascii=False) + "\n")
    print(f"[saved] {out}")
    if qft_target:
        print(f"QFT16 гейтами на {q_target} кубитах ~ "
              f"{qft_target / 1000:.1f} с (экстраполяция)")


if __name__ == "__main__":
    main()
