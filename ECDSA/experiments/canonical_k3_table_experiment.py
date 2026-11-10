#!/usr/bin/env python3
"""Каноническая K3-таблица: одна строка на орбиту вместо трёх.

K3 = x^3 один и тот же для трёх x-представителей одной C6-орбиты
(x, beta*x, beta^2*x), поэтому bucket полной орбиты содержит три строки
(j, lambda*j, lambda^2*j). Экономический эксперимент
(phase_signature_economy_experiment.py) показал, что маршрутизация по всем
трём строкам даёт одну и ту же пару кандидатов {+o, -o} в 100% случаев.

Отсюда следующий шаг: хранить в bucket только одну каноническую строку
(j_min = min(j, lambda*j, lambda^2*j) mod n) и её подпись. Любая строка
орбиты вместе с подписью гиганта восстанавливает ту же пару: если
j_canon = lambda^c * j_head, маршрут даёт r' = a - c (mod 3) и база
lambda^r' * j_canon = lambda^a * j_head совпадает с базой полной таблицы.

Измеряется: сжатие storage (сырые строки -> канонические), корректность
восстановления пары, число точных проверок на hit (baseline 6 ветвей /
наивные 3 строки / каноническая 1 строка).
"""

import argparse
import json
import random
import time
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

from phase_signature_135bit_experiment import (
    P, N, LAMBDA, LAMBDA2,
    point_xy, signature_of_x, route, baseline_factors,
)


def orbit_rows(j):
    """Три x-представителя орбиты: (скаляр, x-координата)."""
    rows = []
    for s in (j, LAMBDA * j % N, LAMBDA2 * j % N):
        x, _ = point_xy(s)
        rows.append((s, x))
    return rows


def build_tables(heads, background, rng):
    """Сырые строки и каноническая таблица K3 -> одна строка."""
    raw = []  # (скаляр, x, подпись)
    for j in heads:
        for s, x in orbit_rows(j):
            raw.append((s, x, signature_of_x(x)))
    for j in background:
        x, _ = point_xy(j)
        raw.append((j, x, signature_of_x(x)))
    rng.shuffle(raw)

    naive = {}  # K3 -> список строк bucket (как хранила бы наивная таблица)
    for s, x, sig in raw:
        naive.setdefault(pow(x, 3, P), []).append((s, x, sig))
    canon = {}  # K3 -> (j, x, подпись): минимум скаляра в bucket
    for s, x, sig in raw:
        key = pow(x, 3, P)
        cur = canon.get(key)
        if cur is None or s < cur[0]:
            canon[key] = (s, x, sig)
    return raw, naive, canon


def run(args, rng):
    heads = [rng.randrange(1, N) for _ in range(args.orbits)]
    background = [rng.randrange(1, N) for _ in range(args.background)]
    t0 = time.perf_counter()
    raw, naive, canon = build_tables(heads, background, rng)
    build_seconds = time.perf_counter() - t0

    full_buckets = sum(
        1 for h in heads if pow(point_xy(h)[0], 3, P) in
        {pow(point_xy(s)[0], 3, P) for s in
         (LAMBDA * h % N, LAMBDA2 * h % N)})
    stats = {
        "orbits": args.orbits, "background": args.background,
        "raw_rows": len(raw), "canonical_rows": len(canon),
        "storage_compression_overall": round(len(raw) / len(canon), 3),
        "full_buckets": full_buckets,
        "full_bucket_compression": 3.0 if full_buckets else None,
        "build_seconds": round(build_seconds, 3),
        "trials": 0, "pair_contains_target": 0, "anomalies": 0,
        "routes": {"0": 0, "1": 0, "2": 0},
        "checks_baseline6": 0, "checks_naive3": 0, "checks_canonical": 0,
        "max_checks_naive3": 0, "max_checks_canonical": 0,
        "canonical_from_min_row": 0,
    }

    for _ in range(args.trials):
        h = heads[rng.randrange(len(heads))]
        a = rng.randrange(3)
        sign = rng.choice((1, -1))
        o = sign * (h, LAMBDA * h % N, LAMBDA2 * h % N)[a] % N
        if o == 0:
            continue
        gx, _ = point_xy(o)
        key = pow(gx, 3, P)
        canon_row = canon.get(key)
        if canon_row is None:
            continue
        stats["trials"] += 1
        gsig = signature_of_x(gx)

        # baseline: шесть ветвей первой попавшейся строки bucket
        checks = 0
        for f in baseline_factors(canon_row[0]):
            checks += 1
            if f == o:
                break
        stats["checks_baseline6"] += checks

        # наивная 3-строчная таблица: маршрутизация по каждой строке bucket
        bucket_rows = naive[key]
        checks = 0
        found = None
        for s, x, rsig in bucket_rows:
            r = route(gsig, rsig)
            if r is None:
                factors = baseline_factors(s)
            else:
                b = (s, LAMBDA * s % N, LAMBDA2 * s % N)[r]
                factors = [b, (-b) % N]
            for f in factors:
                checks += 1
                if f == o:
                    found = f
                    break
            if found is not None:
                break
        stats["checks_naive3"] += checks
        stats["max_checks_naive3"] = max(stats["max_checks_naive3"], checks)

        # каноническая таблица: одна строка
        r = route(gsig, canon_row[2])
        checks = 0
        found = None
        if r is None:
            stats["anomalies"] += 1
            factors = baseline_factors(canon_row[0])
        else:
            stats["routes"][str(r)] += 1
            b = (canon_row[0], LAMBDA * canon_row[0] % N,
                 LAMBDA2 * canon_row[0] % N)[r]
            factors = [b, (-b) % N]
        for f in factors:
            checks += 1
            if f == o:
                found = f
                break
        stats["checks_canonical"] += checks
        stats["max_checks_canonical"] = max(stats["max_checks_canonical"],
                                            checks)
        if found == o:
            stats["pair_contains_target"] += 1
        if canon_row[0] == min(s for s, _ in orbit_rows(h)):
            stats["canonical_from_min_row"] += 1

    t = max(1, stats["trials"])
    stats["mean_checks_baseline6"] = round(stats["checks_baseline6"] / t, 3)
    stats["mean_checks_naive3"] = round(stats["checks_naive3"] / t, 3)
    stats["mean_checks_canonical"] = round(stats["checks_canonical"] / t, 3)
    stats["checks_saved_vs_naive"] = round(
        stats["checks_naive3"] / max(1, stats["checks_canonical"]), 3)
    stats["pair_contains_pct"] = round(
        100.0 * stats["pair_contains_target"] / t, 1)
    return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orbits", type=int, default=400)
    parser.add_argument("--background", type=int, default=400)
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--json", type=Path,
                        default=Path("results/canonical_k3_table.json"))
    args = parser.parse_args()

    rng = random.Random(args.seed)
    print("Каноническая K3-таблица: одна строка на орбиту")
    stats = run(args, rng)
    print(f"  сырых строк: {stats['raw_rows']}, канонических: "
          f"{stats['canonical_rows']}, сжатие storage: "
          f"{stats['storage_compression_overall']}x "
          f"(полные bucket: 3x), build {stats['build_seconds']} с")
    print(f"  trials: {stats['trials']}, пара содержит цель: "
          f"{stats['pair_contains_target']}/{stats['trials']} "
          f"({stats['pair_contains_pct']}%), anomalies: {stats['anomalies']}")
    print(f"  маршруты 0/1/2: {stats['routes']['0']}/{stats['routes']['1']}/"
          f"{stats['routes']['2']}")
    print(f"  точных проверок на hit: baseline6 {stats['mean_checks_baseline6']}, "
          f"наивные 3 строки {stats['mean_checks_naive3']} "
          f"(max {stats['max_checks_naive3']}), каноническая "
          f"{stats['mean_checks_canonical']} (max "
          f"{stats['max_checks_canonical']})")
    print(f"  экономия проверок против наивной 3-строчной: "
          f"{stats['checks_saved_vs_naive']}x")

    payload = {"parameters": vars(args) | {"json": str(args.json)},
               "stats": stats}
    args.json.write_text(json.dumps(payload, indent=1, ensure_ascii=False))
    print(f"  json: {args.json}")


if __name__ == "__main__":
    main()
