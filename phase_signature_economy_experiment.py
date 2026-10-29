#!/usr/bin/env python3
"""Экономика (d01,d12)-маршрутизации в сценарии полного диапазона.

PROGRESS_MEMO.md («Фазовая подпись», практический вывод) предсказывает, что
выгода маршрутизации появляется там, где в K3-bucket реально попадают
несколько партнёров орбиты, то есть при полном диапазоне или таблицах с
настоящими K3-коллизиями. Этот эксперимент моделирует именно такой случай:
таблица представителей разбросана по всему [1, n), а цель выбирается среди
орбитальных партнёров строк.

Два режима:

  single  — таблица из m случайных представителей, bucket содержит одну строку;
            измеряются точные проверки ветвей до/после маршрутизации.
  bucket  — для части строк в таблице лежат все три x-представителя орбиты
            (j, lambda*j, lambda^2*j), bucket содержит 3 строки; проверяется
            согласованность маршрутов по строкам и схлопывание кандидатов
            в одну пару {+o, -o}.

Точная проверка кандидата моделируется сравнением скаляров (как в
phase_signature_135bit_experiment.py): в движке она стоит полное скалярное
умножение, поэтому счёт проверок — главный показатель экономики.
"""

import argparse
import json
import random
import time
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

from phase_signature_135bit_experiment import (
    P, N, LAMBDA, LAMBDA2,
    point_xy, lambda_x, signature_of_x, route, baseline_factors,
)


def k3_of_scalar(scalar):
    x, _ = point_xy(scalar)
    return pow(x, 3, P), x


def build_table(rows, sig_seconds_acc=None):
    """K3 = x^3 -> список (j, подпись)."""
    table = {}
    for j in rows:
        _, x = k3_of_scalar(j)
        t0 = time.perf_counter()
        sig = signature_of_x(x)
        if sig_seconds_acc is not None:
            sig_seconds_acc[0] += time.perf_counter() - t0
        table.setdefault(pow(x, 3, P), []).append((j, sig))
    return table


def run_single_representative(args, rng):
    rows = [rng.randrange(1, N) for _ in range(args.rows)]
    sig_clock = [0.0]
    table = build_table(rows, sig_clock)

    stats = {
        "mode": "single",
        "trials": 0, "recoveries": 0, "anomalies": 0,
        "routes": {"0": 0, "1": 0, "2": 0},
        "checks_baseline": 0, "checks_signature": 0,
        "max_checks_baseline": 0, "max_checks_signature": 0,
        "sig_seconds": 0.0,
    }
    for _ in range(args.trials):
        j = rows[rng.randrange(len(rows))]
        r = rng.randrange(3)
        sign = rng.choice((1, -1))
        base = (j, LAMBDA * j % N, LAMBDA2 * j % N)[r]
        o = sign * base % N
        if o == 0:
            continue
        _, gx = k3_of_scalar(o)
        bucket = table.get(pow(gx, 3, P), [])
        if not bucket:
            continue
        stats["trials"] += 1

        t0 = time.perf_counter()
        gsig = signature_of_x(gx)
        stats["sig_seconds"] += time.perf_counter() - t0

        # маршрут: пара {+lambda^r*j, -lambda^r*j} на строку
        found = None
        checks_sig = 0
        for row_j, rsig in bucket:
            rr = route(gsig, rsig)
            if rr is None:
                stats["anomalies"] += 1
                factors = baseline_factors(row_j)
            else:
                stats["routes"][str(rr)] += 1
                b = (row_j, LAMBDA * row_j % N, LAMBDA2 * row_j % N)[rr]
                factors = [b, (-b) % N]
            for f in factors:
                checks_sig += 1
                if f == o:
                    found = f
                    break
            if found is not None:
                break
        stats["checks_signature"] += checks_sig
        stats["max_checks_signature"] = max(
            stats["max_checks_signature"], checks_sig)
        if found == o:
            stats["recoveries"] += 1

        # базовый путь: шесть ветвей на строку без маршрутизации
        checks = 0
        found_b = None
        for row_j, _ in bucket:
            for f in baseline_factors(row_j):
                checks += 1
                if f == o:
                    found_b = f
                    break
            else:
                continue
            break
        stats["checks_baseline"] += checks
        stats["max_checks_baseline"] = max(stats["max_checks_baseline"], checks)

    t = max(1, stats["trials"])
    stats["mean_checks_baseline"] = round(stats["checks_baseline"] / t, 3)
    stats["mean_checks_signature"] = round(stats["checks_signature"] / t, 3)
    stats["checks_saved_ratio"] = round(
        stats["checks_baseline"] / max(1, stats["checks_signature"]), 3)
    stats["mean_sig_ms"] = round(1000.0 * stats["sig_seconds"] / t, 4)
    stats["rows"] = args.rows
    return stats


def run_full_orbit_bucket(args, rng):
    rows = []
    orbit_heads = []
    for _ in range(args.orbits):
        j = rng.randrange(1, N)
        orbit_heads.append(j)
        rows.extend((j, LAMBDA * j % N, LAMBDA2 * j % N))
    rows.extend(rng.randrange(1, N) for _ in range(args.background))
    sig_clock = [0.0]
    table = build_table(rows, sig_clock)

    stats = {
        "mode": "bucket",
        "trials": 0, "consistent": 0, "anomalies": 0,
        "routes": {"0": 0, "1": 0, "2": 0},
        "single_pair_buckets": 0,
        "checks_baseline": 0, "checks_signature": 0,
        "max_checks_baseline": 0, "max_checks_signature": 0,
        "pair_contains_target": 0,
        "sig_seconds": 0.0, "build_sig_seconds": sig_clock[0],
        "table_rows": len(rows),
    }
    for _ in range(args.trials):
        j = orbit_heads[rng.randrange(len(orbit_heads))]
        a = rng.randrange(3)
        sign = rng.choice((1, -1))
        o = sign * (j, LAMBDA * j % N, LAMBDA2 * j % N)[a] % N
        if o == 0:
            continue
        _, gx = k3_of_scalar(o)
        bucket = table.get(pow(gx, 3, P), [])
        if len(bucket) < 3:
            continue  # интересует только полный орбитальный bucket
        stats["trials"] += 1

        t0 = time.perf_counter()
        gsig = signature_of_x(gx)
        stats["sig_seconds"] += time.perf_counter() - t0

        # маршруты по всем трём строкам одной орбиты
        rs = [route(gsig, rsig) for _, rsig in bucket]
        if any(r is None for r in rs):
            stats["anomalies"] += 1

        # кандидаты после маршрутизации: одна пара на строку, дедуп по множеству
        found = None
        seen_pairs = set()
        checks_sig = 0
        for row_j, rsig in bucket:
            rr = route(gsig, rsig)
            if rr is None:
                factors = baseline_factors(row_j)
            else:
                b = (row_j, LAMBDA * row_j % N, LAMBDA2 * row_j % N)[rr]
                factors = [b, (-b) % N]
                seen_pairs.add((min(factors), max(factors)))
            for f in factors:
                checks_sig += 1
                if f == o:
                    found = f
                    break
            if found is not None:
                break
        stats["checks_signature"] += checks_sig
        stats["max_checks_signature"] = max(
            stats["max_checks_signature"], checks_sig)
        if len(seen_pairs) == 1:
            stats["single_pair_buckets"] += 1
        if found == o:
            stats["pair_contains_target"] += 1

        # согласованность маршрутов: по каждой строке база пары одна и та же
        consensus_ok = len(seen_pairs) == 1 and not any(r is None for r in rs)
        if consensus_ok:
            stats["consistent"] += 1
            stats["routes"][str(rs[0])] += 1

        # базовый путь: строки bucket перебираются по порядку
        checks = 0
        for row_j, _ in bucket:
            stop = False
            for f in baseline_factors(row_j):
                checks += 1
                if f == o:
                    stop = True
                    break
            if stop:
                break
        stats["checks_baseline"] += checks
        stats["max_checks_baseline"] = max(stats["max_checks_baseline"], checks)

    t = max(1, stats["trials"])
    stats["mean_checks_baseline"] = round(stats["checks_baseline"] / t, 3)
    stats["mean_checks_signature"] = round(stats["checks_signature"] / t, 3)
    stats["checks_saved_ratio"] = round(
        stats["checks_baseline"] / max(1, stats["checks_signature"]), 3)
    stats["consistent_pct"] = round(100.0 * stats["consistent"] / t, 1)
    stats["single_pair_pct"] = round(100.0 * stats["single_pair_buckets"] / t, 1)
    stats["mean_sig_ms"] = round(1000.0 * stats["sig_seconds"] / t, 4)
    stats["orbits"] = args.orbits
    stats["background"] = args.background
    return stats


def field_cost_model(single, bucket):
    """Грубая модель полевых операций (modmul-эквиваленты).

    Точная проверка кандидата в движке = полное скалярное умножение:
    ~256 dbl + ~128 add, dbl ~4M, add ~12M -> ~2600 modmul + 1 inv (~390M).
    Подпись гиганта = 3 lambda = 3 modinv (~1170 modmul) + линейные умножения.
    """
    check_cost = 3000      # modmul на одну точную проверку
    inv_cost = 390         # modmul на одну модульную инверсию
    sig_cost = 3 * inv_cost
    saved_single = single["mean_checks_baseline"] - single["mean_checks_signature"]
    saved_bucket = bucket["mean_checks_baseline"] - bucket["mean_checks_signature"]
    return {
        "one_exact_check_modmul": check_cost,
        "one_modinv_modmul": inv_cost,
        "signature_cost_modmul": sig_cost,
        "saved_checks_per_hit_single": round(saved_single, 3),
        "saved_checks_per_hit_bucket": round(saved_bucket, 3),
        "net_win_modmul_per_hit_single": round(
            saved_single * check_cost - sig_cost),
        "net_win_modmul_per_hit_bucket": round(
            saved_bucket * check_cost - sig_cost),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=int, default=2000,
                        help="представителей в таблице режима single")
    parser.add_argument("--orbits", type=int, default=400,
                        help="полных орбит (3 строки) в режиме bucket")
    parser.add_argument("--background", type=int, default=400,
                        help="фоновых случайных строк в режиме bucket")
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--json", type=Path,
                        default=Path("phase_signature_economy.json"))
    args = parser.parse_args()

    rng = random.Random(args.seed)
    print("Экономика (d01,d12)-маршрутизации: сценарий полного диапазона")
    print(f"  seed={args.seed}")

    single = run_single_representative(args, rng)
    print("  single-representative таблица (m случайных из [1, n)):")
    print(f"    trials: {single['trials']}, recoveries: {single['recoveries']}, "
          f"anomalies: {single['anomalies']}")
    print(f"    routes 0/1/2: {single['routes']['0']}/{single['routes']['1']}/"
          f"{single['routes']['2']}")
    print(f"    точных проверок на hit: baseline {single['mean_checks_baseline']} "
          f"(max {single['max_checks_baseline']}), "
          f"подпись {single['mean_checks_signature']} "
          f"(max {single['max_checks_signature']})")
    print(f"    экономия проверок: {single['checks_saved_ratio']}x, "
          f"цена подписи {single['mean_sig_ms']} мс (Python)")

    bucket = run_full_orbit_bucket(args, rng)
    print("  bucket из трёх x-представителей одной орбиты:")
    print(f"    trials: {bucket['trials']}, согласованных маршрутов: "
          f"{bucket['consistent']} ({bucket['consistent_pct']}%), "
          f"anomalies: {bucket['anomalies']}")
    print(f"    bucket схлопнулся в одну пару кандидатов: "
          f"{bucket['single_pair_buckets']} ({bucket['single_pair_pct']}%)")
    print(f"    пара содержит цель: {bucket['pair_contains_target']}")
    print(f"    точных проверок на hit: baseline {bucket['mean_checks_baseline']} "
          f"(max {bucket['max_checks_baseline']}), "
          f"подпись {bucket['mean_checks_signature']} "
          f"(max {bucket['max_checks_signature']})")
    print(f"    экономия проверок: {bucket['checks_saved_ratio']}x")

    model = field_cost_model(single, bucket)
    print("  модель полевых операций (modmul-эквиваленты):")
    print(f"    одна точная проверка ~ {model['one_exact_check_modmul']} modmul, "
          f"инверсия ~ {model['one_modinv_modmul']}, подпись ~ "
          f"{model['signature_cost_modmul']}")
    print(f"    чистый выигрыш на hit: single "
          f"{model['net_win_modmul_per_hit_single']} modmul, bucket "
          f"{model['net_win_modmul_per_hit_bucket']} modmul")

    payload = {
        "parameters": {
            "rows": args.rows, "orbits": args.orbits,
            "background": args.background, "trials": args.trials,
            "seed": args.seed, "json": str(args.json),
        },
        "single_mode": single,
        "bucket_mode": bucket,
        "field_cost_model": model,
    }
    args.json.write_text(json.dumps(payload, indent=1, ensure_ascii=False))
    print(f"  json: {args.json}")


if __name__ == "__main__":
    main()
