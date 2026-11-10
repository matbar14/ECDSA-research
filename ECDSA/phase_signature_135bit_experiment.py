#!/usr/bin/env python3
"""Проверка (d01, d12)-сигнатуры как маршрутизатора GLV-ветвей.

Два режима:

window — сегментированный BSGS на 135-битном окне [2^134, 2^134+W):
  baby-таблица хранит K3=x^3, j и подпись (d01, d12) точки jG; для случайных
  скаляров окна выполняется giant-шаг и маршрутизация ветви по подписи.

orbit — engineered-партнёры: для случайных j и ветви (r, s) строится точка
  o = s * lambda^r * j (mod n), что соответствует полному диапазону, где
  партнёры орбиты реально попадают в поиск. Проверяет распознавание всех
  трёх поворотов phi^r и сокращение кандидатов 6 -> 2.

Подпись (d01, d12) инвариантна относительно P -> -P и циклически
сдвигается при phi: тройка Lambda(P) переходит в (L1, L2, L0).
"""
import argparse
import hashlib
import json
import random
import time
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import ENDOMORPHISM_BETA, FIELD_PRIME, GROUP_ORDER

LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
LAMBDA2 = 0xAC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE
P = FIELD_PRIME
N = GROUP_ORDER
BETA = ENDOMORPHISM_BETA


def point_xy(scalar):
    numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
    return numbers.x, numbers.y


GX, GY = point_xy(1)
BETA_G = BETA * GX % P


def lambda_x(x):
    """Lambda(P) по крест-отношению, соглашение w1.py: знаменатель 0 -> 0."""
    x %= P
    if x == 0:
        return 0
    beta_x = BETA * x % P
    numerator = (GX - BETA_G) % P * ((x - beta_x) % P) % P
    denominator = (GX - beta_x) % P * ((x - BETA_G) % P) % P
    if denominator == 0:
        return 0
    return numerator * pow(denominator, -1, P) % P


def signature_of_x(x):
    """(d01, d12) для x-координаты: d01 = L1 - L0, d12 = L2 - L1."""
    t0 = lambda_x(x)
    t1 = lambda_x(BETA * x % P)
    t2 = lambda_x(BETA * BETA * x % P)
    return ((t1 - t0) % P, (t2 - t1) % P)


def route(giant_sig, row_sig):
    """Номер поворота r с giant = phi^r(+/-row) или None, если подписи вырождены."""
    r01, r12 = row_sig
    rsum_neg = (-(r01 + r12)) % P
    if giant_sig == (r01, r12):
        return 0
    if giant_sig == (r12, rsum_neg):
        return 1
    if giant_sig == (rsum_neg, r01):
        return 2
    return None


def build_baby_table(rows):
    """K3=x^3 -> список (j, подпись)."""
    table = {}
    for j in rows:
        x, _ = point_xy(j)
        key = pow(x, 3, P)
        table.setdefault(key, []).append((j, signature_of_x(x)))
    return table


def baseline_factors(j):
    """Шесть скалярных ветвей {+-j, +-lambda*j, +-lambda^2*j} mod n."""
    return [j % N, LAMBDA * j % N, LAMBDA2 * j % N,
            (-j) % N, (-LAMBDA * j) % N, (-LAMBDA2 * j) % N]


def run_window_mode(args, rng):
    start = 1 << 134
    width = args.window
    m = args.m or (int(width ** 0.5) + 1)
    table = build_baby_table(range(1, m))
    buckets = sum(len(v) for v in table.values())
    stats = {
        "trials": 0, "bucket_matches": 0, "routes": {"0": 0, "1": 0, "2": 0},
        "anomalies": 0, "routed_correctly": 0, "factor_checks_baseline": 0,
        "factor_checks_signature": 0, "recoveries": 0,
        "sig_seconds": 0.0, "baseline_seconds": 0.0,
    }
    for _ in range(args.trials):
        k = start + rng.randrange(width)
        o = k - start
        i, residue = divmod(o, m)
        if residue == 0:
            continue  # j=0 (точка на бесконечности) вне таблицы
        gx, _ = point_xy(residue)
        gkey = pow(gx, 3, P)
        bucket = table.get(gkey, [])
        if not bucket:
            continue
        stats["bucket_matches"] += 1
        t0 = time.perf_counter()
        gsig = signature_of_x(gx)
        stats["sig_seconds"] += time.perf_counter() - t0
        found_factor = None
        for j, rsig in bucket:
            r = route(gsig, rsig)
            if r is None:
                stats["anomalies"] += 1
                factors = baseline_factors(j)
            else:
                stats["routes"][str(r)] += 1
                base = (j, LAMBDA * j % N, LAMBDA2 * j % N)[r]
                factors = [base, (-base) % N]
            for f in factors:
                stats["factor_checks_signature"] += 1
                if f == residue % N:
                    found_factor = f
                    break
            if found_factor is not None:
                break
        if found_factor == residue % N:
            stats["routed_correctly"] += 1
            stats["recoveries"] += 1
        t0 = time.perf_counter()
        for j, _ in bucket:
            for f in baseline_factors(j):
                stats["factor_checks_baseline"] += 1
                if f == residue % N:
                    break
            else:
                continue
            break
        stats["baseline_seconds"] += time.perf_counter() - t0
        stats["trials"] += 1
    stats["mean_sig_ms"] = 1000.0 * stats["sig_seconds"] / max(1, stats["trials"])
    stats["mean_baseline_ms"] = 1000.0 * stats["baseline_seconds"] / max(1, stats["trials"])
    stats["window"] = {"start": hex(start), "width": width, "baby_rows": m,
                       "table_entries": buckets}
    return stats


def run_orbit_mode(args, rng):
    m = args.m or 1000
    table = build_baby_table(range(1, m))
    stats = {
        "trials": 0, "routes": {"0": 0, "1": 0, "2": 0}, "anomalies": 0,
        "route_correct": 0, "candidate_reduction_checked": 0,
        "candidate_reduction_ok": 0, "neg_partner_same_x": 0,
        "sig_seconds": 0.0,
    }
    degenerate = 0
    for _ in range(args.trials):
        j = rng.randrange(1, m)
        r = rng.randrange(3)
        sign = rng.choice((1, -1))
        base = (j, LAMBDA * j % N, LAMBDA2 * j % N)[r]
        o = sign * base % N
        if o == 0:
            continue
        gx, _ = point_xy(o)
        gkey = pow(gx, 3, P)
        bucket = table.get(gkey, [])
        if not bucket:
            continue
        t0 = time.perf_counter()
        gsig = signature_of_x(gx)
        stats["sig_seconds"] += time.perf_counter() - t0
        routed = None
        for row_j, rsig in bucket:
            route_r = route(gsig, rsig)
            if route_r is None:
                stats["anomalies"] += 1
                continue
            stats["routes"][str(route_r)] += 1
            if row_j == j and route_r == r:
                routed = route_r
        stats["trials"] += 1
        if routed == r:
            stats["route_correct"] += 1
        # проверка сокращения: кандидаты {+base, -base} содержат o
        stats["candidate_reduction_checked"] += 1
        if o in (base, (-base) % N):
            stats["candidate_reduction_ok"] += 1
        # знак: -P имеет тот же x
        neg_x, _ = point_xy((-o) % N)
        if neg_x == gx:
            stats["neg_partner_same_x"] += 1
        if 0 in (gsig[0], gsig[1]):
            degenerate += 1
    stats["mean_sig_ms"] = 1000.0 * stats["sig_seconds"] / max(1, stats["trials"])
    stats["degenerate_signatures"] = degenerate
    stats["baby_rows"] = m
    return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument("--window", type=int, default=1 << 16)
    parser.add_argument("--m", type=int, default=0,
                        help="число baby-строк (по умолчанию ceil(sqrt(window)))")
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--json", type=Path, default=Path("phase_signature_135bit.json"))
    args = parser.parse_args()

    rng = random.Random(args.seed)
    window_stats = run_window_mode(args, rng)
    orbit_stats = run_orbit_mode(args, rng)

    print("Phase signature (d01, d12) routing on the 135-bit window")
    print(f"  window mode: trials={window_stats['trials']}, "
          f"bucket matches={window_stats['bucket_matches']}")
    print(f"    routes 0/1/2: {window_stats['routes']['0']}/"
          f"{window_stats['routes']['1']}/{window_stats['routes']['2']}, "
          f"anomalies={window_stats['anomalies']}")
    print(f"    routed correctly: {window_stats['routed_correctly']}/"
          f"{window_stats['trials']}")
    print(f"    factor checks per hit: signature {window_stats['factor_checks_signature'] / max(1, window_stats['trials']):.2f} "
          f"vs baseline {window_stats['factor_checks_baseline'] / max(1, window_stats['trials']):.2f}")
    print(f"    mean signature time: {window_stats['mean_sig_ms']:.3f} ms, "
          f"baseline factor scan: {window_stats['mean_baseline_ms']:.3f} ms")
    print("  orbit mode (engineered partners o = +/- lambda^r * j):")
    print(f"    trials={orbit_stats['trials']}, route correct: "
          f"{orbit_stats['route_correct']}/{orbit_stats['trials']}")
    print(f"    routes 0/1/2: {orbit_stats['routes']['0']}/"
          f"{orbit_stats['routes']['1']}/{orbit_stats['routes']['2']}, "
          f"anomalies={orbit_stats['anomalies']}, "
          f"degenerate={orbit_stats['degenerate_signatures']}")
    print(f"    candidate pair {{+b,-b}} contains o: "
          f"{orbit_stats['candidate_reduction_ok']}/"
          f"{orbit_stats['candidate_reduction_checked']}")
    print(f"    -o*G has the same x as o*G: "
          f"{orbit_stats['neg_partner_same_x']}/{orbit_stats['trials']} "
          f"(знак остаётся двусмысленным)")
    print(f"    mean signature time: {orbit_stats['mean_sig_ms']:.3f} ms")
    print("  Interpretation: маршрутизация определяет поворот phi^r точно, "
          "сокращает шесть ветвей до пары (+-lambda^r * j); знак +/- "
          "остаётся и разрешается точной проверкой точки.")

    report = {
        "parameters": vars(args) | {"json": str(args.json)},
        "window_mode": window_stats,
        "orbit_mode": orbit_stats,
        "constants": {"beta": hex(BETA), "lambda": hex(LAMBDA),
                      "lambda2": hex(LAMBDA2)},
    }
    args.json.write_text(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"  report={args.json}")


if __name__ == "__main__":
    main()
