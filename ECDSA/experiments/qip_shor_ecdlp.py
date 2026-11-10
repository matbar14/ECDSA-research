#!/usr/bin/env python3
"""Двухрегистровый алгоритм Шора для ECDLP на игрушечных кривых.

Постановка: группа <G> порядка n, цель H = x*G. Квантовая схема:
равномерная суперпозиция по (a, b) in Z_n^2, оракул
|a,b,0> -> |a,b,(a + b*x)*G>, измерение точечного регистра z0 = a + b*x,
QFT по Z_n^2, измерение (u, v). Амплитуды концентрируются на двойственной
решётке L^ = {(u, u*x mod n)}: v = u*x (mod n), при взаимно простом u
x = v * u^-1 (mod n) за один выстрел.

Два пути симуляции:
  gate   - гейт-точный на собственном statevector-движке (qip_engine):
           суперсингулярная кривая y^2 = x^3 + ax над F_127 (p+1 = 128 = 2^7),
           21 кубит (7+7+7), оракул - циклический сдвиг точечного регистра,
           QFT гейтами; работает при n = степени двойки.
  fft    - oracle-уровень: косет-состояние строится аналитически,
           QFT по произвольному Z_n через numpy FFT (convention e^{-2pii},
           распределение выборки совпадает с гейт-путём); 2D-FFT по месту
           (ось 0 столбцами, ось 1 строками) - пик памяти одна копия
           состояния, GLV-кривые с простым n: n = 2^13 - 1 GiB состояния
           (песочница 3 GB), n = 2^14 - 4 GiB (целевая машина 32 GB,
           пик ~6 GiB).

Кривые: (i) суперсингулярные на мерсенновских p (n = p+1 = 2^k);
(ii) GLV-кривые y^2 = x^3 + b над p = 1 (mod 3) с простым #E = 1 (mod 3) -
для эксперимента qip_glv_orbit_quantum.py (эндоморфизм phi, lambda).

Запуск из корня репо. Выход: results/qip_shor_ecdlp.json,
results/qip_toy_curves.json (кеш кривых).
"""
import argparse
import json
import math
import random
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qip_engine import Statevector, coset_state_posterior  # noqa: E402

P127 = 127


# ---------- арифметика кривых ----------

def legendre_table(p):
    tab = np.zeros(p, dtype=np.int8)
    tab[0] = 0
    for x in range(1, p):
        tab[x] = 1 if pow(x, (p - 1) // 2, p) == 1 else -1
    return tab


def count_points(p, a, b, legendre=None):
    """#E(F_p) для y^2 = x^3 + a*x + b (сложение точек на бесконечности)."""
    if legendre is None:
        legendre = legendre_table(p)
    xs = np.arange(p)
    rhs = (xs ** 3 + a * xs + b) % p
    return int(p + 1 + int(legendre[rhs].sum()))


def ec_add(p, a4, P, Q):
    """Сложение точек в аффинных координатах, None = O."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + a4) * pow(2 * y1, p - 2, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def ec_mul(p, a4, k, P):
    R = None
    Q = P
    while k:
        if k & 1:
            R = ec_add(p, a4, R, Q)
        Q = ec_add(p, a4, Q, Q)
        k >>= 1
    return R


def point_order(p, a4, P, group_order):
    """Точный порядок точки P | group_order."""
    order = group_order
    for f in prime_factors(group_order):
        while order % f == 0 and ec_mul(p, a4, order // f, P) is None:
            order //= f
    return order


def find_max_order_point(p, a4, b, group_order):
    """Точка максимального порядка (генератор циклической группы)."""
    best, best_order = None, 0
    for gx in range(p):
        rhs = (gx ** 3 + a4 * gx + b) % p
        if not rhs or pow(rhs, (p - 1) // 2, p) != 1:
            continue
        P = (gx, sqrt_mod(rhs, p))
        o = point_order(p, a4, P, group_order)
        if o > best_order:
            best, best_order = P, o
            if o == group_order:
                break
    return best, best_order


def find_point(p, a4, b):
    """Точка на кривой y^2 = x^3 + a4*x + b: первый x с QR-правой частью."""
    for gx in range(p):
        rhs = (gx ** 3 + a4 * gx + b) % p
        if rhs and pow(rhs, (p - 1) // 2, p) == 1:
            return (gx, sqrt_mod(rhs, p))
    return None


def prime_factors(m):
    fs, d = set(), 2
    while d * d <= m:
        while m % d == 0:
            fs.add(d)
            m //= d
        d += 1
    if m > 1:
        fs.add(m)
    return fs


def is_probable_prime(m):
    if m < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if m % sp == 0:
            return m == sp
    d, s = m - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        y = pow(a, d, m)
        if y in (1, m - 1):
            continue
        for _ in range(s - 1):
            y = y * y % m
            if y == m - 1:
                break
        else:
            return False
    return True


def sqrt_mod(a, p):
    """Тоннель-Шенкс для простого p."""
    if a == 0:
        return 0
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t = t * c % p
        r = r * b % p
    return r


# ---------- поиск кривых ----------

def find_toy_curves(max_p=18000, bits_buckets=(6, 8, 10, 12, 13, 14)):
    """(i) суперсингулярные на мерсенновских p (n = p+1 = 2^k);
    (ii) GLV-кривые j=0 над p = 1 (mod 3) с простым #E = 1 (mod 3),
    по одной на корзину битового размера.

    Суперсингулярная ветка ограничена мерсенновскими p = 3 (mod 4):
    127, 8191, 131071 - после n = 2^13 следующая ступень n = 2^17
    (состояние 256 GiB), поэтому рост масштаба идёт по GLV-ветке:
    FFT-путь через косет-состояние не зависит от суперсингулярности.
    """
    curves = []

    for p in (127, 8191, 131071):
        if p > max_p:
            continue
        leg = legendre_table(p)
        for a in range(1, 6):
            n = count_points(p, a, 0, leg)
            if n == p + 1:
                curves.append({"p": p, "a": a, "b": 0, "order": n,
                               "kind": "supersingular",
                               "family": "y^2=x^3+ax, p=3 (mod 4), #E=p+1"})
                break

    xs = None
    for target_bits in bits_buckets:
        found = None
        for p in range(7, max_p + 1, 2):
            if not is_probable_prime(p) or p % 3 != 1:
                continue
            if xs is None or xs.size != p:
                xs = np.arange(p)
                leg = legendre_table(p)
            x3 = xs ** 3 % p
            for b in range(1, 64):
                rhs = (x3 + b) % p
                n = p + 1 + int(leg[rhs].sum())
                if not is_probable_prime(n) or n % 3 != 1:
                    continue
                if n.bit_length() != target_bits:
                    continue
                found = {"p": p, "a": 0, "b": b, "order": n,
                         "kind": "glv",
                         "family": "y^2=x^3+b, j=0, p=1 (mod 3), "
                                   "#E простое = 1 (mod 3)"}
                break
            if found:
                break
        if found:
            curves.append(found)
            print(f"[search] GLV-кривая на {target_bits} бит: "
                  f"p={found['p']} b={found['b']} n={found['order']}")
        else:
            print(f"[search] кривая на {target_bits} бит не найдена "
                  f"до p={max_p}")
    return curves


def glv_lambda(p, n):
    """lambda: lambda^2 + lambda + 1 = 0 (mod n), lambda != 1;
    верифицируется на точках вызывающим кодом."""
    r = sqrt_mod(n - 3, n)  # sqrt(-3) mod n
    inv2 = (n + 1) // 2
    for lam in (((n - 1 + r) * inv2) % n, ((n - 1 - r) * inv2) % n):
        if (lam * lam + lam + 1) % n == 0 and lam != 1:
            return lam
    return None


# ---------- гейт-точный путь (n = 2^w) ----------

def shor_gate_exact(curve, x, rng, shots):
    p, a4, n = curve["p"], curve["a"], curve["order"]
    w = n.bit_length() - 1
    assert n == 2 ** w
    # оракул: |a,b,z> -> |a,b,(z + a + b*x) mod n> - циклический сдвиг
    ab = np.arange(n * n).reshape(n, n)
    shift = (ab // n + (ab % n) * x) % n  # f(a,b) для каждой пары
    z_idx = (np.arange(n)[None, None, :]
             - shift[:, :, None]) % n  # gather: new[(a,b,z)] = old[(a,b,z-f)]
    # плоский индекс (a,b,z): a*n^2 + b*n + z
    flat_idx = ((np.arange(n)[:, None, None] * n * n)
                + (np.arange(n)[None, :, None] * n)
                + z_idx).reshape(-1)

    # Конвенция (важно): плоский индекс i = a*n^2 + b*n + z, т.е. a -
    # старшие биты (3w-1..2w), b - средние (2w-1..w), z - младшие
    # (w-1..0); кубит q соответствует биту (3w-1-q). Поэтому a занимает
    # кубиты 0..w-1 (MSB=q0), b - w..2w-1, z - 2w..3w-1; списки для
    # QFT/измерения идут MSB-первым. Ранняя версия задавала списки в
    # обратной конвенции: QFT и измерение получали поля a/b переставленными
    # и с реверсом значимости битов, выборки ложились мимо двойственной
    # решётки (0/16), x не восстанавливался. Регрессионный тест конвенции
    # оракула - qip_test_engine.py, тест 7.
    q_a = list(range(0, w))          # a: q0(MSB)..q(w-1)
    q_b = list(range(w, 2 * w))      # b: qw(MSB)..q(2w-1)
    q_z = list(range(2 * w, 3 * w))  # z: q2w(MSB)..q(3w-1)
    results = []
    for _ in range(shots):
        sv = Statevector(3 * w)
        for q in q_a + q_b:
            sv.h(q)
        # оракул одним gather: новый вектор по индексам
        sv.psi = sv.psi[flat_idx]
        z0, _ = sv.measure(rng, qubits=q_z)
        sv.qft_pow2(q_a)
        sv.qft_pow2(q_b)
        u, _ = sv.measure(rng, qubits=q_a)
        v, _ = sv.measure(rng, qubits=q_b)
        results.append((u, v))
    return results


def recover_x_2adic(n, samples):
    """CRT-восстановление x по выборкам (u, v) с v = u*x (mod n)."""
    if not samples:
        return None
    x_mod, modulus = 0, 1
    for u, v in samples:
        j = 0
        uu, vv = u, v
        while uu % 2 == 0 and vv % 2 == 0 and j < n.bit_length() - 1:
            uu //= 2
            vv //= 2
            j += 1
        if uu % 2 == 1:
            # x = v/u (mod n/2^j)
            m = n >> j
            cand = vv * pow(uu, -1, m) % m
            # согласование с уже накопленным
            if modulus % m == 0:
                if cand != x_mod % m:
                    return None
            else:
                # CRT по степеням двойки - берём самый точный снимок
                if m > modulus:
                    x_mod, modulus = cand, m
    return x_mod if modulus == n else None


# ---------- FFT-путь (произвольное n) ----------

def shor_fft(curve, x, rng, shots):
    """Косет-состояние + QFT по Z_n^2 через numpy FFT.

    2D-FFT считается по месту: ось 0 столбцами, затем ось 1 строками
    (эквивалент np.fft.fft2 с norm="ortho" - та же нормировка 1/n),
    но пик памяти - одна копия состояния плюс один столбец/строка,
    а не две полные копии: n = 2^13 даёт 1 GiB состояния и проходит
    в песочнице, np.fft.fft2 требовал бы ~3 GiB на тот же n.
    """
    n = curve["order"]
    results, t_build, t_qft = [], 0.0, 0.0
    for _ in range(shots):
        z0 = rng.randrange(n)  # измерение точечного регистра: z = a+bx
        t0 = time.perf_counter()
        psi = coset_state_posterior(n, x, z0)
        t_build += time.perf_counter() - t0
        t0 = time.perf_counter()
        for j in range(n):
            psi[:, j] = np.fft.fft(psi[:, j], norm="ortho")
        for i in range(n):
            psi[i] = np.fft.fft(psi[i], norm="ortho")
        t_qft += time.perf_counter() - t0
        prob = np.abs(psi) ** 2
        del psi
        prob /= prob.sum()
        idx = np.searchsorted(np.cumsum(prob), rng.random())
        u, v = divmod(int(idx), n)
        results.append((u, v))
        del prob
    return results, t_build, t_qft


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shots", type=int, default=16)
    ap.add_argument("--min-n", type=int, default=0,
                    help="пропустить кривые с порядком ниже этого")
    ap.add_argument("--max-n", type=int, default=4096,
                    help="максимальный порядок кривой для FFT-пути")
    ap.add_argument("--seed", type=int, default=20261103)
    ap.add_argument("--out", default="results/qip_shor_ecdlp.json",
                    help="путь выходного JSON")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    cache = Path("results/qip_toy_curves.json")
    if cache.exists():
        curves = json.loads(cache.read_text())
    else:
        curves = find_toy_curves()
        cache.write_text(json.dumps(curves, indent=2) + "\n")

    report = {"shots per run": args.shots, "runs": [], "equivalence": None}
    print(f"кривых в кеше: {len(curves)}")

    # --- верификация GLV-эндоморфизма на кривых kind=glv ---
    for c in curves:
        p, n = c["p"], c["order"]
        if c["kind"] == "supersingular":
            G, mord = find_max_order_point(p, c["a"], c["b"], n)
            c["generator"] = list(G) if G else None
            c["max order"] = mord
            print(f"[supersingular] p={p} a={c['a']}: #E={n}, "
                  f"порядок генератора={mord}")
            continue
        lam = glv_lambda(p, n)
        beta = 1
        for g in range(2, 64):
            beta = pow(g, (p - 1) // 3, p)
            if beta != 1:
                break
        G = find_point(p, c["a"], c["b"])
        c["lambda"] = None
        if G and ec_mul(p, c["a"], n, G) is None:
            phi_P = (beta * G[0] % p, G[1])
            if ec_mul(p, c["a"], lam, G) == phi_P:
                c["lambda"] = lam
            elif ec_mul(p, c["a"], lam * lam % n, G) == phi_P:
                c["lambda"] = lam * lam % n
        c["generator"] = list(G) if G else None
        st = ("lambda ok" if c.get("lambda") else "lambda fail")
        print(f"[glv] p={p} b={c['b']} n={n} ({n.bit_length()} бит) - {st}")

    # --- гейт-точный путь: p=127, n=128 ---
    ss = next(c for c in curves if c["p"] == P127)
    n = ss["order"]
    w = n.bit_length() - 1
    x = rng.randrange(1, n)
    print(f"[gate] p={P127} n={n}=2^{w}, цель x={x}")
    t0 = time.perf_counter()
    samples = shor_gate_exact(ss, x, rng, args.shots)
    gate_time = time.perf_counter() - t0
    rec = recover_x_2adic(n, samples)
    partial = sum(1 for u, v in samples
                  if u and v == u * x % n)
    print(f"[gate] {args.shots} выстрелов за {gate_time:.2f} с; на решётке "
          f"{partial}/{args.shots}; восстановлен x={rec} "
          f"({'PASS' if rec == x else 'FAIL'})")
    report["runs"].append({
        "path": "gate", "p": P127, "n": n, "x": x, "shots": args.shots,
        "on lattice": partial, "recovered": rec, "seconds": round(gate_time, 2),
        "qubits": 3 * w,
    })

    # --- эквивалентность гейт- и FFT-путей при том же n ---
    samples_fft, _, _ = shor_fft(ss, x, rng, args.shots)
    gate_ok = all(v == u * x % n for u, v in samples)
    fft_ok = all(v == u * x % n for u, v in samples_fft)
    report["equivalence"] = {
        "n": n,
        "gate on lattice": f"{partial}/{args.shots}",
        "fft on lattice":
            f"{sum(1 for u, v in samples_fft if v == u * x % n)}/{args.shots}",
        "both on dual lattice": bool(gate_ok and fft_ok),
    }
    print(f"[equivalence] gate={gate_ok}, fft={fft_ok}")

    # --- FFT-путь по всем кривым ---
    for c in curves:
        n = c["order"]
        if n > args.max_n or n < args.min_n:
            continue
        x = rng.randrange(1, n)
        t0 = time.perf_counter()
        samples, t_build, t_qft = shor_fft(c, x, rng, args.shots)
        wall = time.perf_counter() - t0
        ok = sum(1 for u, v in samples
                 if u and v * pow(u, -1, n) % n == x) if is_probable_prime(n) else None
        u0 = sum(1 for u, v in samples if u == 0)
        entry = {
            "path": "fft", "p": c["p"], "n": n, "bits": n.bit_length(),
            "kind": c["kind"], "x": x, "shots": args.shots,
            "u zero": u0, "seconds": round(wall, 3),
            "build s": round(t_build, 3), "qft s": round(t_qft, 3),
            "statevector MB": round(16 * n * n / 2 ** 20, 1),
        }
        if ok is not None:
            entry["recovered"] = f"{ok}/{args.shots}"
            print(f"[fft] p={c['p']} n={n} ({n.bit_length()} бит, "
                  f"{c['kind']}): восстановлен {ok}/{args.shots} за "
                  f"{wall:.2f} с ({entry['statevector MB']} MB)")
        else:
            rec = recover_x_2adic(n, samples)
            entry["recovered"] = rec == x
            print(f"[fft] p={c['p']} n={n} (2^{n.bit_length() - 1}): "
                  f"CRT-восстановление {'PASS' if rec == x else 'FAIL'} "
                  f"за {wall:.2f} с")
        report["runs"].append(entry)

    out = Path(args.out)
    out.write_text(json.dumps({"curves": curves, **report}, indent=2,
                              ensure_ascii=False) + "\n")
    print(f"[saved] {out}")


if __name__ == "__main__":
    main()
