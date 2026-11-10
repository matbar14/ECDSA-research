#!/usr/bin/env python3
"""Проверка точной структурной формулы Тоффоли-счёта предиката на замере.

Формула (из строения схемы, всё считается покомпонентно):
    T(n) = 8*n*(n+m) + 64*(2m-1)*(n+m) + 20*n + 45,   m = ceil(log2 n)

    8*n*(n+m)       - фаза умножений: 2 умножения * n сложений * 2(n+m)
                      Тоффоли на Cuccaro-сложение, x2 на очистку реверсом;
    64*(2m-1)*(n+m) - каскад бинарного сокращения: m cond_sub * 16(n+m)
                      (компаратор 5 + реверс 5 + аддер 2 + реверс 2 + мульти-
                      плексор 2) + реверс (m-1) промежуточных шагов, x2
                      умножения x2 очистка;
    20*n            - канонизация: 4 компаратора [A<K] * 5n Тоффоли/бит
                      (2 прямых + 2 в очистке);
    45              - интервальная логика (2 компаратора на 2 битах, прямые
                      и в очистке: 40) + AND-дерево (5).
Развёрнутый вид: T(n) = 8n^2 + 136*m*n - 44n + 128*m^2 - 64m + 45.
"""
import json
import math
from pathlib import Path

RES = Path("/home/z/my-project/results/qip_predicate_kappa.json")


def m_of(n):
    return max(1, (n - 1).bit_length())


def structural(n):
    m = m_of(n)
    return 8 * n * (n + m) + 64 * (2 * m - 1) * (n + m) + 20 * n + 45


def main():
    d = json.loads(RES.read_text())
    rows = d["scaling"]["rows"]
    print("Проверка структурной формулы на замере:")
    print("n     measured   formula     diff")
    bad = 0
    for r in rows:
        n = r["n"]
        f = structural(n)
        ok = f == r["toffoli"]
        bad += not ok
        print(f"{n:4d} {r['toffoli']:9d} {f:11d} {f - r['toffoli']:6d}  "
              f"{'OK' if ok else 'MISMATCH'}")
    print("mismatches:", bad)

    n = 256
    m = m_of(n)
    total = structural(n)
    p1 = 8 * n * (n + m)
    red = 64 * (2 * m - 1) * (n + m)
    can = 20 * n
    print(f"\nn=256, разложение (итого {total}):")
    print(f"  умножения 8n(n+m):            {p1:8d}  ({p1 / n / n:5.2f} n^2, "
          f"{100 * p1 / total:.1f}%)")
    print(f"  сокращение 64(2m-1)(n+m):     {red:8d}  ({red / n / n:5.2f} n^2, "
          f"{100 * red / total:.1f}%)")
    print(f"  канонизация 20n:              {can:8d}  ({can / n / n:5.2f} n^2, "
          f"{100 * can / total:.1f}%)")
    print(f"  интервал+AND:                     45")
    nl = 136 * m * n
    print(f"член n*log2(n): 136*m*n = {nl} = {100 * nl / total:.1f}% счёта")

    a, b, c = 8.4394, 1007.832, -8947.45
    full = a * n * n + b * n + c
    print(f"\nфит полной формулой при n=256: {full:.0f} против замера "
          f"799277 ({100 * full / 799277 - 100:+.2f}%)")
    print(f"фит только a*n^2 при n=256:    {a * n * n:.0f} "
          f"({100 * a * n * n / 799277 - 100:+.1f}%) - это те самые 45%")
    print(f"b*n = {b * n:.0f} поглощает 136*n*log2(n) = {nl}")

    # порог, за которым член n*log2(n) становится <10% квадратичного
    n_ = 16
    while 136 * m_of(n_) / (0.1 * 8 * n_) >= 1:  # 136*m*n / (0.8*n^2) >= 1
        n_ += 1
    print(f"\n136*n*log2(n) < 10% от 8n^2 начиная с n = {n_}")

    # структурно-осведомлённый фит: базис [n^2, n*m, n, 1] на всех точках
    ns = [float(r["n"]) for r in rows]
    ts = [float(r["toffoli"]) for r in rows]
    cols = [
        [v * v for v in ns],
        [v * m_of(int(v)) for v in ns],
        list(ns),
        [1.0] * len(ns),
    ]
    k = len(cols)
    A = [[sum(cols[i][t] * cols[j][t] for t in range(len(ns)))
          for j in range(k)] for i in range(k)]
    y = [sum(cols[i][t] * ts[t] for t in range(len(ns))) for i in range(k)]
    M = [row[:] + [y[i]] for i, row in enumerate(A)]
    for col in range(k):
        piv = max(range(col, k), key=lambda r_: abs(M[r_][col]))
        M[col], M[piv] = M[piv], M[col]
        for r_ in range(k):
            if r_ == col:
                continue
            factor = M[r_][col] / M[col][col]
            for kk in range(col, k + 1):
                M[r_][kk] -= factor * M[col][kk]
    beta = [M[i][k] / M[i][i] for i in range(k)]
    resid = max(abs(sum(beta[i] * cols[i][t] for i in range(k)) - ts[t]) / ts[t]
                for t in range(len(ns)))
    print(f"\nфит T = a*n^2 + b*n*m + c*n + d по всем 19 точкам:")
    print(f"  a={beta[0]:.4f} b={beta[1]:.4f} c={beta[2]:.4f} d={beta[3]:.1f}, "
          f"макс. отн. остаток {resid:.5f}")
    print(f"  (конструкция: a=8, b=136, c=-44; d ~ 128m^2-64m+45, вне базиса)")


if __name__ == "__main__":
    main()
