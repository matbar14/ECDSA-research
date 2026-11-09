#!/usr/bin/env python3
"""Квантовое восстановление lambda и точек орбиты Lambda(G) на GLV-кривых.

Связка QIP-направления с GLV-темой репозитория. На кривой y^2 = x^3 + b
над p = 1 (mod 3) с простым #E = n эндоморфизм
phi: (x, y) -> (beta*x, y), beta^3 = 1 (beta != 1), действует как
умножение на lambda: phi(P) = lambda*P, lambda^2 + lambda + 1 = 0 (mod n).
Классическая работа репо построена на орбитах Lambda(P) = {lambda^r*P}:
канонические представители, файлы представителей, BSGS-таблицы.

Здесь алгоритм Шора для ECDLP (косет-QFT путь qip_shor_ecdlp.shor_fft)
запускается с целью H = lambda^r*G - точкой орбиты. Восстановленный
x = lambda^r (mod n) извлекается из нескольких выстрелов: каждый даёт
(u, v) с двойственной решётки L^ = {(u, u*x mod n)}, при u != 0 взаимно
просто с n x = v*u^-1 (mod n).

Проверяемое утверждение: квантовая атака не использует ни beta, ни
lambda, ни орбитную структуру - оракулу нужна только целевая точка.
Классическая экономика GLV/BSGS (которую меряет всё предыдущее
направление репо) против квантового атакующего не даёт ничего.

Оракул моделируется на уровне структуры косета: измерение точечного
регистра даёт z0 = a + b*x mod n (равномерно по косету решётки
L = span{(-x, 1)}), QFT по Z_n^2 концентрирует выборку на L^.
Групповые операции (построение H = lambda^r*G, классическая проверка
phi^r(G) == lambda^r*G) - препроцессинг/постобработка.

Запуск из корня репо: python3 experiments/qip_glv_orbit_quantum.py
Выход: results/qip_glv_orbit_quantum.json
"""
import argparse
import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qip_shor_ecdlp import ec_mul, find_point, glv_lambda, shor_fft  # noqa: E402


def find_beta(p):
    """Кубический корень из 1 по модулю p, beta != 1: beta = g^((p-1)/3)."""
    for g in range(2, 64):
        beta = pow(g, (p - 1) // 3, p)
        if beta != 1:
            return beta
    raise RuntimeError(f"beta не найдена для p={p}")


def phi(p, beta, P):
    """Эндоморфизм GLV-кривой: (x, y) -> (beta*x, y)."""
    return (beta * P[0] % p, P[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shots", type=int, default=16)
    ap.add_argument("--orbit-depth", type=int, default=3,
                    help="максимальная степень r в lambda^r*G")
    ap.add_argument("--seed", type=int, default=20261104)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    cache = Path("results/qip_toy_curves.json")
    curves = [c for c in json.loads(cache.read_text())
              if c["kind"] == "glv"]
    print(f"GLV-кривых в кеше: {len(curves)}")

    runs = []
    all_pass = True
    for c in curves:
        p, n = c["p"], c["order"]
        lam0 = glv_lambda(p, n)
        assert lam0 is not None, f"lambda не найдена для p={p}"
        beta = find_beta(p)
        # генератор: первая точка кривой (детерминирован, как в
        # qip_shor_ecdlp.find_toy_curves -> main); кеш хранит только
        # параметры поиска без точек
        G = find_point(p, c["a"], c["b"])
        assert G is not None, f"точка не найдена на p={p}"
        assert ec_mul(p, c["a"], n, G) is None, \
            f"n*G != O на p={p}"

        # дизамбигуация корня: glv_lambda возвращает первый корень
        # lam0 (lambda^2 + lambda + 1 = 0), но конкретному beta
        # отвечает либо lam0, либо lam0^2 - как в qip_shor_ecdlp.main
        lam = lam0
        if ec_mul(p, c["a"], lam, G) != phi(p, beta, G):
            lam2 = lam0 * lam0 % n
            assert ec_mul(p, c["a"], lam2, G) == phi(p, beta, G), \
                f"ни lambda, ни lambda^2 не дают phi(G) на p={p}"
            lam = lam2

        for r in range(1, args.orbit_depth + 1):
            x_target = pow(lam, r, n)
            H = ec_mul(p, c["a"], x_target, G)
            # утверждение орбиты: lambda^r*G == phi^r(G)
            Pr = G
            for _ in range(r):
                Pr = phi(p, beta, Pr)
            assert H == Pr, f"lambda^{r}*G != phi^{r}(G) на p={p}"

            t0 = time.perf_counter()
            samples, t_build, t_qft = shor_fft(c, x_target, rng,
                                               args.shots)
            wall = time.perf_counter() - t0
            rec = sum(1 for u, v in samples
                      if u and v * pow(u, -1, n) % n == x_target)
            u0 = sum(1 for u, v in samples if u == 0)
            on_lattice = all(v == u * x_target % n for u, v in samples)
            ok = rec >= 1
            all_pass &= ok
            mb = 16 * n * n / 2 ** 20
            runs.append({
                "p": p, "n": n, "bits": n.bit_length(),
                "lambda": lam, "r": r, "x target": x_target,
                "shots": args.shots,
                "recovered": f"{rec}/{args.shots}",
                "u zero": u0,
                "on dual lattice": on_lattice,
                "seconds": round(wall, 3),
                "statevector MB": round(mb, 1),
            })
            print(f"[orbit] p={p:5d} n={n:5d}: r={r} -> "
                  f"x=lambda^{r}={x_target:5d}: {rec}/{args.shots} "
                  f"за {wall:.2f} с ({mb:.1f} MB) - "
                  f"{'PASS' if ok else 'FAIL'}")

    report = {
        "description": "квантовое восстановление точек орбиты Lambda(G): "
                       "цель lambda^r*G, выборка с двойственной решётки "
                       "восстанавливает x = lambda^r без знания beta/lambda",
        "shots per run": args.shots,
        "orbit depth": args.orbit_depth,
        "runs": runs,
        "verdict": "PASS" if all_pass else "FAIL",
    }
    out = Path("results/qip_glv_orbit_quantum.json")
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"[saved] {out} - вердикт {report['verdict']}")


if __name__ == "__main__":
    main()
