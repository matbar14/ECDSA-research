#!/usr/bin/env python3
"""Верификация движка qip_engine против матричных определений гейтов.

Гейт-путь Шора (qip_shor_ecdlp.py) доверяет движку полностью: рассинхрон
между конвенцией гейтов/измерения и оракулом ломает выборки тихо, без
ошибок (реальный случай: списки q_a/q_b/q_z, заданные в обратной
конвенции, давали 0/16 на двойственной решётке при внешне исправной
схеме). Поэтому каждый примитив проверяется здесь против явной матрицы:

  1. QFT-контур (qft_pow2) == DFT_{2^w} с положительной экспонентой,
     F[k, j] = e^{2*pi*i*j*k/N} / sqrt(N), w = 1..4;
  2. CX против матрицы перестановки для ВСЕХ пар control != target
     на w = 3 (ранняя reshape-версия была корректна только при
     target == control + 1, на остальных парах падала или молча
     делала частичный своп);
  3. Toffoli против матрицы для всех различных троек, w = 3;
  4. X, H, phase(θ) против их матриц на каждом кубите, w = 2;
  5. apply_1q со случайным унитаром против вложенной матрицы, w = 3;
  6. измерение: детерминированный коллапс базисного состояния,
     частичное измерение |+> на одном кубите (частоты и поддержка
     после коллапса), полное измерение равномерной суперпозиции;
  7. конвенция регистров оракула Шора: gather
     new[(a, b, z)] = old[(a, b, (z - f) mod n)], f = (a + b*x) mod n,
     переводит |a, b, 0> в |a, b, f>; поля читаются кубитами
     q_a = [0..w-1], q_b = [w..2w-1], q_z = [2w..3w-1] (MSB первым) -
     регрессионный тест бага конвенции;
  8. coset_state_posterior: ненулевые амплитуды точно на косете
     {(a, b): a + b*x = z0 (mod n)}, вес 1/sqrt(n).

Запуск из корня репо: python3 experiments/qip_test_engine.py
Выход: строки [PASS]/[FAIL], код возврата 0/1.
"""
import sys
from functools import reduce
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qip_engine import Statevector, coset_state_posterior  # noqa: E402

RNG = np.random.default_rng(20261109)


def circuit_matrix(w, apply_circuit):
    """Матрица схемы: колонка j = действие на базисное состояние |j>."""
    N = 2 ** w
    C = np.zeros((N, N), dtype=np.complex128)
    for j in range(N):
        sv = Statevector(w)
        sv.psi = np.zeros(N, dtype=np.complex128)
        sv.psi[j] = 1.0
        apply_circuit(sv)
        C[:, j] = sv.psi
    return C


def embed1(w, U, qubit):
    """Вложение 2x2 унитария на кубит qubit (кубит 0 - старший фактор)."""
    mats = [np.eye(2, dtype=np.complex128) for _ in range(w)]
    mats[qubit] = np.asarray(U, dtype=np.complex128)
    return reduce(np.kron, mats)


def dft_positive(N):
    j = np.arange(N)
    return np.exp(2j * np.pi * np.outer(j, j) / N) / np.sqrt(N)


def assert_close(A, B, tol=1e-10, msg=""):
    d = float(np.max(np.abs(A - B)))
    assert d < tol, f"{msg}: max|A-B| = {d:.3e} >= {tol:.0e}"


# ---------- тесты ----------

def test_qft_matrix():
    for w in (1, 2, 3, 4):
        C = circuit_matrix(w, lambda sv, w=w: sv.qft_pow2(list(range(w))))
        F = dft_positive(2 ** w)
        assert_close(C, F, msg=f"qft_pow2 w={w}")


def test_cx_all_pairs():
    w = 3
    N = 2 ** w
    mc_bit = lambda q: 1 << (w - 1 - q)
    for control in range(w):
        for target in range(w):
            if control == target:
                continue
            C = circuit_matrix(
                w, lambda sv, c=control, t=target: sv.cx(c, t))
            E = np.zeros((N, N))
            for i in range(N):
                j = i ^ mc_bit(target) if i & mc_bit(control) else i
                E[j, i] = 1.0
            assert_close(C, E, msg=f"cx({control},{target})")


def test_ccx_all_triples():
    w = 3
    N = 2 ** w
    bit = lambda q: 1 << (w - 1 - q)
    for c1 in range(w):
        for c2 in range(w):
            for t in range(w):
                if len({c1, c2, t}) != 3:
                    continue
                C = circuit_matrix(
                    w, lambda sv, a=c1, b=c2, q=t: sv.ccx(a, b, q))
                E = np.zeros((N, N))
                for i in range(N):
                    j = (i ^ bit(t)
                         if (i & bit(c1)) and (i & bit(c2)) else i)
                    E[j, i] = 1.0
                assert_close(C, E, msg=f"ccx({c1},{c2},{t})")


def test_1q_gates():
    w = 2
    r = 1.0 / np.sqrt(2.0)
    gates = {
        "x": (np.array([[0, 1], [1, 0]], dtype=complex),
              lambda sv, q: sv.x(q)),
        "h": (np.array([[r, r], [r, -r]], dtype=complex),
              lambda sv, q: sv.h(q)),
        "phase": (np.diag([1.0, np.exp(0.7j)]),
                  lambda sv, q: sv.phase(0.7, q)),
    }
    for name, (U, fn) in gates.items():
        for q in range(w):
            C = circuit_matrix(w, lambda sv, q=q, fn=fn: fn(sv, q))
            assert_close(C, embed1(w, U, q), msg=f"{name}(q{q})")


def test_apply_1q_random_unitary():
    w = 3
    G = RNG.normal(size=(2, 2)) + 1j * RNG.normal(size=(2, 2))
    Q, R = np.linalg.qr(G)
    Q = Q * (np.diag(R) / np.abs(np.diag(R)))  # фазовая нормализация
    assert_close(Q.conj().T @ Q, np.eye(2), tol=1e-12,
                 msg="QR не дал унитарий")
    for q in range(w):
        C = circuit_matrix(w, lambda sv, q=q: sv.apply_1q(Q, q))
        assert_close(C, embed1(w, Q, q), msg=f"apply_1q(q{q})")


def test_measurement():
    rng = np.random.default_rng(20261110)

    # базисное состояние: полный регист
    for _ in range(5):
        sv = Statevector(4)
        sv.psi = np.zeros(16, dtype=np.complex128)
        sv.psi[9] = 1.0
        val, psi = sv.measure(rng)
        assert val == 9, f"полное измерение |9> дало {val}"
        assert (np.nonzero(psi)[0] == [9]).all(), "коллапс не точечный"

    # частичное измерение кубита 0 в состоянии |+>|0>
    ones = 0
    trials = 400
    for _ in range(trials):
        sv = Statevector(2)
        sv.psi = np.zeros(4, dtype=np.complex128)
        sv.psi[0] = sv.psi[2] = 1.0 / np.sqrt(2.0)  # |00> + |10>
        val, psi = sv.measure(rng, qubits=[0])
        ones += val
        support = set(np.nonzero(psi)[0])
        subspace = {0, 1} if val == 0 else {2, 3}
        assert support <= subspace, \
            (f"после измерения q0={val} поддержка {support} "
             f"вне подпространства {subspace}")
        # у |+>|0> амплитуда только на |q1=0>: коллапс даёт базисное
        # состояние |00> при val=0 и |10> при val=1
        assert support == ({0} if val == 0 else {2}), \
            f"ожидаемое базисное состояние, получено {support}"
        assert abs(np.linalg.norm(psi) - 1) < 1e-12, \
            "норма после коллапса != 1"
    assert abs(ones / trials - 0.5) < 0.15, \
        f"частота q0=1: {ones / trials:.2f} != 0.5"

    # полное измерение равномерной суперпозиции w=2
    counts = np.zeros(4)
    trials = 800
    for _ in range(trials):
        sv = Statevector(2)
        sv.psi = np.ones(4, dtype=np.complex128) / 2.0
        val, _ = sv.measure(rng)
        counts[val] += 1
    assert np.max(np.abs(counts / trials - 0.25)) < 0.1, \
        f"частоты {counts / trials}"


def test_shor_oracle_convention():
    """Регрессионный тест бага конвенции: gather оракула + списки кубитов.

    Плоский индекс i = a*n^2 + b*n + z; кубит q соответствует биту
    (3w-1-q); q_a = [0..w-1], q_b = [w..2w-1], q_z = [2w..3w-1] MSB-первым.
    """
    n, x = 8, 5
    w = 3
    ab = np.arange(n * n).reshape(n, n)
    shift = (ab // n + (ab % n) * x) % n
    z_idx = (np.arange(n)[None, None, :]
             - shift[:, :, None]) % n
    flat_idx = ((np.arange(n)[:, None, None] * n * n)
                + (np.arange(n)[None, :, None] * n)
                + z_idx).reshape(-1)

    # |a=2, b=1, z=0> -> |2, 1, (2 + 1*5) mod 8 = 7>
    sv = Statevector(3 * w)
    sv.psi = np.zeros(2 ** (3 * w), dtype=np.complex128)
    sv.psi[2 * n * n + 1 * n + 0] = 1.0
    sv.psi = sv.psi[flat_idx]
    expect = 2 * n * n + 1 * n + 7
    got = int(np.nonzero(sv.psi)[0][0])
    assert got == expect, f"оракул: {got} != {expect}"

    # поля читаются своими списками кубитов (детерминированно)
    q_a = list(range(0, w))
    q_b = list(range(w, 2 * w))
    q_z = list(range(2 * w, 3 * w))
    rng = np.random.default_rng(1)
    val_a, _ = sv.measure(rng, qubits=q_a)
    val_b, _ = sv.measure(rng, qubits=q_b)
    val_z, _ = sv.measure(rng, qubits=q_z)
    assert (val_a, val_b, val_z) == (2, 1, 7), \
        f"поля (a,b,z) = ({val_a},{val_b},{val_z}) != (2,1,7)"


def test_coset_state_posterior():
    n, x, z0 = 8, 5, 3
    psi = coset_state_posterior(n, x, z0)
    rows, cols = np.nonzero(psi)  # 2D: кортеж (rows, cols)
    amplitudes = psi[rows, cols]
    expected = {(a, b) for b in range(n)
                for a in [(z0 - b * x) % n]}
    got = {(int(r), int(c)) for r, c in zip(rows, cols)}
    assert got == expected, f"косет {got} != {expected}"
    assert len(got) == n, f"мощность косета {len(got)} != {n}"
    assert_close(np.abs(amplitudes),
                 np.full(n, 1.0 / np.sqrt(n)), msg="веса косета")


TESTS = [
    ("1. qft_pow2 == DFT+ (w=1..4)", test_qft_matrix),
    ("2. cx все пары (w=3)", test_cx_all_pairs),
    ("3. ccx все тройки (w=3)", test_ccx_all_triples),
    ("4. x/h/phase против матриц", test_1q_gates),
    ("5. apply_1q случайный унитарий", test_apply_1q_random_unitary),
    ("6. измерение и коллапс", test_measurement),
    ("7. конвенция оракула Шора", test_shor_oracle_convention),
    ("8. coset_state_posterior", test_coset_state_posterior),
]


def main():
    failures = 0
    for name, fn in TESTS:
        try:
            fn()
            print(f"[PASS] {name}")
        except AssertionError as e:
            failures += 1
            print(f"[FAIL] {name}: {e}")
    total = len(TESTS)
    print(f"{total - failures}/{total} PASS")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
