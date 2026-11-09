#!/usr/bin/env python3
"""Минимальный statevector-симулятор квантового компьютера на NumPy.

Базис: |q_{w-1} ... q_0>, старший бит - старший индекс амплитуды.
Вектор состояния - complex128 длины 2^w (16 байт на амплитуду).

Гейты: H, X, Z, фаза, CX, Toffoli, произвольный 1-кубитный унитар;
QFT по Z_{2^width} гейтами (H + controlled-фазы). Плюс точный QFT по
Z_m для произвольного m через numpy FFT (oracle-уровень) - эквивалентность
гейт-пути и FFT-пути проверяется в qip_shor_ecdlp.py.

Запуск: из корня репо (импорт через sys.path, как w1-эксперименты).
"""
import numpy as np


class Statevector:
    """Вектор состояния на w кубитах, complex128."""

    def __init__(self, width):
        self.w = width
        self.psi = np.zeros(2 ** width, dtype=np.complex128)
        self.psi[0] = 1.0

    def norm(self):
        return float(np.vdot(self.psi, self.psi).real)

    # ---------- гейты ----------

    def apply_1q(self, U, qubit):
        """Произвольный унитар 2x2 на кубите qubit."""
        psi = self.psi.reshape(-1, 2, 2 ** (self.w - qubit - 1))
        a = psi[:, 0, :].copy()
        b = psi[:, 1, :].copy()
        psi[:, 0, :] = U[0, 0] * a + U[0, 1] * b
        psi[:, 1, :] = U[1, 0] * a + U[1, 1] * b

    def h(self, qubit):
        r = 1.0 / np.sqrt(2.0)
        self.apply_1q(np.array([[r, r], [r, -r]], dtype=np.complex128),
                      qubit)

    def x(self, qubit):
        self.apply_1q(np.array([[0, 1], [1, 0]], dtype=np.complex128),
                      qubit)

    def phase(self, theta, qubit):
        U = np.diag([1.0, np.exp(1j * theta)]).astype(np.complex128)
        self.apply_1q(U, qubit)

    def cx(self, control, target):
        """CNOT: инверсия целевого бита там, где управление = 1.

        Общая реализация через маску индексов (как ccx/_swap): корректна
        для любой пары control != target. Ранняя reshape-версия была
        корректна только при target == control + 1.
        """
        mc = 1 << (self.w - control - 1)
        mt = 1 << (self.w - target - 1)
        idx = np.arange(self.psi.size)
        sel = (idx & mc) != 0
        other = idx[sel] ^ mt
        tmp = self.psi[idx[sel]].copy()
        self.psi[idx[sel]] = self.psi[other]
        self.psi[other] = tmp

    def ccx(self, c1, c2, target):
        """Toffoli через маску индексов (универсальный, но медленный путь)."""
        m1 = 1 << (self.w - c1 - 1)
        m2 = 1 << (self.w - c2 - 1)
        mt = 1 << (self.w - target - 1)
        idx = np.arange(self.psi.size)
        sel = (((idx & m1) != 0) & ((idx & m2) != 0))
        other = idx[sel] ^ mt
        tmp = self.psi[idx[sel]].copy()
        self.psi[idx[sel]] = self.psi[other]
        self.psi[other] = tmp

    def qft_pow2(self, qubits):
        """Гейт-точный QFT по Z_{2^len} на младших битах списка qubits.

        qubits[0] - старший значимый бит регистра (стандартная схема QFT:
        H на MSB, затем controlled-фазы вниз по цепочке).
        """
        L = len(qubits)
        for j in range(L):
            self.h(qubits[j])
            for k in range(j + 1, L):
                # фаза pi/2^(k-j) при контроле qubits[j]=1 на qubits[k]
                theta = np.pi / 2 ** (k - j)
                self._controlled_phase(qubits[j], qubits[k], theta)
        # реверс порядка битов регистра (стандартный завершающий swap-слой)
        for j in range(L // 2):
            self._swap(qubits[j], qubits[L - 1 - j])

    def _controlled_phase(self, control, target, theta):
        U = np.diag([1.0, np.exp(1j * theta)]).astype(np.complex128)
        m = 1 << (self.w - control - 1)
        idx = np.arange(self.psi.size)
        sel = (idx & m) != 0
        # применяем фазу к целевому кубиту внутри подмножества control=1:
        # через apply_1q на target, но только на выбранной половине -
        # проще умножением на диагональ по глобальным индексам
        tmask = 1 << (self.w - target - 1)
        factor = np.ones(idx.size, dtype=np.complex128)
        factor[sel & ((idx & tmask) != 0)] = U[1, 1]
        self.psi *= factor

    def _swap(self, qa, qb):
        ma = 1 << (self.w - qa - 1)
        mb = 1 << (self.w - qb - 1)
        idx = np.arange(self.psi.size)
        sel = (((idx & ma) != 0) & ((idx & mb) == 0))
        other = idx[sel] ^ ma ^ mb
        tmp = self.psi[idx[sel]].copy()
        self.psi[idx[sel]] = self.psi[other]
        self.psi[other] = tmp

    # ---------- измерение ----------

    def probabilities(self):
        return np.abs(self.psi) ** 2

    def measure(self, rng, qubits=None):
        """Выборочное измерение: возвращает (битовую строку qubits,
        обновлённый вектор после коллапса). qubits=None - весь регистр."""
        p = self.probabilities()
        cum = np.cumsum(p)
        total = cum[-1] if cum[-1] > 0 else 1.0
        idx = min(int(np.searchsorted(cum, rng.random() * total)),
                  p.size - 1)
        if qubits is None:
            self.psi = np.zeros_like(self.psi)
            self.psi[idx] = 1.0
            return idx, self.psi
        mask = sum(1 << (self.w - q - 1) for q in qubits)
        width = len(qubits)
        value = 0
        for pos, q in enumerate(qubits):
            if (idx >> (self.w - q - 1)) & 1:
                value |= 1 << (width - 1 - pos)
        keep = (np.arange(p.size) & mask) == ((idx & mask))
        self.psi = np.where(keep, self.psi, 0)
        nrm = np.linalg.norm(self.psi)
        if nrm > 0:
            self.psi /= nrm
        return value, self.psi


def qft_z_m_exact(vec, m):
    """Точный QFT по Z_m для произвольного m: |j> -> (1/sqrt(m)) *
    sum_k e^{2pi i jk/m} |k>. vec - длины m (или (m, m) для Z_m^2)."""
    arr = np.asarray(vec, dtype=np.complex128)
    if arr.ndim == 1:
        assert arr.shape[0] == m
        return np.fft.fft(arr, norm="ortho")
    assert arr.shape == (m, m)
    return np.fft.fft2(arr, norm="ortho")


def coset_state_posterior(n, x, z0):
    """Косет-состояние (a, b) после измерения точечного регистра.

    Измеренный результат z0 = a + b*x mod n; состояние равномерно по
    косету {(z0 - t*x, t) : t in Z_n} решётки L = span{(-x, 1)}.
    Возвращает (n, n) комплексную матрицу амплитуд.
    """
    psi = np.zeros((n, n), dtype=np.complex128)
    t = np.arange(n)
    rows = (z0 - t * x) % n
    psi[rows, t] = 1.0 / np.sqrt(n)
    return psi
