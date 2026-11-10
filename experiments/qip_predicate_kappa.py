#!/usr/bin/env python3
"""kappa-модель предиката «k in S»: реверсивная схема, Тоффоли/глубина, масштабирование.

S routed-economy - не таблица (QRAM-модель снята, см. memo 10.11), а предикат:
    f(k) = 1  <=>  k - канонический представитель своей Lambda-орбиты
                   (орбита {k, lambda*k, lambda^2*k} mod N, lambda^3 = 1)
                   и сигнатура k (старшие b бит) попадает в интервал [lo, hi).
Предикат вычисляется арифметикой, а не ищется: 2 модульных умножения на
константу, 2 сравнения (канонизация - ровно 2 сравнения, порядок lambda = 3),
интервальная логика на b битах.

Схема (гейт-уровень, все ancilla очищаются, кроме K и выхода m):
- Cuccaro-рипл-аддер: 2 Toffoli + 4 CNOT на бит, 1 ancilla (cin). Выведен
  здесь и бит-верифицирован исчерпывающе; опубликованный счёт Cuccaro et al.
  2004 - 2n+1 Toffoli, совпадает;
- компаратор [X < Y] заёмной цепочкой: 5 Toffoli на бит;
- умножение на константу mod p: широкий аккумулятор (n + ceil(log2 n) бит),
  n inplace-сложений (аддер сам чистит tmp и cin), затем бинарное сокращение
  mod p за m = ceil(log2 n) шагов compare-and-conditional-subtract; каскад
  промежуточных регистров сокращения очищается внутри (rev шагов 1..m-1),
  широкий аккумулятор остаётся грязным и убирается реверсом у вызывающего;
- предикат: умножение lambda, умножение lambda^2, 2 компаратора [A < K],
  интервальные сравнения, AND; очистка - реверс под-схем в обратном
  хронологическом порядке (каждый sub - биекция, реверс восстанавливает).

Верификация:
- бит-уровень: все 2^n входов при n = 4..8 - K сохранён, m совпадает с
  классическим референсом на каждом входе, прочие провода нулеvые, обратный
  прогон возвращает исходное состояние (перестановка);
- движок qip_engine: таблицы истинности аддера и компаратора при n = 4 на
  базисных состояниях. Полный предикат и reduce на движке не прогоняются:
  >45 проводов (2^45 амплитуд) - их корректность покрывается бит-уровнем.

Замер: Toffoli/глубина/провода при n = 4..256. n - ширина арифметики
(log2 N): для secp256k1 берутся реальный порядок N и лямбды из
phase_signature_135bit.json; ширина окна W входит только в число итераций
Гровера, не в стоимость оракула.

Структурная модель счёта (точна на всех точках сетки, проверено
покомпонентно и на полном предикате):
    T(n) = 8*n*(n+m) + 64*(2m-1)*(n+m) + 20*n + 45,  m = ceil(log2 n),
т.е. 8*n^2 + 136*n*log2(n) - 44*n + 128*log2(n)^2 - 64*log2(n) + 45.
Квадратичный фит T(n) = a*n^2 + b*n + c - только сводка: его член b*n
поглощает структурный член 136*n*log2(n), поэтому a*n^2 по отдельности
не предиктор (при n=256: a*n^2 = 553k против замера 799277).

kappa = глубина * логический такт, tau из {10 ns, 100 ns, 1 mks}; глубина
берётся прямым замером, фит в kappa и W* не входит.
Префикс-вариант - механический пересчёт от замера: глубина аддера 6n
заменяется на Brent-Kung 2*log2(n)-1. Древесный параллелизм - оценка,
не замер (см. json). W* = ((pi/4)*rate*M*kappa)^2, rate/M - замеры из
h_curve_summary.json (canonical H rows, H = 400k).

Фальсификатор (пре-регистрирован в memo 10.11 «Дальше», коммит edd802d,
до кода): ветка O(N) - канонизация сравнением со всеми элементами орбиты
вместо двух соседей. В схеме канонизация - ровно 2 сравнения, ветка
не сработала; проверяется структурой схемы и бит-верификацией.

Запуск из корня репо: python3 experiments/qip_predicate_kappa.py
Выход: results/qip_predicate_kappa.json
"""
import json
import math
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qip_engine import Statevector  # noqa: E402

# реальный secp256k1: порядок N и лямбды (phase_signature_135bit.json)
N_SECP = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAM_SECP = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
LAM2_SECP = 0xac9c52b33fa3cf1f5ad9e3fd77ed9ba4a880b9fc8ec739c2e0cfc810b51283ce

B_SIG = 2          # бит сигнатуры (маршрутизация rho=3 укладывается в 2 бита)
LO, HI = 0, 2      # принятый интервал сигнатуры
TAUS = [(1e-8, "10 ns"), (1e-7, "100 ns"), (1e-6, "1 mks")]


# --------------------------------------------------------------------------
# гейт-IR: x / cx / ccx, бит-симуляция, счёт, глубина (ASAP), реверс
# --------------------------------------------------------------------------

class Circuit:
    def __init__(self):
        self.gates = []

    def x(self, q):
        self.gates.append(("x", q))

    def cx(self, c, t):
        self.gates.append(("cx", c, t))

    def ccx(self, c1, c2, t):
        self.gates.append(("ccx", c1, c2, t))

    def extend(self, other):
        self.gates.extend(other.gates)

    def reversed(self):
        r = Circuit()
        r.gates = list(reversed(self.gates))
        return r

    def counts(self):
        d = {"x": 0, "cx": 0, "ccx": 0}
        for g in self.gates:
            d[g[0]] += 1
        return d

    def toffoli(self):
        return sum(1 for g in self.gates if g[0] == "ccx")

    def depth(self):
        last = {}
        mx = 0
        for g in self.gates:
            t = max(last.get(q, 0) for q in g[1:]) + 1
            for q in g[1:]:
                last[q] = t
            if t > mx:
                mx = t
        return mx

    def run_bits(self, bits):
        for g in self.gates:
            if g[0] == "x":
                bits[g[1]] ^= 1
            elif g[0] == "cx":
                bits[g[2]] ^= bits[g[1]]
            else:
                bits[g[3]] ^= bits[g[1]] & bits[g[2]]

    def wires(self):
        return 1 + max((q for g in self.gates for q in g[1:]), default=-1)


class Wires:
    """Аллокатор виртуальных проводов (каждый alloc - свежие, нулевые)."""

    def __init__(self):
        self.next = 0

    def alloc(self, k):
        ws = list(range(self.next, self.next + k))
        self.next += k
        return ws

    def one(self):
        return self.alloc(1)[0]


# --------------------------------------------------------------------------
# примитивы
# --------------------------------------------------------------------------

def ripple_add(circ, A, B, cin):
    """B <- A + B (mod 2^len(B)); A сохранён, cin восстановлен.

    Cuccaro-структура: forward MAJ' (cx, cx, ccx) - перенос живёт на A[i]
    и модифицирует carry-провод; backward UMA (ccx, cx, cx) восстанавливает
    A[i] и carry-провод и дописывает бит суммы в B[i]. 2 Toffoli + 4 CNOT
    на бит, 1 ancilla. Верифицирован исчерпывающе (см. check_*).
    """
    n = len(B)
    for i in range(n):
        w = cin if i == 0 else A[i - 1]
        circ.cx(A[i], w)
        circ.cx(A[i], B[i])
        circ.ccx(B[i], w, A[i])
    for i in range(n - 1, -1, -1):
        w = cin if i == 0 else A[i - 1]
        circ.ccx(B[i], w, A[i])
        circ.cx(A[i], w)
        circ.cx(w, B[i])


def cmp_lt(X, Y, al):
    """(sub, flag): flag = [X < Y]; X, Y не изменяются.

    Заёмная цепочка b_{i+1} = maj(!X[i], Y[i], b_i); инверсия X[i] - X-гейтами
    вокруг Toffoli. 5 Toffoli + 2 CNOT + 4 X на бит. u, v очищаются в каждом
    бите; цепочка очищается реверсом sub (X, Y нетронуты).
    """
    n = len(X)
    s = Circuit()
    chain = al.alloc(n + 1)          # chain[0] = b_0 = 0; chain[n] = [X < Y]
    u, v = al.alloc(2)
    for i in range(n):
        s.x(X[i])
        s.ccx(X[i], Y[i], u)         # u = !X*Y
        s.ccx(X[i], chain[i], v)     # v = !X*b
        s.x(X[i])
        s.cx(u, chain[i + 1])
        s.cx(v, chain[i + 1])
        s.ccx(Y[i], chain[i], chain[i + 1])   # + Y*b -> maj(!X, Y, b)
        s.x(X[i])
        s.ccx(X[i], Y[i], u)         # uncompute u
        s.ccx(X[i], chain[i], v)     # uncompute v
        s.x(X[i])
    return s, chain[n]


def cond_sub(W, c, al):
    """(sub, R): R = W - c*[W >= c]; W не изменяется.

    1) компаратор f = [W < c] (цепочка читает W, W нетронут);
    2) g = [W >= c] = !f;
    3) D = W - c: копия W, затем D += !c + 1 (дополнение до двух);
    4) мультиплексор R = !g*W + g*D;
    5) очистка: g, h, реверс компаратора, реверс D-части.
    """
    nw = len(W)
    s = Circuit()
    ct = al.alloc(nw)                # |c>
    D = al.alloc(nw)
    R = al.alloc(nw)
    cin2 = al.one()
    for k in range(nw):
        if (c >> k) & 1:
            s.x(ct[k])
    cmp_s, f = cmp_lt(W, ct, al)
    s.extend(cmp_s)
    g = al.one()
    s.x(g)
    s.cx(f, g)                       # g = [W >= c]

    dpart = Circuit()
    for i in range(nw):
        dpart.cx(W[i], D[i])         # D = W
    for k in range(nw):
        dpart.x(ct[k])               # ct = !c
    dpart.x(cin2)
    ripple_add(dpart, ct, D, cin2)   # D = W + !c + 1 = W - c (mod 2^nw)
    dpart.x(cin2)
    for k in range(nw):
        dpart.x(ct[k])               # ct = c
    s.extend(dpart)

    h = al.one()
    s.x(h)
    s.cx(g, h)                       # h = !g
    for i in range(nw):
        s.ccx(h, W[i], R[i])
        s.ccx(g, D[i], R[i])         # R = !g*W + g*D

    s.cx(g, h)
    s.x(h)                           # h = 0
    s.cx(f, g)
    s.x(g)                           # g = 0
    s.extend(cmp_s.reversed())       # f, chain = 0 (W, ct нетронуты)
    s.extend(dpart.reversed())       # D = 0
    for k in range(nw):              # снять материализацию константы ct -> 0
        if (c >> k) & 1:
            s.x(ct[k])
    return s, R


def build_reduce(ACC, p, al, m):
    """(sub, out): out = ACC mod p; ACC не изменяется.

    Бинарное сокращение: шаги j = m-1..0, R_j = R_{j+1} - p*2^j*[R_{j+1} >= p*2^j],
    R_m := ACC. Каскад R_1..R_{m-1} неизбежно жив до конца (реверс шага j
    требует R_{j+1}); очищается реверсом шагов 1..m-1 в хронологическом
    порядке. ACC остаётся грязным (широкое произведение).
    """
    red = Circuit()
    steps = []
    inp = ACC
    for j in range(m - 1, -1, -1):
        s, Rj = cond_sub(inp, p << j, al)
        red.extend(s)
        steps.append((j, s))
        inp = Rj
    out = inp
    for j, s in reversed(steps):     # j = 0, 1, ..., m-1
        if j >= 1:
            red.extend(s.reversed())  # R_j -> 0; вход R_{j+1} ещё жив
    return red, out


def build_mult(K, coeffs, p, al):
    """(sub, out, ACC): out = (sum_i K[i]*coeffs[i]) mod p; K не изменяется.

    Фаза 1: широкий аккумулятор ACC (n + m бит, m = ceil(log2 n)),
    n inplace-сложений tmp_i = K[i]*c_i (аддер сам восстанавливает tmp, cin).
    Фаза 2: reduce_mod. ACC остаётся грязным - очистка реверсом sub у
    вызывающего (sub - биекция).
    """
    n = len(K)
    m = max(1, (n - 1).bit_length())
    nw = n + m
    s = Circuit()
    ACC = al.alloc(nw)
    tmp = al.alloc(nw)               # старшие m бит всегда 0 (паддинг аддера)
    cin = al.one()
    for i in range(n):
        ci = coeffs[i]
        for j in range(n):
            if (ci >> j) & 1:
                s.cx(K[i], tmp[j])
        ripple_add(s, tmp, ACC, cin)
        for j in range(n):
            if (ci >> j) & 1:
                s.cx(K[i], tmp[j])
    red, out = build_reduce(ACC, p, al, m)
    s.extend(red)
    return s, out, ACC


def build_predicate(n, p, lam, lam2, b_sig=B_SIG, lo=LO, hi=HI):
    """(C, m_wire, K): m = f(k); все прочие провода очищены.

    f(k) = [k - канонический представитель орбиты] AND [sig(k) in [lo, hi)],
    канонизация: k <= lambda*k и k <= lambda^2*k - ровно 2 сравнения.
    """
    C = Circuit()
    al = Wires()
    K = al.alloc(n)
    coeffs1 = [lam * pow(2, i, p) % p for i in range(n)]
    coeffs2 = [lam2 * pow(2, i, p) % p for i in range(n)]

    m1, A1, _ = build_mult(K, coeffs1, p, al)
    C.extend(m1)
    m2, A2, _ = build_mult(K, coeffs2, p, al)
    C.extend(m2)

    c1s, f1 = cmp_lt(A1[:n], K, al)  # f1 = [A1 < K]; A1 < p < 2^n - старшие m бит 0
    C.extend(c1s)
    c2s, f2 = cmp_lt(A2[:n], K, al)
    C.extend(c2s)

    sig = K[n - b_sig:]
    lo_w = al.alloc(b_sig)
    hi_w = al.alloc(b_sig)
    xset = Circuit()
    for k in range(b_sig):
        if (lo >> k) & 1:
            xset.x(lo_w[k])
        if (hi >> k) & 1:
            xset.x(hi_w[k])
    los, flo = cmp_lt(sig, lo_w, al)   # flo = [sig < lo]
    his, fhi = cmp_lt(sig, hi_w, al)   # fhi = [sig < hi]
    xset.extend(los)
    xset.extend(his)
    C.extend(xset)

    # AND: c1 = !f1 = [K <= A1], c2 = !f2, ta = !flo, tb = ta&fhi,
    #      t1 = c1&c2, m = t1&tb
    c1 = al.one()
    c2 = al.one()
    ta = al.one()
    tb = al.one()
    t1 = al.one()
    mm = al.one()
    C.x(c1)
    C.cx(f1, c1)
    C.x(c2)
    C.cx(f2, c2)
    C.x(ta)
    C.cx(flo, ta)
    C.ccx(ta, fhi, tb)
    C.ccx(c1, c2, t1)
    C.ccx(t1, tb, mm)

    # очистка в обратном хронологическом порядке
    C.ccx(c1, c2, t1)                # t1 = 0 (c1, c2 ещё живы)
    C.cx(f1, c1)
    C.x(c1)                          # c1 = 0
    C.cx(f2, c2)
    C.x(c2)                          # c2 = 0
    C.ccx(ta, fhi, tb)               # tb = 0 (ta ещё жив)
    C.cx(flo, ta)
    C.x(ta)                          # ta = 0
    C.extend(xset.reversed())        # rev his + rev los + снятие X-set lo_w/hi_w
    C.extend(c2s.reversed())         # f2, цепочка, A2 нетронуты
    C.extend(c1s.reversed())
    C.extend(m2.reversed())          # A2, ACC2, каскад -> 0
    C.extend(m1.reversed())
    return C, mm, K


def f_classical(x, n, p, lam, lam2, b_sig=B_SIG, lo=LO, hi=HI):
    """Классический референс предиката."""
    y = lam * x % p
    z = lam2 * x % p
    canonical = (x <= y) and (x <= z)
    sig = x >> (n - b_sig)
    return 1 if (canonical and lo <= sig < hi) else 0


# --------------------------------------------------------------------------
# параметры: игрушечные кривые и реальный secp256k1
# --------------------------------------------------------------------------

def is_prime_mr(v, rounds=24):
    """Миллер-Рабин с детерминированным seedом на кандидата (воспроизводимость)."""
    if v < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if v % sp == 0:
            return v == sp
    rng = random.Random(20261110 + v)
    d = v - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for _ in range(rounds):
        a = rng.randrange(2, v - 1)
        x = pow(a, d, v)
        if x == 1 or x == v - 1:
            continue
        for _ in range(r - 1):
            x = x * x % v
            if x == v - 1:
                break
        else:
            return False
    return True


def toy_params(n):
    """(p, lambda): простое p ~ 2^n, p = 1 (mod 3), lambda порядка 3 mod p.

    Случайный поиск кандидатов + Миллер-Рабин; lambda = a^((p-1)/3) mod p
    для случайного a (если != 1, порядок ровно 3). Перебор lambda подряд
    и пробное деление невозможны при n >= 48 (корни куба разрежены,
    sqrt(p) делений - часы/дни).
    """
    rng = random.Random(20261110)
    while True:
        p = rng.randrange(2 ** (n - 1) + 1, 2 ** n) | 1
        if p % 3 != 1 or not is_prime_mr(p):
            continue
        for _ in range(100):
            a = rng.randrange(2, p - 1)
            lam = pow(a, (p - 1) // 3, p)
            if lam != 1:
                assert pow(lam, 3, p) == 1
                return p, lam
        raise RuntimeError("lambda не найдена за 100 попыток")


def params_for(n):
    if n == 256:
        assert (N_SECP - 1) % 3 == 0
        assert pow(LAM_SECP, 3, N_SECP) == 1 and LAM_SECP != 1
        assert LAM2_SECP == LAM_SECP * LAM_SECP % N_SECP
        return N_SECP, LAM_SECP, LAM2_SECP
    p, lam = toy_params(n)
    return p, lam, lam * lam % p


# --------------------------------------------------------------------------
# верификация
# --------------------------------------------------------------------------

def init_bits(C, K, x):
    bits = {q: 0 for q in range(C.wires())}
    for j in range(len(K)):
        bits[K[j]] = (x >> j) & 1
    return bits


def check_reduce_bits(n, p, lam):
    """Все значения ACC из области определения (ACC < n*p, вне её частное
    q >= 2^m и m шагов сокращения не хватает - контракт аккумулятора)."""
    m = max(1, (n - 1).bit_length())
    nw = n + m
    fails = 0
    total = 0
    for v in range(n * p):
        al = Wires()
        ACC = al.alloc(nw)
        red, out = build_reduce(ACC, p, al, m)
        rev = red.reversed()
        bits = init_bits(red, ACC, v)
        red.run_bits(bits)
        got = sum(bits[out[j]] << j for j in range(nw))
        if got != v % p:
            fails += 1
        total += 1
        rev.run_bits(bits)
        # после реверса все провода, кроме ACC (вход reduce, по дизайну
        # не очищается - чистится реверсом у вызывающего), должны быть нулевыми
        if any(bits[q] for q in range(red.wires()) if q not in ACC):
            fails += 1
    return fails, total


def check_predicate_bits(n, p, lam, lam2, exhaustive_rev=True):
    """Все 2^n входов: m == f_cl(x), K сохранён, прочие провода нулевые,
    обратный прогон возвращает исходное состояние."""
    C, mm, K = build_predicate(n, p, lam, lam2)
    rev = C.reversed()
    wires = C.wires()
    fails = []
    for x in range(2 ** n):
        bits = init_bits(C, K, x)
        C.run_bits(bits)
        got = bits[mm]
        exp = f_classical(x, n, p, lam, lam2)
        k_back = sum(bits[K[j]] << j for j in range(n))
        dirty = [q for q in range(wires) if bits[q] and q not in K and q != mm]
        if got != exp or k_back != x or dirty:
            fails.append((x, got, exp, k_back, dirty[:5]))
        rev.run_bits(bits)
        # после реверса остались только биты K (значение x)
        clean = sum(bits[q] for q in range(wires))
        if clean != bin(x).count("1"):
            fails.append((x, "reverse", clean))
    return fails, 2 ** n, wires, C.toffoli(), C.depth()


def check_engine_primitives(n=4):
    """Таблицы истинности аддера и компаратора на движке (базисные состояния)."""
    p, lam = toy_params(n)
    res = {"n": n, "p": p, "lambda": lam}

    # аддер: все (a, b) - бит-уровень и движок
    fails_bit = 0
    add_gates = None
    add_wires = None
    for a in range(2 ** n):
        for b in range(2 ** n):
            al = Wires()
            C = Circuit()
            A = al.alloc(n)
            B = al.alloc(n)
            cin = al.one()
            ripple_add(C, A, B, cin)
            add_gates = C
            add_wires = C.wires()
            bits = {q: 0 for q in range(add_wires)}
            for j in range(n):
                bits[A[j]] = (a >> j) & 1
                bits[B[j]] = (b >> j) & 1
            C.run_bits(bits)
            got = sum(bits[B[j]] << j for j in range(n))
            a_back = sum(bits[A[j]] << j for j in range(n))
            if got != (a + b) % 2 ** n or a_back != a or bits[cin] != 0:
                fails_bit += 1
    res["adder bit-level pairs"] = 2 ** (2 * n)
    res["adder bit-level fails"] = fails_bit

    # движок: каждый (a, b) - базисное состояние, результат - базисное
    fails_eng = 0
    for a in range(2 ** n):
        for b in range(2 ** n):
            al = Wires()
            C = Circuit()
            A = al.alloc(n)
            B = al.alloc(n)
            cin = al.one()
            ripple_add(C, A, B, cin)
            w = C.wires()
            sv = Statevector(w)
            start = 0
            target = 0
            bits = {}
            for j in range(n):
                bits[A[j]] = (a >> j) & 1
                bits[B[j]] = (b >> j) & 1
            bits[cin] = 0
            for q, v in bits.items():
                start |= v << (w - 1 - q)
            exp_bits = dict(bits)
            for j in range(n):
                exp_bits[B[j]] = ((a + b) >> j) & 1
            for q, v in exp_bits.items():
                target |= v << (w - 1 - q)
            sv.psi[:] = 0
            sv.psi[start] = 1.0
            for g in C.gates:
                if g[0] == "x":
                    sv.x(g[1])
                elif g[0] == "cx":
                    sv.cx(g[1], g[2])
                else:
                    sv.ccx(g[1], g[2], g[3])
            dev = abs(abs(sv.psi[target]) - 1.0)
            rest = abs(sv.psi) ** 2
            rest[target] = 0.0
            if dev > 1e-9 or rest.max() > 1e-9:
                fails_eng += 1
    res["adder engine fails"] = fails_eng

    # компаратор: все (x, y) - бит-уровень; движок - выборка
    fails_cmp = 0
    cmp_gates = None
    cmp_wires = None
    for x in range(2 ** n):
        for y in range(2 ** n):
            al = Wires()
            C = Circuit()
            X = al.alloc(n)
            Y = al.alloc(n)
            cmp_s, f = cmp_lt(X, Y, al)
            C.extend(cmp_s)
            cmp_gates = C
            cmp_wires = C.wires()
            bits = {q: 0 for q in range(cmp_wires)}
            for j in range(n):
                bits[X[j]] = (x >> j) & 1
                bits[Y[j]] = (y >> j) & 1
            C.run_bits(bits)
            got = bits[f]
            exp = 1 if x < y else 0
            cmp_bad = got != exp
            # промежуточные провода цепочки очищаются реверсом sub - это часть
            # дизайна (значение флага уже снято); проверяем что реверс чистит
            C.reversed().run_bits(bits)
            if any(bits[q] for q in range(cmp_wires) if q not in X and q not in Y):
                cmp_bad = True
            x_back = sum(bits[X[j]] << j for j in range(n))
            y_back = sum(bits[Y[j]] << j for j in range(n))
            if x_back != x or y_back != y:
                cmp_bad = True
            if cmp_bad:
                fails_cmp += 1
    res["comparator bit-level pairs"] = 2 ** (2 * n)
    res["comparator bit-level fails"] = fails_cmp

    fails_eng_cmp = 0
    for x in range(2 ** n):
        for y in range(2 ** n):
            al = Wires()
            C = Circuit()
            X = al.alloc(n)
            Y = al.alloc(n)
            cmp_s, f = cmp_lt(X, Y, al)
            C.extend(cmp_s)
            w = C.wires()
            sv = Statevector(w)
            bits = {q: 0 for q in range(w)}
            for j in range(n):
                bits[X[j]] = (x >> j) & 1
                bits[Y[j]] = (y >> j) & 1
            start = sum(v << (w - 1 - q) for q, v in bits.items())
            exp_bits = dict(bits)
            exp_bits[f] = 1 if x < y else 0
            target = sum(v << (w - 1 - q) for q, v in exp_bits.items())
            sv.psi[:] = 0
            sv.psi[start] = 1.0
            for g in C.gates:
                if g[0] == "x":
                    sv.x(g[1])
                elif g[0] == "cx":
                    sv.cx(g[1], g[2])
                else:
                    sv.ccx(g[1], g[2], g[3])
            # промежуточные провода цепочки очищаются реверсом - гоняем реверс
            # и требуем возврат в исходное базисное состояние
            for g in reversed(C.gates):
                if g[0] == "x":
                    sv.x(g[1])
                elif g[0] == "cx":
                    sv.cx(g[1], g[2])
                else:
                    sv.ccx(g[1], g[2], g[3])
            dev = abs(abs(sv.psi[start]) - 1.0)
            rest = abs(sv.psi) ** 2
            rest[start] = 0.0
            if dev > 1e-9 or rest.max() > 1e-9:
                fails_eng_cmp += 1
    res["comparator engine fails"] = fails_eng_cmp
    return res


# --------------------------------------------------------------------------
# замер и масштабирование
# --------------------------------------------------------------------------

N_GRID = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 24, 32, 48, 64, 96, 128, 192, 256]


def measure_all():
    rows = []
    for n in N_GRID:
        p, lam, lam2 = params_for(n)
        t0 = time.perf_counter()
        C, mm, K = build_predicate(n, p, lam, lam2)
        cnt = C.counts()
        row = {
            "n": n,
            "modulus bits": n,
            "toffoli": cnt["ccx"],
            "cnot": cnt["cx"],
            "x gates": cnt["x"],
            "gates total": len(C.gates),
            "depth": C.depth(),
            "wires": C.wires(),
            "toffoli / n^2": round(cnt["ccx"] / n / n, 2),
            "build seconds": round(time.perf_counter() - t0, 2),
        }
        rows.append(row)
        print(f"[n={n:3d}] Toffoli {cnt['ccx']:9d} | CNOT {cnt['cx']:9d} | "
              f"глубина {C.depth():9d} | провода {C.wires():6d} | "
              f"билд {row['build seconds']:5.2f} с")
    return rows


def fit_quadratic(rows):
    """T(n) = a*n^2 + b*n + c (МНК) на n >= 24: на малых n скачки
    m = ceil(log2 n) (2 -> 3 -> 4 ...) и константные накладки маскируют
    квадратичный член (на 4..13 фит даёт отрицательный a - мусор).

    ВНИМАНИЕ: a - не предиктор по отдельности. Строение схемы даёт точную
    модель T(n) = 8*n*(n+m) + 64*(2m-1)*(n+m) + 20*n + 45, m = ceil(log2 n),
    т.е. 8n^2 + 136*n*log2(n) + O(n). Член 136*n*log2(n) фит поглощает в
    b*n (при n=256: 136*8*256 = 278528, 34.8% счёта; b*n = 258005), поэтому
    a*n^2 занижает замер на 31%, а полная фит-формула попадает с точностью
    ~0.4%. См. fit_structural - МНК по структурному базису восстанавливает
    коэффициенты конструкции."""
    pts = [r for r in rows if r["n"] >= 24]
    ns = [float(r["n"]) for r in pts]
    ts = [float(r["toffoli"]) for r in pts]
    # базис [1, n, n^2], beta = (c, b, a)
    A = [[float(len(ns)), sum(ns), sum(v * v for v in ns)],
         [sum(ns), sum(v * v for v in ns), sum(v ** 3 for v in ns)],
         [sum(v * v for v in ns), sum(v ** 3 for v in ns), sum(v ** 4 for v in ns)]]
    y = [sum(ts),
         sum(t * v for t, v in zip(ts, ns)),
         sum(t * v * v for t, v in zip(ts, ns))]
    # Гаусс
    Mm = [row[:] + [y[i]] for i, row in enumerate(A)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda r: abs(Mm[r][col]))
        Mm[col], Mm[piv] = Mm[piv], Mm[col]
        for r in range(3):
            if r == col:
                continue
            factor = Mm[r][col] / Mm[col][col]
            for k in range(col, 4):
                Mm[r][k] -= factor * Mm[col][k]
    c_, b_, a_ = (Mm[0][3] / Mm[0][0], Mm[1][3] / Mm[1][1], Mm[2][3] / Mm[2][2])
    residuals = []
    for r in pts:
        pred = a_ * r["n"] ** 2 + b_ * r["n"] + c_
        residuals.append(abs(pred - r["toffoli"]) / r["toffoli"])
    return {"a": round(a_, 4), "b": round(b_, 3), "c": round(c_, 2),
            "max rel residual": round(max(residuals), 5),
            "points": "n = 24..256 (малые n исключены: скачки m маскируют n^2)",
            "caveat": "a*n^2 - не предиктор: член 136*n*log2(n) поглощён "
                      "b*n; предиктор - структурная модель (см. structural model)"}


def fit_structural(rows):
    """МНК по структурному базису [n^2, n*m, n, m^2, m, 1] на ВСЕХ точках.

    Структурная формула T(n) = 8*n*(n+m) + 64*(2m-1)*(n+m) + 20*n + 45
    линейна по этому базису с коэффициентами (8, 136, -44, 128, -64, 45);
    фит должен их восстановить с нулевым остатком."""
    ns = [float(r["n"]) for r in rows]
    ts = [float(r["toffoli"]) for r in rows]
    ms = [float(max(1, (int(v) - 1).bit_length())) for v in ns]
    cols = [[v * v for v in ns],
            [v * mm for v, mm in zip(ns, ms)],
            list(ns),
            [mm * mm for mm in ms],
            ms,
            [1.0] * len(ns)]
    k = len(cols)
    A = [[sum(cols[i][t] * cols[j][t] for t in range(len(ns)))
          for j in range(k)] for i in range(k)]
    y = [sum(cols[i][t] * ts[t] for t in range(len(ns))) for i in range(k)]
    M = [row[:] + [y[i]] for i, row in enumerate(A)]
    for col in range(k):
        piv = max(range(col, k), key=lambda r: abs(M[r][col]))
        M[col], M[piv] = M[piv], M[col]
        for r in range(k):
            if r == col:
                continue
            factor = M[r][col] / M[col][col]
            for kk in range(col, k + 1):
                M[r][kk] -= factor * M[col][kk]
    beta = [M[i][k] / M[i][i] for i in range(k)]
    resid = max(abs(sum(beta[i] * cols[i][t] for i in range(k)) - ts[t]) / ts[t]
                for t in range(len(ns)))
    names = ["a (n^2)", "b (n*m)", "c (n)", "d (m^2)", "e (m)", "f (1)"]
    construction = [8.0, 136.0, -44.0, 128.0, -64.0, 45.0]
    dev = max(abs(be - co) for be, co in zip(beta, construction))
    return {"coefficients": {nm: round(be, 6) for nm, be in zip(names, beta)},
            "construction": {nm: co for nm, co in zip(names, construction)},
            "max |fit-construction|": round(dev, 6),
            "max rel residual": float(f"{resid:.2e}"),
            "points": "все 19 точек сетки n = 4..256",
            "verdict": "фит по структурному базису восстанавливает "
                       "коэффициенты конструкции: модель точна, квадратичная "
                       "сводка (a=8.44) - не предиктор по отдельности"}


def decompose_point(n, p, lam):
    """Покомпонентный разбор Тоффоли-счёта при ширине n.

    Строит компоненты отдельно и сверяет каждый с формулой строения:
    умножение = фаза 1 (n сложений по 2(n+m)) + сокращение ((2m-1) cond_sub
    по 16(n+m)); предикат = 4*умножение (2 умножения x очистка) + 4*компа-
    ратор [A<K] по 5n (2 прямых + 2 в очистке) + интервальная логика
    (2 компаратора на B_SIG битах, прямые + реверс) + AND-дерево (5).
    Счёт Тоффоли не зависит от значений констант (только от ширин), поэтому
    сверка компонент - сверка структуры, значения p/lam роли не играют.
    """
    m = max(1, (n - 1).bit_length())
    nw = n + m
    al = Wires()
    K = al.alloc(n)
    ph = Circuit()
    ACC = al.alloc(nw)
    tmp = al.alloc(nw)
    cin = al.one()
    coeffs = [lam * pow(2, i, p) % p for i in range(n)]
    for i in range(n):
        ci = coeffs[i]
        for j in range(n):
            if (ci >> j) & 1:
                ph.cx(K[i], tmp[j])
        ripple_add(ph, tmp, ACC, cin)
        for j in range(n):
            if (ci >> j) & 1:
                ph.cx(K[i], tmp[j])
    mult_phase1 = ph.toffoli()
    red, _ = build_reduce(ACC, p, al, m)
    mult_reduce = red.toffoli()
    mult_one = mult_phase1 + mult_reduce

    A1 = al.alloc(n)
    cs, _ = cmp_lt(A1, K, al)
    cmp_one = cs.toffoli()

    sig = al.alloc(B_SIG)
    lo_w = al.alloc(B_SIG)
    hi_w = al.alloc(B_SIG)
    los, _ = cmp_lt(sig, lo_w, al)
    his, _ = cmp_lt(sig, hi_w, al)
    interval_cmp = los.toffoli() + his.toffoli()

    and_tree = Circuit()
    fhi = al.one()          # флаг интервального компаратора (значение неважно)
    c1, c2, ta, tb, t1, mm = (al.one() for _ in range(6))
    and_tree.ccx(ta, fhi, tb)
    and_tree.ccx(c1, c2, t1)
    and_tree.ccx(t1, tb, mm)
    and_tree.ccx(c1, c2, t1)   # очистка
    and_tree.ccx(ta, fhi, tb)  # очистка
    and_t = and_tree.toffoli()

    checks = {
        "mult phase-1 (x1)": (mult_phase1, 2 * n * nw),
        "mult reduce (x1)": (mult_reduce, 16 * (2 * m - 1) * nw),
        "comparator [A<K] (x1)": (cmp_one, 5 * n),
        "interval cmps (both)": (interval_cmp, 2 * 5 * B_SIG),
        "AND tree": (and_t, 5),
    }
    comp = {}
    ok_all = True
    for name, (got, exp) in checks.items():
        comp[name] = {"measured": got, "formula": exp, "match": got == exp}
        ok_all &= got == exp
    pred = 4 * mult_one + 4 * cmp_one + 2 * interval_cmp + and_t
    comp["predicate (total)"] = {"measured": pred,
                                 "formula": "4*mult + 4*cmp + 2*interval + AND",
                                 "match": ok_all}
    return comp


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def load_rate_m():
    hs = json.loads(Path("results/h_curve_summary.json").read_text())
    for pt in hs["points"]:
        if pt["heads"] == 400000:
            ch = pt["throughput"]["canonical H rows"]
            return ch["giant steps/second"], ch["table rows stored"]
    raise RuntimeError("canonical H rows @400k не найдены")


def main():
    print("== Верификация примитивов (движок + бит-уровень) ==")
    eng = check_engine_primitives(4)
    print(f"аддер: бит-уровень {eng['adder bit-level pairs']} пар, "
          f"фейлов {eng['adder bit-level fails']}; движок фейлов "
          f"{eng['adder engine fails']}")
    print(f"компаратор: бит-уровень {eng['comparator bit-level pairs']} пар, "
          f"фейлов {eng['comparator bit-level fails']}; движок фейлов "
          f"{eng['comparator engine fails']}")

    print("== Верификация reduce (бит-уровень, все значения ACC) ==")
    red_checks = {}
    for n in (4, 5, 6):
        p, lam, _ = params_for(n)
        fails, total = check_reduce_bits(n, p, lam)
        red_checks[str(n)] = {"values": total, "fails": fails}
        print(f"n={n}: {total} значений ACC, фейлов {fails}")

    print("== Верификация полного предиката (бит-уровень, все 2^n входов) ==")
    pred_checks = {}
    for n in (4, 5, 6, 7, 8):
        p, lam, lam2 = params_for(n)
        fails, total, wires, toff, depth = check_predicate_bits(n, p, lam, lam2)
        pred_checks[str(n)] = {
            "inputs": total, "fails": len(fails),
            "wires": wires, "toffoli": toff, "depth": depth,
        }
        print(f"n={n}: {total} входов, фейлов {len(fails)} | "
              f"провода {wires} | Toffoli {toff} | глубина {depth}")
        if fails:
            for f in fails[:5]:
                print("  FAIL:", f)

    print("== Замер Toffoli/глубина/провода, n = 4..256 ==")
    rows = measure_all()

    fit = fit_quadratic(rows)
    print(f"фит T(n) = {fit['a']}*n^2 + {fit['b']}*n + {fit['c']} "
          f"({fit['points']}), макс. относительный остаток "
          f"{fit['max rel residual']}")

    print("== Структурная модель и покомпонентный разбор ==")
    # точность структурной формулы на ВСЕХ точках сетки
    exact = 0
    for r in rows:
        n = r["n"]
        m_ = max(1, (n - 1).bit_length())
        exact += (8 * n * (n + m_) + 64 * (2 * m_ - 1) * (n + m_)
                  + 20 * n + 45 == r["toffoli"])
    print(f"структурная формула T(n) = 8n(n+m)+64(2m-1)(n+m)+20n+45: "
          f"точна на {exact}/{len(rows)} точках сетки n=4..256")

    fit_st = fit_structural(rows)
    print(f"фит по структурному базису: макс. |коэф-конструкция| "
          f"{fit_st['max |fit-construction|']}, остаток "
          f"{fit_st['max rel residual']}")

    dec = []
    for n in (8, 24, 64, 128, 192, 256):
        p, lam, _ = params_for(n)
        comp = decompose_point(n, p, lam)
        meas = next(r["toffoli"] for r in rows if r["n"] == n)
        pred = comp["predicate (total)"]["measured"]
        match = pred == meas
        dec.append({"n": n, "components": comp,
                    "predicate measured": meas, "match": match})
        print(f"n={n:3d}: разбор {pred:8d} против замера {meas:8d} - "
              f"{'совпало' if match else 'MISMATCH'}")

    # глубина аддера на бит (для префикс-пересчёта)
    al = Wires()
    C = Circuit()
    A = al.alloc(256)
    B = al.alloc(256)
    cin = al.one()
    ripple_add(C, A, B, cin)
    d_add = C.depth()
    d_add_bit = d_add / 256
    bk = 2 * math.log2(256) - 1
    factor = bk / d_add
    print(f"аддер n=256: глубина {d_add} ({d_add_bit:.2f} на бит); "
          f"Brent-Kung {bk:.0f}; коэффициент префикс-пересчёта "
          f"{factor:.5f}")

    n256 = next(r for r in rows if r["n"] == 256)
    rate, M = load_rate_m()
    rm = rate * M
    log2_rows = []
    for tau, tau_name in TAUS:
        for kind, dep in (("эталонная схема", n256["depth"]),
                          ("префикс-суммы", n256["depth"] * factor),
                          ("древесный параллелизм (оценка)",
                           6 * math.log2(256) ** 2)):
            kappa = dep * tau
            w_star = ((math.pi / 4) * rm * kappa) ** 2
            log2_rows.append({
                "tau": tau_name, "variant": kind,
                "kappa seconds": float(f"{kappa:.3e}"),
                "log2 W*": round(math.log2(w_star), 1) if w_star > 0 else None,
            })
            print(f"tau={tau_name:7s} | {kind:32s} | kappa = {kappa:.3e} с | "
                  f"W* = 2^{math.log2(w_star):.1f}")

    # контекст итераций: |S| = 1.167*W/M (routed), итерации (pi/4)*sqrt(M/1.167)
    iter_mark = (math.pi / 4) * math.sqrt(M / 1.167)
    print(f"итераций маркировки (pi/4)*sqrt(M/1.167) = {iter_mark:.0f} "
          f"(константа, не зависит от W)")

    # W* при ожиданиях из постановки
    ref = []
    for kappa, tag in ((1e-6, "kappa = 1 mks (ранний угол memo)"),
                       (1e-3, "kappa = 1 ms"),
                       (11e-6, "нижняя граница бенда [2^40, 2^55] из постановки"),
                       (2e-3, "верхняя граница бенда [2^40, 2^55] из постановки")):
        w_star = ((math.pi / 4) * rm * kappa) ** 2
        ref.append({"kappa": tag, "log2 W*": round(math.log2(w_star), 1)})
    print("справка W*(kappa):", ref)

    out = {
        "predicate": "k in S routed-economy: k - канонический представитель "
                     "Lambda-орбиты (порядок lambda = 3, ровно 2 сравнения) "
                     "и сигнатура (старшие 2 бита) в интервале [0, 2)",
        "falsifier": {
            "pre-registered": "memo 10.11 «Дальше», коммит edd802d (до кода)",
            "branch O(N)": "канонизация сравнением со всеми элементами "
                           "орбиты вместо двух соседей",
            "outcome": "не сработал: в схеме канонизация - ровно 2 сравнения; "
                       "предикат 8*n*(n+m) + 64*(2m-1)*(n+m) + 20*n + 45 "
                       "Тоффоли (структурная формула, точна на всех точках "
                       "сетки), глубина O(n^2) на эталонной схеме (замер "
                       "886590 слоёв при n=256) и O(n log n) с префикс-суммами",
        },
        "construction": {
            "adder": "Cuccaro-структура: 2 Toffoli + 4 CNOT на бит, 1 ancilla; "
                     "выведена и бит-верифицирована исчерпывающе",
            "comparator": "заёмная цепочка, 5 Toffoli на бит",
            "mult": "широкий аккумулятор n+ceil(log2 n) бит, n inplace-сложений, "
                    "бинарное сокращение mod p за m шагов (каскад промежуточных "
                    "регистров неизбежен: реверс шага j требует R_{j+1})",
            "cleanup": "реверс под-схем в обратном хронологическом порядке; "
                       "каждая под-схема - биекция",
            "secp256k1 constants": "порядок N и lambda/lambda^2 из "
                                   "phase_signature_135bit.json (lambda^3 = 1 "
                                   "mod N проверена в рантайме)",
        },
        "verification": {
            "engine primitives": eng,
            "reduce bit-level": red_checks,
            "predicate bit-level": pred_checks,
            "note": "полный предикат и reduce на движке не прогоняются "
                    "(>45 проводов = 2^45 амплитуд); их корректность "
                    "покрывается бит-уровнем на всех 2^n входах",
        },
        "scaling": {
            "rows": rows,
            "fit toffoli": fit,
            "structural model": {
                "formula": "T(n) = 8*n*(n+m) + 64*(2m-1)*(n+m) + 20*n + 45, "
                           "m = ceil(log2 n)",
                "expanded": "8*n^2 + 136*n*log2(n) - 44*n + 128*log2(n)^2 "
                            "- 64*log2(n) + 45",
                "components": "8n(n+m): фаза умножений (2 умножения x n "
                              "сложений x 2(n+m) Тоффоли, x2 очистка); "
                              "64(2m-1)(n+m): каскад сокращения ((2m-1) "
                              "cond_sub x 16(n+m), x2 умножения x2 очистка); "
                              "20n: канонизация (4 компаратора x 5n); 45: "
                              "интервальная логика и AND",
                "exact match": f"{exact}/{len(rows)} точек сетки n=4..256, "
                               "0 расхождений",
                "component check": dec,
                "fit structural basis": fit_st,
                "n*log2(n) share at n=256": round(
                    136 * 8 * 256 / next(r["toffoli"] for r in rows
                                         if r["n"] == 256), 3),
                "asymptote": "член 136*n*log2(n) < 10% от 8n^2 только при "
                             "n >= 1871 - за пределами практической ширины; "
                             "для secp256k1 (n=256) рабочее число - прямой "
                             "замер 799277 = 12.2*n^2, не коэффициент фита",
                "kappa independence": "kappa = глубина x такт; глубина - "
                                      "прямой замер n=256 (886590 слоёв), "
                                      "префикс/древесная - от того же замера; "
                                      "фит Тоффоли в kappa и W* не входит, "
                                      "W* бенд не меняется",
            },
            "adder n=256 depth": d_add,
            "adder depth per bit": round(d_add_bit, 3),
            "brent-kung depth": round(bk, 1),
            "prefix factor": round(factor, 5),
        },
        "kappa_w_star": {
            "formula": "W* = ((pi/4)*rate*M*kappa)^2 - определение из memo 10.11; "
                       "rate/M - замеры canonical H rows @400k",
            "rate steps/s": rate,
            "M rows": M,
            "rows": log2_rows,
            "reference": ref,
            "note n": "n в схеме - ширина арифметики log2(N) = 256 для "
                      "secp256k1; ширина окна W входит только в число "
                      "итераций Гровера, не в стоимость оракула",
        },
        "marking_iterations": {
            "S density": "routed 1.167 кандидата на гигантский шаг "
                         "(h_curve_summary.json economy)",
            "S(W)": "1.167*W/M",
            "grover iterations": round(iter_mark, 1),
            "grover iterations note": "(pi/4)*sqrt(W/|S|) = (pi/4)*sqrt(M/1.167) "
                                      "- константа, не зависит от W; но Гровер "
                                      "с предикатом возвращает равномерного "
                                      "кандидата из S, а не ключ - как финишёр "
                                      "он неполон без QRAM или EC-оракула",
        },
    }
    out_path = Path("results/qip_predicate_kappa.json")
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
    print(f"[saved] {out_path}")


if __name__ == "__main__":
    main()
