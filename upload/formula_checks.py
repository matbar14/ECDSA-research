import random

random.seed(2026)

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
BETA = 0x7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE
LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
B = 7
W = BETA * BETA % P

results = {}
cnt = {"M": 0, "S": 0}


def check(name, ok):
    r = results.setdefault(name, [0, 0])
    r[0] += bool(ok)
    r[1] += 1


def inv(v):
    return pow(v % P, -1, P)


def aadd(p, q):
    if p is None:
        return q
    if q is None:
        return p
    (x1, y1), (x2, y2) = p, q
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    s = (3 * x1 * x1 * inv(2 * y1) if x1 == x2 else (y2 - y1) * inv(x2 - x1)) % P
    x3 = (s * s - x1 - x2) % P
    return x3, (s * (x1 - x3) - y1) % P


def amul(k, p=(GX, GY)):
    r = None
    k %= N
    while k:
        if k & 1:
            r = aadd(r, p)
        p = aadd(p, p)
        k >>= 1
    return r


def aneg(p):
    return p[0], -p[1] % P


def proj(p):
    return (0, 1, 0) if p is None else (p[0], p[1], 1)


def aff(pt):
    x, y, z = pt
    if z % P == 0:
        return None
    iz = inv(z)
    return x * iz % P, y * iz % P


def on_curve(pt, b=B):
    x, y, z = pt
    return (y * y * z - x ** 3 - b * z ** 3) % P == 0 and (x % P, y % P, z % P) != (0, 0, 0)


def same(p1, p2):
    return aff(p1) == aff(p2)


def rz(p):
    x, y, z = proj(p)
    l = random.randrange(1, P)
    return x * l % P, y * l % P, z * l % P


def M(a, c):
    cnt["M"] += 1
    return a * c % P


def S(a):
    cnt["S"] += 1
    return a * a % P


def prefix(p, q):
    x1, y1, z1 = p
    x2, y2, z2 = q
    t0 = M(x1, x2)
    t1 = M(y1, y2)
    t2 = M(z1, z2)
    t3 = (M(x1 + y1, x2 + y2) - t0 - t1) % P
    t4 = (M(x1 + z1, x2 + z2) - t0 - t2) % P
    t5 = (M(y1 + z1, y2 + z2) - t1 - t2) % P
    b3 = 3 * B
    return t3, (t1 - b3 * t2) % P, t5, b3 * t4 % P, 3 * t0 % P, (t1 + b3 * t2) % P


def add_rcb(p, q):
    a, b, c, d, e, f = prefix(p, q)
    return (M(a, b) - M(c, d)) % P, (M(e, d) + M(f, b)) % P, (M(c, f) + M(a, e)) % P


def add_5m(p, q):
    a, b, c, d, e, f = prefix(p, q)
    s = (f + d) % P
    r = (b - c) % P
    m0 = M(s, r)
    m1 = M(f, c)
    m2 = M(d, (e - r) % P)
    m3 = M((a - s) % P, b)
    m4 = M(a, e)
    return (m0 + m1 + m3) % P, (m0 + m1 + m2) % P, (m1 + m4) % P


def psi(pt, b):
    x, y, z = pt
    y2 = S(y)
    z2 = S(z)
    return M(x, (y2 + 3 * b * z2) % P), M(y, (y2 - 9 * b * z2) % P), M(z, (y2 - b * z2) % P)


def tpl(pt):
    x, y, z = psi(psi(pt, B), -27 * B)
    return 3 * x % P, y, 27 * z % P


def tau(pt):
    x, y, z = psi(pt, B)
    o = (1 - W) % P
    return o * x % P, y, pow(o, 3, P) * z % P


def dbl(pt):
    x, y, z = pt
    y2 = y * y % P
    z2 = z * z % P
    return (2 * x * y * (y2 - 9 * B * z2) % P,
            ((y2 + 9 * B * z2) ** 2 - 108 * B * B * z2 * z2) % P,
            8 * y2 * y * z % P)


def omega(p, s=1, j=1):
    x, y = p
    return x * pow(W, j, P) % P, y if s == 1 else -y % P


def test_curve_formulas():
    pts = [amul(random.randrange(1, N)) for _ in range(60)]
    cases = [(pts[i], pts[i + 1]) for i in range(40)]
    for i in range(10):
        cases += [(pts[i], pts[i]), (pts[i], aneg(pts[i])), (None, pts[i]), (pts[i], None)]
    cases.append((None, None))
    for p, q in cases:
        ref = proj(aadd(p, q))
        pp, qq = rz(p), rz(q)
        r1 = add_rcb(pp, qq)
        r2 = add_5m(pp, qq)
        check("add RCB 12M", same(r1, ref) and on_curve(r1))
        check("add 5M block 11M", same(r2, ref) and on_curve(r2))
        if q is not None:
            check("mixed add 10M", same(add_5m(pp, proj(q)), ref))
    for p in pts[:30] + [None]:
        pp = rz(p)
        three = None if p is None else aadd(aadd(p, p), p)
        two = None if p is None else aadd(p, p)
        check("psi lands on E_-27b", on_curve(psi(pp, B), -27 * B))
        check("complete tripling", same(tpl(pp), proj(three)))
        check("complete doubling", same(dbl(pp), proj(two)))
    for p in pts[:30]:
        check("tau = 1 - omega", same(tau(rz(p)), proj(aadd(p, aneg(omega(p))))))
        x = p[0]
        xt = aff(tau(proj(p)))[0]
        check("x(tau P) = -b^2(x^3+28)/(3x^2)", xt == (-W * W * (pow(x, 3, P) + 28) * inv(3 * x * x)) % P)
    check("omega P = lambda^2 P", all(omega(p) == amul(LAM * LAM, p) for p in pts[:3]))
    check("no F_p points with x = 0", pow(7, (P - 1) // 2, P) != 1)
    costs = []
    for f in (add_rcb, add_5m, tpl):
        cnt.update(M=0, S=0)
        if f is tpl:
            f(rz(pts[0]))
        else:
            f(rz(pts[0]), rz(pts[1]))
        costs.append(dict(cnt))
    return costs


def zm(p, q):
    x1, y1 = p
    x2, y2 = q
    return x1 * x2 - 3 * y1 * y2, x1 * y2 + x2 * y1 + 3 * y1 * y2


def znorm(p):
    return p[0] ** 2 + 3 * p[0] * p[1] + 3 * p[1] ** 2


def zunit(s, j):
    u = (1, 0)
    for _ in range(j):
        u = zm(u, (1, -1))
    return u if s == 1 else (-u[0], -u[1])


BASE = [(1, 0), (-2, 0), (-1, -1)]
DIG = {}
for ci, c in enumerate(BASE):
    for s in (1, -1):
        for j in range(3):
            d = zm(zunit(s, j), c)
            DIG[(d[0] % 9, d[1] % 3)] = (d, ci, s, j)


def zdiv(e, dv):
    conj = (dv[0] + 3 * dv[1], -dv[1])
    q = zm(e, conj)
    n = znorm(dv)
    return q[0] % n == 0 and q[1] % n == 0


def recode(a, b):
    ds = []
    while (a, b) != (0, 0):
        if a % 3:
            d = DIG[(a % 9, b % 3)]
            ds.append(d)
            a -= d[0][0]
            b -= d[0][1]
        else:
            ds.append(None)
        a, b = a + b, -a // 3
    return ds


def glv(k):
    a1, b1 = 0x3086D221A7D46BCDE86C90E49284EB15, -0xE4437ED6010E88286F547FA90ABFE4C3
    a2, b2 = 0x114CA50F7A8E2F3F657C1108D9D44CFD8, a1
    c1 = (b2 * k + N // 2) // N
    c2 = (-b1 * k + N // 2) // N
    k1 = k - c1 * a1 - c2 * a2
    k2 = -c1 * b1 - c2 * b2
    assert (k1 + k2 * LAM - k) % N == 0
    return k1, k2


def tnaf_mul(k, pt):
    k1, k2 = glv(k)
    ds = recode(k1 - 2 * k2, k2)
    tab = [pt, amul(N - 2, pt), aadd(amul(N - 1, pt), aneg(aff(tau(proj(pt)))))]
    ttab = [aff(tau(proj(q))) for q in tab]
    r = (0, 1, 0)
    n = len(ds)
    ntpl = nadd = 0
    for m in range((n - 1) // 2, -1, -1):
        r = tpl(r)
        ntpl += 1
        for dl in (1, 0):
            i = 2 * m + dl
            if i < n and ds[i]:
                _, ci, s, j = ds[i]
                q = omega((ttab if dl else tab)[ci], s * (-1) ** m, (j + m) % 3)
                r = add_5m(r, (q[0], q[1], 1))
                nadd += 1
    return r, n, ntpl, nadd


def test_tnaf():
    classes = {(x, y) for x in range(9) for y in range(3) if x % 3}
    check("digits cover 18 classes", set(DIG) == classes)
    check("digit norms 1, 4, 7", sorted({znorm(c) for c in BASE}) == [1, 4, 7])
    t3 = zm((0, 1), zm((0, 1), (0, 1)))
    check("tau^3 lemma", all(zdiv((x, y), t3) == (x % 9 == 0 and y % 3 == 0)
                             for x in range(-30, 30) for y in range(-30, 30)))
    ls, tp, ad = [], [], []
    for _ in range(300):
        k = random.randrange(1, N)
        r, n, nt, na = tnaf_mul(k, (GX, GY))
        check("tau-NAF kP", same(r, proj(amul(k))))
        ls.append(n)
        tp.append(nt)
        ad.append(na)
    for k in [1, 2, 3, N - 1, N - 2, LAM, N - LAM, (1 - LAM) % N]:
        check("tau-NAF special scalars", same(tnaf_mul(k, (GX, GY))[0], proj(amul(k))))
    avg = lambda v: sum(v) / len(v)
    return avg(ls), avg(tp), avg(ad)


Q = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
H = pow(2, -1, Q)


def f2add(a, b):
    return (a[0] + b[0]) % Q, (a[1] + b[1]) % Q


def f2sub(a, b):
    return (a[0] - b[0]) % Q, (a[1] - b[1]) % Q


def f2mul(a, b):
    return (a[0] * b[0] - a[1] * b[1]) % Q, (a[0] * b[1] + a[1] * b[0]) % Q


def f2hlv(a):
    return a[0] * H % Q, a[1] * H % Q


def f2conj(a):
    return a[0], -a[1] % Q


def f2k(n):
    return n % Q, 0


def f2rnd():
    return random.randrange(Q), random.randrange(Q)


def f2pow(a, e):
    r = (1, 0)
    while e:
        if e & 1:
            r = f2mul(r, a)
        a = f2mul(a, a)
        e >>= 1
    return r


def imul(w):
    return -w[1] % Q, w[0]


def isub(w1, w2):
    return (w2[1] - w1[1]) % Q, (w1[0] - w2[0]) % Q


def jmul(w):
    w0 = (w[0] + w[1]) * H % Q
    w1 = (w[0] - w[1]) * H % Q
    return (w1, w0), (w0, -w1 % Q)


def toom_eval(a):
    a0, a1, a2 = a
    e5 = f2add(a0, a2)
    return a0, f2add(e5, a1), f2sub(e5, a1), f2add(f2sub(a0, a2), imul(a1)), a2


def toom_interp(x, y):
    e0 = f2mul(x[0], y[0])
    e1 = f2hlv(f2mul(x[1], y[1]))
    e2 = f2hlv(f2mul(x[2], y[2]))
    e3 = f2hlv(f2mul(x[3], y[3]))
    e4 = f2mul(x[4], y[4])
    e5 = f2add(e0, e4)
    e6, e7 = jmul(e1)
    e8, e9 = jmul(e2)
    c2 = f2sub(f2add(e1, e2), e5)
    e5 = isub(e5, e3)
    return [e0, f2add(f2sub(e7, e8), e5), c2, f2sub(f2sub(e6, e9), e5), e4]


def school(a, b):
    c = [(0, 0)] * 5
    for i in range(3):
        for j in range(3):
            c[i + j] = f2add(c[i + j], f2mul(a[i], b[j]))
    return c


def test_toom():
    for _ in range(200):
        a = tuple(f2rnd() for _ in range(3))
        b = tuple(f2rnd() for _ in range(3))
        check("imaginary Toom-3 mul", toom_interp(toom_eval(a), toom_eval(b)) == school(a, b))
        check("imaginary Toom-3 sqr", toom_interp(toom_eval(a), toom_eval(a)) == school(a, a))


Z = next((0, c) for c in range(1, 50) if f2pow((0, c), (Q * Q - 1) // 3) != (1, 0))
G3 = f2pow(Z, (Q ** 3 - 1) // 3)
ONE6 = ((1, 0), (0, 0), (0, 0))


def f6mul(x, y):
    c = [(0, 0)] * 5
    for i, a in enumerate(x):
        for j, b in enumerate(y):
            c[i + j] = f2add(c[i + j], f2mul(a, b))
    return f2add(c[0], f2mul(Z, c[3])), f2add(c[1], f2mul(Z, c[4])), c[2]


def f6pow(x, e):
    r = ONE6
    while e:
        if e & 1:
            r = f6mul(r, x)
        x = f6mul(x, x)
        e >>= 1
    return r


def frob3(x):
    a, b, c = x
    return f2conj(a), f2mul(f2conj(b), G3), f2mul(f2conj(c), f2mul(G3, G3))


def cube_formula(x):
    a, b, c = x
    ab = f2conj(a)
    a2 = f2mul(a, a)
    t = f2sub(f2mul(f2k(3), a2), f2mul(f2k(2), ab))
    nb = f2mul(b, f2conj(b))
    nc = f2mul(c, f2conj(c))
    c0 = f2add(f2mul(f2k(9), f2mul(a, f2sub(a2, ab))), (1, 0))
    c1 = f2mul(f2k(3), f2sub(f2mul(b, t), f2mul(Z, nc)))
    c2 = f2mul(f2k(3), f2add(f2mul(c, t), nb))
    return c0, c1, c2


def test_cube():
    check("conjugation maps s to -s", G3 == (Q - 1, 0))
    for _ in range(10):
        x = tuple(f2rnd() for _ in range(3))
        a = f6mul(frob3(x), f6pow(x, Q ** 6 - 2))
        ab = (f2conj(a[0]), f2sub((0, 0), f2conj(a[1])), f2conj(a[2]))
        check("a^(p^3+1) = 1", f6mul(a, ab) == ONE6)
        check("cube, only a^(p^3+1) = 1", cube_formula(a) == f6pow(a, 3))
        a = f6pow(a, Q + 1)
        check("a^Phi6(p) = 1", f6pow(a, Q * Q - Q + 1) == ONE6)
        check("cube, a in Phi6(p) subgroup", cube_formula(a) == f6pow(a, 3))


costs = test_curve_formulas()
ell, ntpl, nadd = test_tnaf()
test_toom()
test_cube()

for name, (ok, total) in results.items():
    print(f"{'OK  ' if ok == total else 'FAIL'} {name}: {ok}/{total}")
print("cost RCB:", costs[0], "5M:", costs[1], "tpl:", costs[2])
print(f"tau-NAF: digits {ell:.1f}, tpl {ntpl:.1f}, add {nadd:.1f}, density 1/{ell / nadd:.2f}, "
      f"model {10 * ntpl + 10 * nadd + 170:.0f}M (paper 1434M)")
