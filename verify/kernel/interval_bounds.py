import math
import time
import numpy as np
from flint import arb, ctx
ctx.prec = 100
INF = float('inf')
U = 2.0 ** (-53)

def up(x):
    return np.nextafter(x, INF)

def down(x):
    return np.nextafter(x, -INF)

def au(x):
    return float(up(float(arb(x).upper())))

def al(x):
    return float(down(float(arb(x).lower())))
PI = arb.pi()
K = arb(49581) / 500000
KLO, KHI = (al(K), au(K))
TWO_PI_LO = al(2 * PI)
EXPCUT = au(arb(-25).exp())
FACT_LO = [al(arb(1) / math.factorial(j)) for j in range(17)]
MPEAK = au(1 / (216 * K * K))
XPEAK_LO = al(1 / (6 * K))

def ia(x):
    return (x, x)

def add(a, b):
    return (down(a[0] + b[0]), up(a[1] + b[1]))

def neg(a):
    return (-a[1], -a[0])

def sub(a, b):
    return add(a, neg(b))

def mul(a, b):
    vals = [a[i] * b[j] for i in (0, 1) for j in (0, 1)]
    return (down(np.minimum.reduce(vals)), up(np.maximum.reduce(vals)))

def divpos(a, b):
    return (down(a[0] / b[1]), up(a[1] / b[0]))

def sqpos(a):
    return (down(a[0] * a[0]), up(a[1] * a[1]))

def cubepos(a):
    return mul(sqpos(a), a)

def exp_upper_from_nonnegative_lower(p):
    p = np.maximum(0, p)
    x = np.minimum(p, 25.0) / 64.0
    q = np.full_like(x, FACT_LO[16])
    for j in range(15, -1, -1):
        q = down(down(q * x) + FACT_LO[j])
    y = up(1 / q)
    for _ in range(6):
        y = up(y * y)
    return np.minimum(1.0, np.where(p >= 25.0, EXPCUT, y))

def one_minus_cos_lower(z):
    zz = sqpos(z)
    co = [arb((-1) ** (j + 1)) / math.factorial(2 * j) for j in range(1, 7)]
    p = (al(co[-1]), au(co[-1]))
    for c in co[-2::-1]:
        p = add(mul(p, zz), (al(c), au(c)))
    return np.maximum(0, mul(p, zz)[0])

def psi_lower(t, V, E):
    t2 = sqpos(t)
    t3 = mul(t2, t)
    cubic = sub(mul(ia(V / 2), t2), mul((KLO, KHI), mul(ia(E), t3)))[0]
    if E <= 0:
        return np.maximum(0, mul(ia(V / 2), t2)[0])
    if V <= 0:
        return np.zeros_like(t[0])
    arg = divpos(mul(ia(E), t), ia(V))
    valid = (arg[0] >= 4.0) & (arg[1] <= TWO_PI_LO)
    z = (np.maximum(0, down(TWO_PI_LO - arg[1])), np.minimum(2.3, up(au(2 * PI) - arg[0])))
    z = (np.minimum(z[0], 2.3), np.clip(z[1], 0.0, 2.3))
    coslow = one_minus_cos_lower(z)
    pref = divpos(cubepos(ia(V)), sqpos(ia(E)))[0]
    cosine = down(pref * coslow)
    return np.maximum(0, np.where(valid, np.maximum(cubic, cosine), cubic))

def f_upper(t, V, E, single=False):
    p = psi_lower(t, V, E)
    if single:
        return np.minimum(1.0, up(np.sqrt(np.maximum(0, up(1 - down(2 * p))))))
    return exp_upper_from_nonnegative_lower(p)

def f_cells(t, V, E, single=False):
    f = f_upper(t, V, E, single)
    return np.maximum(f[:-1], f[1:])

def m_upper(t_hi, d):
    if d <= 0:
        return np.zeros_like(t_hi)
    sq = arb(d).sqrt()
    x = mul(ia(t_hi), (al(sq), au(sq)))
    x2 = sqpos(x)
    x3 = mul(x2, x)
    q = sub(mul(ia(0.5), x2), mul((al(2 * K), au(2 * K)), x3))[1]
    return np.where(x[1] <= XPEAK_LO, np.maximum(0, q), MPEAK)

def exp_cubic_cells(tlo, thi, V, E, d):
    positive = mul((KLO, KHI), mul(ia(E), cubepos(ia(thi))))[1]
    negative = mul(ia(V / 2), sqpos(ia(tlo)))[0]
    exponent = up(up(positive + m_upper(thi, d)) - negative)
    return exp_upper_from_nonnegative_lower(np.maximum(0, -exponent))

def log_product_cells(tlo, thi, Vlo, Vhi, Bhi, dhi, tauhi, deleted=False):
    if Vlo <= 0:
        return np.ones_like(tlo)
    t2lo = sqpos(ia(tlo))[0]
    t3lo = cubepos(ia(tlo))[0]
    t3hi = cubepos(ia(thi))[1]
    k = up(KHI * t3hi)
    sqrt_d = au(arb(dhi).sqrt())
    sqrt_V = au(arb(Vhi).sqrt())
    A = down(down(t2lo / 2) - up(up(2 * k) * sqrt_d))
    a = down(down(np.maximum(0, A) / sqrt_V) + down(KLO * t3lo))
    c = up(k * Bhi)
    m = m_upper(thi, dhi) if deleted else np.zeros_like(thi)
    valid = (A >= 0) & (a > 0)
    safe_a = np.where(valid, a, 1.0)
    x = (down(down(k / 2) / safe_a), up(up(k / 2) / safe_a))
    lo = down(down(c + np.maximum(m, x[0])) / safe_a)
    hi = up(up(c + np.maximum(m, x[1])) / safe_a)
    tl = np.maximum(0, np.minimum(tauhi, lo))
    th = np.maximum(0, np.minimum(tauhi, hi))
    q = np.maximum(0, down(down(safe_a * tl) - c))
    penalty = np.maximum(0, down(down(q * q) - up(m * m)))
    scalar = up(up(k * th) - penalty)
    scalar = np.where(valid, scalar, up(k * tauhi))
    negative = down(Vlo / 2 * t2lo)
    exponent = up(up(up(c + m) + scalar) - negative)
    return exp_upper_from_nonnegative_lower(np.maximum(0, -exponent))

def positive_dot(a, b):
    n = len(a)
    return float(up(np.dot(a, b) * (1 + 4 * n * U)))

def regularized_kernel(u):
    a = PI * u * (1 - u)
    b = (1 - u) * (PI * u).cos() / u.sinc_pi() + u
    return (a.abs_upper() ** 2 + b.abs_upper() ** 2).sqrt()

def two_kernel(u):
    if u.lower() >= arb(1) / 2:
        z = 1 - u
        imag = (1 - (PI * z).cos() / z.sinc_pi()) / PI
        return (z.abs_upper() ** 2 + imag.abs_upper() ** 2).sqrt()
    return regularized_kernel(u) / (PI * u)

class Weights:

    def __init__(self, Lmid, N=1024, log_correction=False):
        self.N = N
        self.log_correction = log_correction
        self.T = float(round(math.pi / Lmid, 8))
        self.s = float(min(0.14 + 0.52 * Lmid, 5 / self.T))
        T, s = (arb(self.T), arb(self.s))
        ts = [T * s * j / N for j in range(N + 1)]
        self.t = (np.array([al(x) for x in ts]), np.array([au(x) for x in ts]))
        self.t[0][0] = self.t[1][0] = 0.0
        self.tlo = self.t[0][:-1]
        self.thi = self.t[1][1:]
        self.gauss = np.array([au((-x * x / 2).exp()) for x in ts[:-1]])
        primitive = [(x * (x * x / 2).exp() - (PI / 2).sqrt() * (x / arb(2).sqrt()).erfi()) / 2 for x in ts]
        self.inner = np.array([au(primitive[j + 1] - primitive[j]) for j in range(N)])
        assert np.all(self.inner > 0)
        self.low = np.zeros(N)
        for j in range(1, N):
            ui = (s * j / N).union(s * (j + 1) / N)
            self.low[j] = au(regularized_kernel(ui) / PI * (arb(j + 1) / j).log())
        self.first = au(regularized_kernel(arb(0).union(s / N)) * ts[1] ** 3 / (18 * PI))
        us = [s + (1 - s) * j / N for j in range(N + 1)]
        high_ts = [T * u for u in us]
        self.ht = (np.array([al(x) for x in high_ts]), np.array([au(x) for x in high_ts]))
        self.high = np.array([au(two_kernel(us[j].union(us[j + 1])) * (us[j + 1] - us[j])) for j in range(N)])
        z = T * s / arb(2).sqrt()
        e = (-z * z).exp()
        er = z.erf()
        i0 = (PI / 2).sqrt() / T * er
        i1 = (1 - e) / (T * T)
        i2 = (PI / 2).sqrt() / T ** 3 * er - s * e / (T * T)
        self.G = au(i0 - i1 + PI * PI / 18 * i2 + (z * z).expint(1) / (2 * PI))
        assert all((np.isfinite(x).all() for x in [self.gauss, self.inner, self.low, self.high]))
        assert math.isfinite(self.first) and math.isfinite(self.G)

    def integrate(self, B, fcf, Llo, lowcf=None):
        cumul = up(np.cumsum(up(B * self.inner)) * (1 + 4 * self.N * U))
        delta = up(self.gauss * cumul)
        if lowcf is not None:
            delta = np.minimum(delta, up(up(lowcf + self.gauss) / Llo))
        delta[0] = 0.0
        i1 = up(self.first + positive_dot(self.low, delta))
        i2 = positive_dot(self.high, fcf)
        return float(up(i1 + up(up(i2 + self.G) / Llo)))

    def uniform(self, Llo, Lhi):
        d = au(arb(Lhi) ** (arb(2) / 3))
        E = au(2 * arb(Lhi))
        B = exp_cubic_cells(self.tlo, self.thi, 1.0, E, d)
        F = f_cells(self.ht, 1.0, E)
        lowcf = f_cells(self.t, 1.0, E)
        return self.integrate(B, F, Llo, lowcf)

    def box(self, Llo, Lhi, dlo, dhi, blo, bhi):
        dl = arb(dlo)
        dh = arb(dhi)
        if bhi < al(dl * dl.sqrt()) or blo > Lhi:
            return None
        blo = max(blo, al(dl * dl.sqrt()))
        bhi = min(bhi, Lhi)
        if blo > bhi:
            return None
        whi = au(1 - dl)
        wlo = max(0, al(1 - dh))
        dh3 = dh * dh.sqrt()
        dl3 = dl * dl.sqrt()
        tau = min(au(arb(Lhi) - arb(blo) + dh3), au(dh.sqrt()))
        E = au(arb(Lhi) + arb(tau))
        ER = max(0, au(arb(Lhi) - arb(blo) + arb(tau) - dl3))
        Ej = au(arb(bhi) + dh3)
        d2 = min(dhi, whi, au((arb(Lhi) - arb(blo)) ** (arb(2) / 3)))
        lamlo = max(0, al(arb(blo) / arb(Lhi)))
        lamhi = min(1, au(arb(bhi) / arb(Llo)))
        R = f_cells(self.t, wlo, ER)
        A = f_cells(self.t, dlo, Ej, True)
        D = exp_cubic_cells(self.tlo, self.thi, wlo, ER, d2)
        tauR = max(0, au(arb(tau) - dl3))
        BR = max(0, au(arb(Lhi) - arb(blo)))
        if self.log_correction:
            R = np.minimum(R, log_product_cells(self.tlo, self.thi, wlo, whi, BR, d2, tauR))
            D = np.minimum(D, log_product_cells(self.tlo, self.thi, wlo, whi, BR, d2, tauR, True))
        other = np.minimum(1.0, up(A * D))

        def mixture(lam):
            return up(up(lam * R) + up(up(1 - lam) * other))
        B = np.minimum(1.0, np.maximum(mixture(lamlo), mixture(lamhi)))
        B = np.minimum(B, exp_cubic_cells(self.tlo, self.thi, 1.0, E, dhi))
        if self.log_correction:
            B = np.minimum(B, log_product_cells(self.tlo, self.thi, 1.0, 1.0, Lhi, dhi, tau, True))
        F = np.minimum(f_cells(self.ht, 1.0, E), up(f_cells(self.ht, dlo, Ej, True) * f_cells(self.ht, wlo, ER)))
        lowcf = np.minimum(f_cells(self.t, 1.0, E), up(A * R))
        if self.log_correction:
            hl, hh = (self.ht[0][:-1], self.ht[1][1:])
            Rhi = np.minimum(f_cells(self.ht, wlo, ER), log_product_cells(hl, hh, wlo, whi, BR, d2, tauR))
            F = np.minimum(F, np.minimum(up(f_cells(self.ht, dlo, Ej, True) * Rhi), log_product_cells(hl, hh, 1.0, 1.0, Lhi, dhi, tau)))
            lowcf = np.minimum(lowcf, log_product_cells(self.tlo, self.thi, 1.0, 1.0, Lhi, dhi, tau))
        return self.integrate(B, F, Llo, lowcf)

def elementary_checks():
    rng = np.random.default_rng(260906)
    p = np.r_[0.0, np.geomspace(1e-20, 100, 300), rng.uniform(0, 40, 100)]
    e = exp_upper_from_nonnegative_lower(p)
    for x, y in zip(p, e):
        assert arb(float(y)) >= (-arb(float(x))).exp(), (x, y)
    z = rng.uniform(0, 2.3, 300)
    lower = one_minus_cos_lower((z, z))
    for x, y in zip(z, lower):
        assert arb(float(y)) <= 1 - arb(float(x)).cos(), (x, y)
    for _ in range(80):
        V = float(rng.uniform(0.001, 1))
        E = float(rng.uniform(2 * V ** 1.5, 2)) if V < 1 else 2.0
        t = np.r_[rng.uniform(0.001, 15, 16), 4 * V / E, 2 * math.pi * V / E]
        for single in [False, True]:
            out = f_upper((t, t), V, E, single)
            for x, y in zip(t, out):
                aa = arb(E) * arb(float(x)) / arb(V)
                if aa < 4:
                    psi = arb(V) * arb(float(x)) ** 2 / 2 - K * arb(E) * arb(float(x)) ** 3
                elif aa < 2 * PI:
                    psi = arb(V) ** 3 / arb(E) ** 2 * (1 - aa.cos())
                else:
                    psi = arb(0)
                if psi < 0:
                    psi = arb(0)
                ref = (1 - 2 * psi).sqrt() if single else (-psi).exp()
                assert y == 1.0 or arb(float(y)) >= ref, str((V, E, x, single, y, ref))
    assert np.finfo(np.float64).bits == 64 and np.nextafter(0.0, 1.0) > 0
    log_count = 0
    for _ in range(120):
        V = float(rng.uniform(0.05, 1))
        d = float(rng.uniform(0.001, V))
        B = float(rng.uniform(d ** 1.5, 1.7))
        tc = float(min(B, math.sqrt(d) * V))
        tl = float(rng.uniform(0.001, 5))
        th = tl + float(rng.uniform(0, 0.02))
        for deleted in [False, True]:
            y = float(log_product_cells(np.array([tl]), np.array([th]), V, V, B, d, tc, deleted)[0])
            for tval in [tl, (tl + th) / 2, th]:
                t = arb(tval)
                vv = arb(V)
                bb = arb(B)
                dd = arb(d)
                cap = arb(tc)
                mm = min(t * dd.sqrt(), 1 / (6 * K))
                mm = mm * mm / 2 - 2 * K * mm ** 3 if deleted else arb(0)
                AA = t * t / 2 - 2 * K * dd.sqrt() * t ** 3
                kk = K * t ** 3
                if AA >= 0:
                    aa = AA / vv.sqrt() + kk
                    cc = kk * bb
                    tt = min(cap, (cc + max(mm, kk / (2 * aa))) / aa)
                    qq = max(arb(0), aa * tt - cc) ** 2
                    penalty = max(arb(0), qq - mm * mm)
                else:
                    tt = cap
                    penalty = arb(0)
                ref = (-vv * t * t / 2 + kk * bb + mm + kk * tt - penalty).exp()
                assert y == 1.0 or arb(y) >= ref, str((V, B, d, tc, tl, th, deleted, y, ref))
                log_count += 1
    return {'exp_enclosures': len(p), 'cosine_enclosures': len(z), 'cf_enclosures': 80 * 18 * 2, 'log_enclosures': log_count, 'all_passed': True, 'arb_bits': ctx.prec, 'numpy': np.__version__}
