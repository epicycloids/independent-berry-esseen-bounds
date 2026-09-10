from functools import lru_cache
import numpy as np
from flint import arb
from flint import ctx
import math
ctx.prec = 100
K = arb(49581)/500000
def al(x): return math.nextafter(float(arb(x).lower()), -math.inf)
def au(x): return math.nextafter(float(arb(x).upper()), math.inf)
SUPPORT_GRID = 1024
LAMBDAS = (0.75, 1.0, 1.5, 2.0, 3.0)

def g2_arb(q):
    return -(1 - 2 * q).log() / 2 - q

def phi_arb(r, lam):
    lam = arb(lam)
    qstar = lam / (2 * (1 + lam))
    if r <= qstar:
        return g2_arb(r)
    return lam * r - (lam - (1 + lam).log()) / 2

@lru_cache(maxsize=1)
def support_table():
    a = 1 / (32 * K * K)
    assert 0 < a < arb(27) / 8
    uf = np.arange(1, 548, dtype=float) / SUPPORT_GRID
    rf = float(a) * uf * uf * (1 - uf)
    gf = -np.log1p(-2 * rf) / 2 - rf
    gpf = 2 * rf / (1 - 2 * rf)
    af = (gpf * float(a) * uf * uf * (2 - 3 * uf) - 2 * gf) / uf ** 3
    bf = af * uf - gf / uf ** 2
    rows = []
    for price in LAMBDAS:
        qf = price / (2 * (1 + price))
        ids = np.flatnonzero(rf < qf)
        u = uf[ids]
        alpha, beta = (af[ids], bf[ids])
        lo, hi = (u.copy(), np.ones_like(u))
        conjugate = (price - np.log1p(price)) / 2
        for _ in range(44):
            middle = (lo + hi) / 2
            r = float(a) * middle ** 2 * (1 - middle)
            g = -np.log1p(-2 * r) / 2 - r
            phi = np.where(r <= qf, g, price * r - conjugate)
            residual = phi / middle ** 2 - alpha * middle + beta
            lo = np.where(residual >= 0, middle, lo)
            hi = np.where(residual >= 0, hi, middle)
        candidates = np.maximum(u, lo * (1 - 2.0 ** (-24)))
        us, aa, bb, domains = ([0.0], [0.0], [0.0], [1.0])
        lam = arb(price)
        qstar = lam / (2 * (1 + lam))
        for j, value in enumerate(u):
            v = arb(float(value))
            r = a * v * v * (1 - v)
            g = g2_arb(r)
            gp = 2 * r / (1 - 2 * r)
            gpp = 2 / (1 - 2 * r) ** 2
            zrp = a * v * v * (2 - 3 * v)
            D = gpp * zrp ** 2 - 6 * r * gp + 6 * g
            if not (r <= qstar and D >= 0):
                continue
            alpha = (gp * zrp - 2 * g) / v ** 3
            beta = alpha * v - g / v ** 2
            assert alpha >= 0 and beta >= 0
            endpoint = arb(float(candidates[j]))
            assert v <= endpoint <= 1
            re = a * endpoint ** 2 * (1 - endpoint)
            assert phi_arb(re, lam) / endpoint ** 2 >= alpha * endpoint - beta
            us.append(value)
            aa.append(max(0.0, al(64 * K ** 3 * alpha)))
            bb.append(max(0.0, au(16 * K * K * beta)))
            domains.append(float(candidates[j]))
        us, aa, bb, domains = map(np.asarray, (us, aa, bb, domains))
        assert np.all(domains[1:] <= domains[:-1])
        assert np.all(us[1:] > us[:-1])
        for x in (us, aa, bb, domains):
            x.flags.writeable = False
        rows.append((price, us, aa, bb, domains, au(qstar), max(0.0, al((lam - (1 + lam).log()) / 2))))
    return tuple(rows)
