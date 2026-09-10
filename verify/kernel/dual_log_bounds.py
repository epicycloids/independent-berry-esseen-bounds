from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
ROUND27 = HERE.parent
for path in (ROUND27.parent / '26', ROUND27, ROUND27 / 'logarithm', ROUND27 / 'entropy_moments', ROUND27 / 'higher_variance'):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
import stable_bounds
from interval_bounds import K, KHI, al, au, up, down, sqpos, cubepos, m_upper, exp_upper_from_nonnegative_lower
from higher_moments import higher_moment_range
ORIGINAL = stable_bounds.log_range_cells
SUPPORT_GRID = 1024
LOG_GRID = 4096
LAMBDAS = (0.125, 0.1875, 0.25, 0.375, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0)

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

@lru_cache(maxsize=1)
def logarithm_upper_table():
    last = int(np.ceil(au(1 / (216 * K * K)) * LOG_GRID))
    values, derivatives = (np.empty(last + 1), np.empty(last + 1))
    for j in range(last + 1):
        left, right = (arb(j) / LOG_GRID, arb(j + 1) / LOG_GRID)
        assert right < arb(1) / 2
        values[j] = max(0.0, au(g2_arb(left)))
        derivatives[j] = au(2 * right / (1 - 2 * right))
    values[0] = 0.0
    values.flags.writeable = derivatives.flags.writeable = False
    return (values, derivatives)

def g2_upper(m):
    m = np.asarray(m, dtype=float)
    values, derivatives = logarithm_upper_table()
    index = np.floor(m * LOG_GRID).astype(np.intp)
    assert np.all((index >= 0) & (index < len(values)))
    left = index / LOG_GRID
    difference = np.maximum(0.0, up(m - left))
    return np.maximum(0.0, up(values[index] + up(derivatives[index] * difference)))

def phi_upper(m, row, log_upper=None):
    price, _, _, _, _, qstar_hi, conjugate_lo = row
    if log_upper is None:
        log_upper = g2_upper(m)
    affine = np.maximum(0.0, up(up(price * m) - conjugate_lo))
    return np.where(m >= qstar_hi, np.minimum(log_upper, affine), log_upper)

def dual_exponent_lower(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False, prices=None, return_details=False):
    tlo, thi = np.broadcast_arrays(np.asarray(tlo, dtype=float), np.asarray(thi, dtype=float))
    assert np.all((0 <= tlo) & (tlo <= thi))
    zero = np.zeros_like(tlo)
    if Vlo <= 0 or Vhi <= 0 or dhi <= 0:
        details = dict(price=zero, support=zero, enabled=np.zeros_like(tlo, bool))
        return (zero, details) if return_details else zero
    tau_lo, tau_hi = (max(0.0, taulo), min(Bhi, max(0.0, tauhi)))
    assert tau_lo <= tau_hi
    t2lo, t2hi = sqpos((tlo, thi))
    t3lo, t3hi = cubepos((tlo, thi))
    t2lo, t3lo = (np.maximum(0.0, t2lo), np.maximum(0.0, t3lo))
    t2hi, t3hi = (np.maximum(0.0, t2hi), np.maximum(0.0, t3hi))
    khi = up(KHI * t3hi)
    sigma = np.maximum(0.0, down(down(Vlo / 2) * t2lo))
    m = m_upper(thi, dhi) if deleted else zero
    d32hi = au(arb(dhi) * arb(dhi).sqrt()) if deleted else 0.0
    argument = up(up(au(4 * K) * thi) * au(arb(dhi).sqrt()))
    best = np.full_like(tlo, -np.inf)
    chosen_price, chosen_u = (zero.copy(), zero.copy())
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        guide_tau = max(0.0, tau_lo - dhi ** 1.5) if deleted else tau_lo
        guide_V = max(0.0, Vhi - dhi) if deleted else Vhi
        mean = np.divide(4 * float(K) * tlo ** 3 * guide_tau, thi ** 2 * guide_V, out=zero.copy(), where=thi > 0)
        mean = np.where(np.isfinite(mean), mean, 0.0)
    for row in support_table():
        price, u, aa, bb, domain, _, _ = row
        if prices is not None and price not in prices:
            continue
        wanted = np.searchsorted(u, mean, side='right') - 1
        permitted = np.searchsorted(-domain, -argument, side='right') - 1
        index = np.minimum(np.clip(wanted, 0, len(u) - 1), np.maximum(0, permitted))
        valid = argument <= domain[index]
        positive = np.maximum(0.0, down(aa[index] * t3lo))
        negative = up(bb[index] * t2hi)
        coefficient = down(positive + down((price - 1.0) * khi))
        chosen_tau = np.where(coefficient >= 0, tau_lo, tau_hi)
        shared = down(coefficient * chosen_tau)
        constant = up(up(up((1 + price) * khi) * Bhi) + up(negative * Vhi))
        delta = np.maximum(0.0, up(up(positive * d32hi) - down(negative * dhi))) if deleted else zero
        loss = up(m + delta) if deleted else zero
        value = down(down(down(sigma + shared) - constant) - loss)
        value = np.where(valid, value, -np.inf)
        improved = value > best
        best = np.maximum(best, value)
        chosen_price = np.where(improved, price, chosen_price)
        chosen_u = np.where(improved, u[index], chosen_u)
    best = np.maximum(0.0, best)
    if return_details:
        return (best, dict(price=chosen_price, support=chosen_u, enabled=best > 0))
    return best

def dual_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False, baseline='higher4', prices=None):
    exponent = dual_exponent_lower(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted, prices=prices)
    result = exp_upper_from_nonnegative_lower(exponent)
    if baseline == 'none':
        return result
    if baseline == 'higher4':
        old = higher_moment_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted, p=4)
    else:
        assert baseline == 'original'
        old = ORIGINAL(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted)
    return np.minimum(old, result)

@contextmanager
def using(baseline='higher4', prices=None):
    previous = stable_bounds.log_range_cells
    stable_bounds.log_range_cells = lambda *a, **k: dual_range(*a, **k, baseline=baseline, prices=prices)
    try:
        yield
    finally:
        stable_bounds.log_range_cells = previous

def warmup():
    return dict(support_rows=sum((len(row[1]) for row in support_table())), logarithm_intervals=len(logarithm_upper_table()[0]))
