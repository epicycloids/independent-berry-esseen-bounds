from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
from flint import arb
_round26 = Path(__file__).resolve().parent.parent.parent / '26'
if _round26.is_dir() and str(_round26) not in sys.path:
    sys.path.insert(0, str(_round26))
import stable_bounds
from interval_bounds import K, KLO, KHI, al, au, up, down, ia, sqpos, cubepos, m_upper, exp_upper_from_nonnegative_lower
ORIGINAL = stable_bounds.log_range_cells
GRID = 512
MAX_INDEX = 504
R_CAP = MAX_INDEX / GRID
METHOD = 'tangent'
PASSES = 2

@lru_cache(maxsize=1)
def tangent_table():
    alpha = np.empty(MAX_INDEX + 1)
    beta = np.empty(MAX_INDEX + 1)
    alpha[0], beta[0] = (1.0, 0.0)
    for j in range(1, MAX_INDEX + 1):
        r = arb(j) / GRID
        phi = 2 * (-(1 - r).log() - r) / (r * r)
        alpha[j] = max(1.0, al(2 / (1 - r) - phi))
        beta[j] = max(0.0, au(r * (1 / (1 - r) - phi)))
    assert np.all(np.isfinite(alpha)) and np.all(np.isfinite(beta))
    alpha.flags.writeable = beta.flags.writeable = False
    return (alpha, beta)

def _components(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted):
    tlo, thi = np.broadcast_arrays(np.asarray(tlo, dtype=float), np.asarray(thi, dtype=float))
    tl, th = (max(0.0, min(Bhi, taulo)), max(0.0, min(Bhi, tauhi)))
    assert 0 <= tl <= th and 0 < Vlo <= Vhi and (Bhi >= 0) and (dhi >= 0)
    t2lo = np.maximum(0.0, sqpos(ia(tlo))[0])
    t2hi = sqpos(ia(thi))[1]
    t3lo = np.maximum(0.0, cubepos(ia(tlo))[0])
    t3hi = cubepos(ia(thi))[1]
    klo = np.maximum(0.0, down(KLO * t3lo))
    khi = up(KHI * t3hi)
    sigma_lo = np.maximum(0.0, down(down(Vlo / 2) * t2lo))
    sigma_hi = up(up(Vhi / 2) * t2hi)
    A = down(down(t2lo / 2) - up(up(2 * khi) * au(arb(dhi).sqrt())))
    a = down(down(np.maximum(0.0, A) / au(arb(Vhi).sqrt())) + klo)
    c = up(khi * Bhi)
    m = m_upper(thi, dhi) if deleted else np.zeros_like(thi)
    valid = (A >= 0) & (a > 0)
    a = np.where(valid, a, 1.0)
    base = up(up(c + m) - sigma_lo)
    return (a, c, klo, khi, sigma_hi, m, valid, base, tl, th)

def _quadratic_exponent(a, c, k, m, valid, base, tl, th):
    xlo = down(down(k / 2) / a)
    xhi = up(up(k / 2) / a)
    arglo = down(down(c + np.maximum(m, xlo)) / a)
    arghi = up(up(c + np.maximum(m, xhi)) / a)
    low = np.maximum(tl, np.minimum(th, arglo))
    high = np.maximum(tl, np.minimum(th, arghi))
    q = np.maximum(0.0, down(down(a * low) - c))
    charge = np.maximum(0.0, down(down(q * q) - up(m * m)))
    scalar = np.where(valid, up(up(k * high) - charge), up(k * th))
    return (up(base + scalar), (low + high) / 2)

def _quadratic_max(a, c, w, alpha, tl, th):
    positive = w > 0
    wp = np.maximum(0.0, w)
    denlo = np.maximum(0.0, down(down(2 * alpha) * a))
    denhi = up(up(2 * alpha) * a)
    xlo = down(wp / denhi)
    xhi = up(np.divide(wp, denlo, out=np.full_like(wp, np.inf), where=denlo > 0))
    arglo = down(down(c + xlo) / a)
    arghi = up(up(c + xhi) / a)
    low = np.where(positive, np.maximum(tl, np.minimum(th, arglo)), tl)
    high = np.where(positive, np.maximum(tl, np.minimum(th, arghi)), tl)
    q = np.maximum(0.0, down(down(a * low) - c))
    charge = np.maximum(0.0, down(alpha * np.maximum(0.0, down(q * q))))
    value = up(up(w * high) - charge)
    return (value, (low + high) / 2)

def _choose_tangent(a, c, klo, sigma_hi, m, guide):
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        q = np.maximum(0.0, a * guide - c)
        Q = np.maximum(0.0, q * q - m * m)
        S = sigma_hi - 2 * klo * guide - m
        r = np.divide(2 * Q, S, out=np.zeros_like(Q), where=S > 0)
        r = np.nan_to_num(r, nan=0.0, posinf=R_CAP, neginf=0.0)
    index = np.rint(GRID * np.clip(r, 0.0, R_CAP)).astype(np.intp)
    alpha, beta = tangent_table()
    return (alpha[index], beta[index])

def _tangent_exponent(parts, original, guide, passes):
    a, c, klo, khi, sigma_hi, m, valid, base, tl, th = parts
    best = original
    for _ in range(passes):
        alpha, beta = _choose_tangent(a, c, klo, sigma_hi, m, guide)
        w = up(khi - down(down(2 * beta) * klo))
        scalar, guide = _quadratic_max(a, c, w, alpha, tl, th)
        loss = np.maximum(0.0, up(up(alpha * up(m * m)) - down(beta * m)))
        offset = up(up(base + loss) + up(beta * sigma_hi))
        candidate = up(offset + scalar)
        best = np.minimum(best, np.where(valid, candidate, original))
    return best

def _phi_lower(r, terms=12):
    assert isinstance(terms, int) and terms >= 1
    r = np.asarray(r, dtype=float)
    assert np.all((r >= 0) & (r < 1))
    d1 = up(1 + up(np.sqrt(up(1 - r))))
    r1 = np.maximum(0.0, down(r / d1))
    d2 = up(1 + up(np.sqrt(up(1 - r1))))
    r2 = np.maximum(0.0, down(r1 / d2))
    den = up(2 - r2)
    v = np.maximum(0.0, down(r2 / den))
    v2 = np.maximum(0.0, down(v * v))
    coefficients = _atanh_coefficients(terms)
    p = np.full_like(r, coefficients[-1])
    for coefficient in coefficients[-2::-1]:
        p = np.maximum(0.0, down(down(p * v2) + coefficient))
    den3 = up(up(den * den) * den)
    psi = down(down(1 / den) + down(down(down(2 * r2) / den3) * p))
    psi1 = down(down(down(2 * psi) + 1) / up(d2 * d2))
    psi0 = down(down(down(2 * psi1) + 1) / up(d1 * d1))
    return np.maximum(1.0, down(2 * psi0))

@lru_cache(maxsize=16)
def _atanh_coefficients(terms):
    return tuple((al(arb(1) / (2 * j + 3)) for j in range(terms)))

def weighted_charge_lower(Qlo, Shi, terms=12):
    Qlo, Shi = np.broadcast_arrays(np.asarray(Qlo, dtype=float), np.asarray(Shi, dtype=float))
    valid = (Qlo > 0) & (Shi > 0)
    den = np.where(valid, Shi, 1.0)
    Q = np.where(valid, Qlo, 0.0)
    r = np.maximum(0.0, down(down(2 * Q) / den))
    r = np.minimum(r, R_CAP)
    return np.maximum(0.0, down(Q * _phi_lower(r, terms)))

def _direct_exponent(parts, original):
    a, c, klo, khi, sigma_hi, m, valid, base, tl, th = parts
    q = np.maximum(0.0, down(down(a * tl) - c))
    Q = np.maximum(0.0, down(down(q * q) - up(m * m)))
    Q = np.where(valid, Q, 0.0)
    S = up(sigma_hi - down(down(2 * klo) * tl))
    charge = weighted_charge_lower(Q, S)
    candidate = up(up(base + up(khi * th)) - charge)
    return np.minimum(original, np.where(valid, candidate, original))

def weighted_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False, method=None, passes=None):
    if Vlo <= 0:
        return np.ones_like(np.asarray(tlo, dtype=float))
    method = METHOD if method is None else method
    passes = PASSES if passes is None else passes
    assert method in ('tangent', 'direct', 'both', 'quadratic')
    assert isinstance(passes, int) and 1 <= passes <= 4
    parts = _components(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted)
    a, c, _, khi, _, m, valid, base, tl, th = parts
    exponent, guide = _quadratic_exponent(a, c, khi, m, valid, base, tl, th)
    if method in ('tangent', 'both'):
        exponent = _tangent_exponent(parts, exponent, guide, passes)
    if method in ('direct', 'both'):
        exponent = _direct_exponent(parts, exponent)
    return exp_upper_from_nonnegative_lower(np.maximum(0.0, -exponent))

@contextmanager
def using(method='tangent', passes=2):
    old = stable_bounds.log_range_cells

    def selected(*args, **kwargs):
        return weighted_range(*args, **kwargs, method=method, passes=passes)
    stable_bounds.log_range_cells = selected
    try:
        yield
    finally:
        stable_bounds.log_range_cells = old

def install(method='tangent', passes=2):
    global METHOD, PASSES
    assert method in ('tangent', 'direct', 'both', 'quadratic')
    assert isinstance(passes, int) and 1 <= passes <= 4
    METHOD, PASSES = (method, passes)
    tangent_table()
    stable_bounds.log_range_cells = weighted_range
    import stable_cover
    if getattr(stable_cover, '_weighted_log_installed', False):
        return
    original_factory = stable_cover.weight_factory

    def factory(hi, N, kind, target=None):
        if kind != 'weighted_log':
            return original_factory(hi, N, kind, target)
        from gaussian_batch import GaussianBatchWeights
        cache = {}

        def choose(tauhi):
            bucket = min(20, max(1, round(20 * tauhi / hi)))
            if bucket not in cache:
                cache[bucket] = GaussianBatchWeights(hi, N, hi * bucket / 20)
            return cache[bucket]
        return choose
    stable_cover.weight_factory = factory
    stable_cover._weighted_log_installed = True
