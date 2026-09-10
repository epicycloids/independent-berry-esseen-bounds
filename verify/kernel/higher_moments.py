from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
ROUND27 = HERE.parent
for path in (ROUND27.parent / '26', ROUND27, ROUND27 / 'logarithm', ROUND27 / 'entropy_moments'):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
import stable_bounds
from interval_bounds import K, KLO, KHI, al, au, up, down, sqpos, cubepos, m_upper, exp_upper_from_nonnegative_lower
from moment_fast import square_sum_lower, tangent_charge
from weighted_log_bounds import weighted_charge_lower
from combined_moments import joint_charge
ORIGINAL = stable_bounds.log_range_cells
SUPPORT_GRID = 1024
LOG_GRID = 1024
LOG_MAX_INDEX = 504
COEFFICIENT = {p: al(arb(2) ** (p - 1) / p) for p in range(3, 9)}

def power_lower(x, p):
    result = np.asarray(x, dtype=float)
    for _ in range(p - 1):
        result = np.maximum(0.0, down(result * x))
    return result

def power_upper(x, p):
    result = np.asarray(x, dtype=float)
    for _ in range(p - 1):
        result = up(result * x)
    return result

@lru_cache(maxsize=2)
def variance_support_table(p=3):
    assert p in (3, 4)
    a, b = (2 * p - 2, p)
    u0 = (arb(a) - (arb(a * b) / (a + b - 1)).sqrt()) / (a + b)
    last = int(np.floor(al(u0) * SUPPORT_GRID))
    u = np.arange(last + 1, dtype=float) / SUPPORT_GRID
    alpha = u ** (a - 1) * (1 - u) ** (b - 1) * (a - (a + b) * u)
    beta = u ** a * (1 - u) ** (b - 1) * (a - 1 - (a + b - 1) * u)
    lo, hi = (u.copy(), np.ones_like(u))
    for _ in range(44):
        middle = (lo + hi) / 2
        residual = middle ** a * (1 - middle) ** b - alpha * middle + beta
        lo = np.where(residual >= 0, middle, lo)
        hi = np.where(residual >= 0, hi, middle)
    domain = np.maximum(u, lo * (1 - 2.0 ** (-28)))
    domain[0] = 1.0
    positive = np.empty_like(u)
    negative = np.empty_like(u)
    for i, value in enumerate(u):
        v = arb(float(value))
        assert v <= u0
        aa = v ** (a - 1) * (1 - v) ** (b - 1) * (a - (a + b) * v)
        bb = v ** a * (1 - v) ** (b - 1) * (a - 1 - (a + b - 1) * v)
        assert aa >= 0 and bb >= 0
        endpoint = arb(float(domain[i]))
        assert v <= endpoint <= 1
        assert endpoint ** a * (1 - endpoint) ** b - aa * endpoint + bb >= 0
        positive[i] = max(0.0, al(aa / (arb(2) ** p * (4 * K) ** (2 * p - 3))))
        negative[i] = max(0.0, au(bb / (arb(2) ** p * (4 * K) ** (2 * p - 2))))
    assert np.all(domain[1:] <= domain[:-1])
    for value in (u, positive, negative, domain):
        value.flags.writeable = False
    return (u, positive, negative, domain)

def root_lower(value, p):
    value = np.maximum(0.0, np.asarray(value, dtype=float))
    assert p in (3, 4)
    proposed = np.cbrt(value) if p == 3 else np.sqrt(np.sqrt(value))
    proposed = np.maximum(0.0, down(proposed * (1 - 2.0 ** (-48))))
    valid = np.isfinite(proposed) & (power_upper(proposed, p) <= value)
    return np.where(valid, proposed, 0.0)

def power_sum_lower(tlo, thi, Vhi, Bhi, dhi, taulo, p=3, return_details=False):
    assert p in (3, 4)
    tlo, thi = np.broadcast_arrays(np.asarray(tlo, dtype=float), np.asarray(thi, dtype=float))
    assert np.all((0 <= tlo) & (tlo <= thi))
    if Vhi <= 0 or dhi <= 0:
        zero = np.zeros_like(tlo)
        details = dict(u=zero, valid=np.zeros_like(tlo, dtype=bool), R=zero)
        return (zero, details) if return_details else zero
    tau = min(Bhi, max(0.0, taulo))
    u, aa, bb, domain = variance_support_table(p)
    argument = up(up(au(4 * K) * thi) * au(arb(dhi).sqrt()))
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        mean = np.divide(4 * float(K) * tlo ** 3 * tau, thi ** 2 * Vhi, out=np.zeros_like(tlo), where=thi > 0)
        selection_finite = np.isfinite(mean)
        mean = np.where(selection_finite, mean, 0.0)
        wanted = np.floor(np.clip(mean, 0.0, u[-1]) * SUPPORT_GRID).astype(np.intp)
    permitted = np.searchsorted(-domain, -argument, side='right') - 1
    index = np.minimum(wanted, np.maximum(0, permitted))
    valid = selection_finite & (argument <= domain[index])
    t2hi = np.maximum(0.0, sqpos((thi, thi))[1])
    t3lo, t3hi = cubepos((tlo, thi))
    t3lo = np.maximum(0.0, t3lo)
    positive = np.maximum(0.0, down(t3lo * down(aa[index] * tau)))
    negative = up(t2hi * up(bb[index] * Vhi))
    R = np.where(valid, np.maximum(0.0, down(positive - negative)), 0.0)
    excess = max(0.0, au(arb(Bhi) - arb(tau)))
    loss = up(up(KHI * t3hi) * excess)
    norm = np.maximum(0.0, down(root_lower(R, p) - loss))
    Qp = power_lower(norm, p)
    if return_details:
        return (Qp, dict(u=u[index], valid=valid, R=R))
    return Qp

@lru_cache(maxsize=2)
def logarithm_support_table(p=3):
    assert p in (3, 4)
    aa, bb = (np.empty(LOG_MAX_INDEX + 1), np.empty(LOG_MAX_INDEX + 1))
    aa[0], bb[0] = (COEFFICIENT[p], 0.0)
    for i in range(1, LOG_MAX_INDEX + 1):
        q = arb(i) / LOG_GRID
        g = -(1 - 2 * q).log() / 2 - q
        derivative = 2 * q / (1 - 2 * q)
        for k in range(2, p):
            coefficient = arb(2) ** (k - 1) / k
            g -= coefficient * q ** k
            derivative -= k * coefficient * q ** (k - 1)
        a = (q * derivative - g) / ((p - 1) * q ** p)
        b = q ** (p - 1) * a - g / q
        assert a >= 0 and b >= 0
        aa[i], bb[i] = (max(COEFFICIENT[p], al(a)), max(0.0, au(b)))
    aa.flags.writeable = bb.flags.writeable = False
    return (aa, bb)

def higher_tangent_charge(Q2, Qp, S, m, p=3, Q3=None):
    assert p in (3, 4)
    Q2, Qp, S, m = np.broadcast_arrays(*[np.asarray(x, dtype=float) for x in (Q2, Qp, S, m)])
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        guide = np.divide(np.maximum(0.0, Qp - m ** p), S - m, out=np.zeros_like(Qp), where=S > m)
        location = np.sqrt(guide) if p == 3 else np.cbrt(guide)
        index = np.rint(np.clip(location * LOG_GRID, 0.0, LOG_MAX_INDEX)).astype(np.intp)
    aa, bb = logarithm_support_table(p)
    a, b = (aa[index], bb[index])
    positive = np.maximum(0.0, down(Q2 + down(a * Qp)))
    removed = up(up(m * m) + up(a * power_upper(m, p)))
    if p == 4:
        assert Q3 is not None
        positive = np.maximum(0.0, down(positive + down(COEFFICIENT[3] * Q3)))
        removed = up(removed + up(COEFFICIENT[3] * power_upper(m, 3)))
    loss = np.maximum(0.0, up(removed - down(b * m)))
    return np.maximum(0.0, down(down(positive - up(b * S)) - loss))

def cubic_gap_charge(Q2, Q3, S, m=0.0, last=3):
    assert last in (3, 8)
    Q2, Q3, S, m = np.broadcast_arrays(*[np.asarray(x, dtype=float) for x in (Q2, Q3, S, m)])
    Q = np.maximum(0.0, down(Q2 - power_upper(m, 2)))
    T = np.maximum(0.0, down(Q3 - power_upper(m, 3)))
    charge = weighted_charge_lower(Q, S)
    valid = S > 0
    den = np.where(valid, S, 1.0)
    base = up(up(Q * Q) / den)
    alternative = np.where(valid, T, 0.0)
    if last > 3:
        x2hi = up(Q / den)
        x3lo = np.maximum(0.0, down(np.sqrt(np.maximum(0.0, down(T / den)))))
    for p in range(3, last + 1):
        gap = np.maximum(0.0, down(alternative - base))
        charge = np.maximum(0.0, down(charge + down(COEFFICIENT[p] * gap)))
        if p < last:
            alternative = np.maximum(0.0, down(alternative * x3lo))
            base = up(base * x2hi)
    return charge

def higher_moment_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False, p=3, charge='tangent'):
    assert p in (3, 4) and charge in ('tangent', 'gap3', 'tail8')
    assert p == 3 or charge == 'tangent'
    tlo, thi = np.broadcast_arrays(np.asarray(tlo, dtype=float), np.asarray(thi, dtype=float))
    old = ORIGINAL(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted)
    if Vlo <= 0:
        return old
    Q2 = square_sum_lower(tlo, thi, Vhi, Bhi, dhi, taulo)
    m = m_upper(thi, dhi) if deleted else np.zeros_like(thi)
    t2lo, t2hi = sqpos((tlo, thi))
    t3lo, t3hi = cubepos((tlo, thi))
    t2lo, t3lo = (np.maximum(0.0, t2lo), np.maximum(0.0, t3lo))
    khi = up(KHI * t3hi)
    tau = min(Bhi, max(0.0, taulo))
    S = np.maximum(0.0, up(up(up(Vhi / 2) * t2hi) - down(down(down(2 * KLO) * t3lo) * tau)))
    penalty = np.maximum(tangent_charge(Q2, S, m), joint_charge(Q2, S, m, mode='continuous', joint_delete=True))
    Q3 = power_sum_lower(tlo, thi, Vhi, Bhi, dhi, taulo, p=3)
    if charge == 'tangent':
        Qp = Q3 if p == 3 else power_sum_lower(tlo, thi, Vhi, Bhi, dhi, taulo, p=4)
        extra = higher_tangent_charge(Q2, Qp, S, m, p=p, Q3=Q3)
    else:
        extra = cubic_gap_charge(Q2, Q3, S, m, last=3 if charge == 'gap3' else 8)
    penalty = np.maximum(penalty, extra)
    positive = up(khi * au(arb(Bhi) + arb(tauhi)))
    negative = down(down(Vlo / 2) * t2lo)
    exponent = up(up(up(positive + m) - negative) - penalty)
    return np.minimum(old, exp_upper_from_nonnegative_lower(np.maximum(0.0, -exponent)))

@contextmanager
def using(p=3, charge='tangent'):
    previous = stable_bounds.log_range_cells
    stable_bounds.log_range_cells = lambda *args, **kwargs: higher_moment_range(*args, **kwargs, p=p, charge=charge)
    try:
        yield
    finally:
        stable_bounds.log_range_cells = previous

def warmup(include_quartic=True):
    for p in (3, 4) if include_quartic else (3,):
        variance_support_table(p)
        logarithm_support_table(p)
