from contextlib import contextmanager
from pathlib import Path
import sys
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
DUAL_DIRECTORY = HERE.parent / 'dual_logarithm'
if str(DUAL_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(DUAL_DIRECTORY))
import dual_log_bounds as dual
from interval_bounds import K, KLO, KHI, MPEAK, XPEAK_LO, al, au, up, down, sqpos, cubepos, m_upper, exp_upper_from_nonnegative_lower
XPEAK_HI = au(1 / (6 * K))

def endpoint_polynomial_upper(A, U, dhi):
    d32hi = au(arb(dhi) * arb(dhi).sqrt())
    return up(up(A * d32hi) - down(U * dhi))

def cubic_max_upper(a, b, dhi):
    a, b = np.broadcast_arrays(np.asarray(a, dtype=float), np.asarray(b, dtype=float))
    assert np.isfinite(dhi) and dhi >= 0
    assert np.all(np.isfinite(a) & np.isfinite(b))
    if dhi == 0:
        return np.zeros_like(a)
    d = arb(dhi)
    root_hi = au(d.sqrt())
    d32lo, d32hi = (max(0.0, al(d * d.sqrt())), au(d * d.sqrt()))
    signed_power = np.where(b >= 0, d32hi, d32lo)
    endpoint = up(up(a * dhi) + up(b * signed_power))
    candidate = (a > 0) & (b < 0)
    safe_a = np.where(candidate, a, 0.0)
    safe_b = np.where(candidate, b, -1.0)
    denominator_lo = np.maximum(0.0, down(-3.0 * safe_b))
    denominator_hi = up(-3.0 * safe_b)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        stationary_lo = np.maximum(0.0, down(down(2.0 * safe_a) / denominator_hi))
        stationary_hi = np.divide(up(2.0 * safe_a), denominator_lo, out=np.full_like(a, np.inf), where=denominator_lo > 0)
        stationary_hi = np.minimum(root_hi, up(stationary_hi))
        stationary_value = up(up(safe_a * up(stationary_hi * stationary_hi)) / 3.0)
    include = candidate & (stationary_lo <= root_hi)
    return np.maximum(0.0, np.maximum(endpoint, np.where(include, stationary_value, 0.0)))

def boundary_loss_upper(t, A, U, dhi):
    t, A, U = np.broadcast_arrays(np.asarray(t, dtype=float), A, U)
    _, t2hi = sqpos((t, t))
    t3lo, _ = cubepos((t, t))
    t3lo = np.maximum(0.0, t3lo)
    a = up(up(t2hi / 2.0) - U)
    negative_cubic_lo = np.maximum(0.0, down(down(2.0 * KLO) * t3lo))
    b = up(A - negative_cubic_lo)
    return cubic_max_upper(a, b, dhi)

def joint_loss_upper(tlo, thi, A, U, dhi):
    tlo, thi, A, U = np.broadcast_arrays(np.asarray(tlo, dtype=float), np.asarray(thi, dtype=float), A, U)
    assert np.all((0 <= tlo) & (tlo <= thi) & (A >= 0) & (U >= 0))
    if dhi <= 0:
        return np.zeros_like(tlo)
    left = boundary_loss_upper(tlo, A, U, dhi)
    right = boundary_loss_upper(thi, A, U, dhi)
    sqrt_lo = max(0.0, al(arb(dhi).sqrt()))
    sqrt_hi = au(arb(dhi).sqrt())
    ridge_possible = (down(tlo * sqrt_lo) <= XPEAK_HI) & (up(thi * sqrt_hi) >= XPEAK_LO)
    ridge = up(MPEAK + endpoint_polynomial_upper(A, U, dhi))
    return np.maximum(left, np.maximum(right, np.where(ridge_possible, ridge, 0.0)))

def joint_exponent_lower(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False, prices=None, return_details=False):
    if not deleted:
        return dual.dual_exponent_lower(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, False, prices, return_details)
    tlo, thi = np.broadcast_arrays(np.asarray(tlo, dtype=float), np.asarray(thi, dtype=float))
    assert np.all((0 <= tlo) & (tlo <= thi))
    zero = np.zeros_like(tlo)
    if Vlo <= 0 or Vhi <= 0 or dhi <= 0:
        details = dict(price=zero, support=zero, enabled=np.zeros_like(tlo, bool), positive=zero, negative=zero, old_loss=zero, joint_loss=zero)
        return (zero, details) if return_details else zero
    tau_lo, tau_hi = (max(0.0, taulo), min(Bhi, max(0.0, tauhi)))
    assert tau_lo <= tau_hi
    t2lo, t2hi = sqpos((tlo, thi))
    t3lo, t3hi = cubepos((tlo, thi))
    t2lo, t3lo = (np.maximum(0.0, t2lo), np.maximum(0.0, t3lo))
    t2hi, t3hi = (np.maximum(0.0, t2hi), np.maximum(0.0, t3hi))
    khi = up(KHI * t3hi)
    sigma = np.maximum(0.0, down(down(Vlo / 2) * t2lo))
    m = m_upper(thi, dhi)
    d32hi = au(arb(dhi) * arb(dhi).sqrt())
    argument = up(up(au(4 * K) * thi) * au(arb(dhi).sqrt()))
    best = np.full_like(tlo, -np.inf)
    chosen = {key: zero.copy() for key in ('price', 'support', 'positive', 'negative', 'old_loss', 'joint_loss')}
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        guide_tau = max(0.0, tau_lo - dhi ** 1.5)
        guide_V = max(0.0, Vhi - dhi)
        mean = np.divide(4 * float(K) * tlo ** 3 * guide_tau, thi ** 2 * guide_V, out=zero.copy(), where=thi > 0)
        mean = np.where(np.isfinite(mean), mean, 0.0)
    for row in dual.support_table():
        price, u, aa, bb, domain, _, _ = row
        if prices is not None and price not in prices:
            continue
        wanted = np.searchsorted(u, mean, side='right') - 1
        permitted = np.searchsorted(-domain, -argument, side='right') - 1
        index = np.minimum(np.clip(wanted, 0, len(u) - 1), np.maximum(0, permitted))
        valid = argument <= domain[index]
        A = np.maximum(0.0, down(aa[index] * t3lo))
        U = up(bb[index] * t2hi)
        coefficient = down(A + down((price - 1.0) * khi))
        chosen_tau = np.where(coefficient >= 0, tau_lo, tau_hi)
        shared = down(coefficient * chosen_tau)
        constant = up(up(up((1 + price) * khi) * Bhi) + up(U * Vhi))
        delta = np.maximum(0.0, up(up(A * d32hi) - down(U * dhi)))
        old_loss = up(m + delta)
        loss = np.minimum(old_loss, joint_loss_upper(tlo, thi, A, U, dhi))
        value = down(down(down(sigma + shared) - constant) - loss)
        value = np.where(valid, value, -np.inf)
        if return_details:
            improved = value > best
            values = dict(price=price, support=u[index], positive=A, negative=U, old_loss=old_loss, joint_loss=loss)
            for key, entry in values.items():
                chosen[key] = np.where(improved, entry, chosen[key])
        best = np.maximum(best, value)
    best = np.maximum(0.0, best)
    if return_details:
        return (best, dict(**chosen, enabled=best > 0))
    return best

def joint_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted=False, baseline='dual15', prices=None):
    if baseline not in ('dual15', 'original', 'none'):
        raise ValueError('baseline must be dual15, original, or none')
    if not deleted and baseline == 'dual15':
        return dual.dual_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, False, baseline='higher4')
    exponent = joint_exponent_lower(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted, prices=prices)
    result = exp_upper_from_nonnegative_lower(exponent)
    if baseline == 'none':
        return result
    if baseline == 'dual15':
        old = dual.dual_range(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted, baseline='higher4')
    else:
        old = dual.ORIGINAL(tlo, thi, Vlo, Vhi, Bhi, dhi, taulo, tauhi, deleted)
    return np.minimum(old, result)

@contextmanager
def using(baseline='dual15', prices=None):
    previous = dual.stable_bounds.log_range_cells
    dual.stable_bounds.log_range_cells = lambda *args, **kwargs: joint_range(*args, **kwargs, baseline=baseline, prices=prices)
    try:
        yield
    finally:
        dual.stable_bounds.log_range_cells = previous

def warmup():
    return dual.warmup()
