from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
from flint import arb
_round27 = Path(__file__).resolve().parent.parent
for _directory in (_round27.parent / '26', _round27 / 'logarithm'):
    if str(_directory) not in sys.path:
        sys.path.insert(0, str(_directory))
import weighted_log_bounds as weighted
from interval_bounds import al, up, down
CONTINUOUS = weighted.weighted_charge_lower
MAX_SUPPORT = 1 << 20

def _extremal_coordinates(Qlo, Shi):
    Qlo, Shi = np.broadcast_arrays(np.asarray(Qlo, dtype=float), np.asarray(Shi, dtype=float))
    valid = (Qlo > 0) & (Shi > 0) & np.isfinite(Qlo) & np.isfinite(Shi)
    Q = np.where(valid, Qlo, 1.0)
    S = np.where(valid, Shi, 1.0)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        S2lo = np.maximum(0.0, down(S * S))
        S2hi = up(S * S)
        nlo = np.floor(down(S2lo / Q))
        nhi = np.floor(up(S2hi / Q))
    valid &= (nlo == nhi) & (nlo >= 1) & (nlo <= MAX_SUPPORT)
    m = np.where(valid, nlo, 1.0)
    mp1 = m + 1
    Q = np.where(valid, Qlo, 0.1)
    S = np.where(valid, Shi, 0.4)
    S2lo = np.maximum(0.0, down(S * S))
    S2hi = up(S * S)
    Dlo = np.maximum(0.0, down(down(down(mp1 * Q) - S2hi) / m))
    Dhi = np.maximum(0.0, up(up(up(mp1 * Q) - S2lo) / m))
    rootlo = np.maximum(0.0, down(np.sqrt(Dlo)))
    roothi = up(np.sqrt(Dhi))
    alo = np.maximum(0.0, down(down(S + rootlo) / mp1))
    ahi = up(up(S + roothi) / mp1)
    blo = np.maximum(0.0, down(S - up(m * ahi)))
    valid &= np.isfinite(ahi) & (ahi < 0.5)
    return (m, np.where(valid, alo, 0.0), np.where(valid, blo, 0.0), np.where(valid, Dlo, 0.0), valid)

@lru_cache(maxsize=16)
def _gap_coefficients(order):
    return tuple((al(arb(2) ** (p - 1) / p) for p in range(3, order + 1)))

def _gap_lower(Qlo, Shi, m, a, b, D, order):
    assert isinstance(order, int) and 3 <= order <= 16
    den = np.where(Shi > 0, Shi, 1.0)
    mean = np.maximum(0.0, down(np.maximum(0.0, Qlo) / den))
    gap3 = np.maximum(0.0, down(down(down(down(m * a) * b) * D) / den))
    degree = order - 3
    h = [np.ones_like(a)] + [np.zeros_like(a) for _ in range(degree)]
    for x in (a, b, mean):
        for j in range(1, degree + 1):
            h[j] = np.maximum(0.0, down(h[j] + down(x * h[j - 1])))
    coefficients = _gap_coefficients(order)
    factor = np.zeros_like(a)
    for coefficient, polynomial in zip(coefficients, h):
        factor = np.maximum(0.0, down(factor + down(coefficient * polynomial)))
    return np.maximum(0.0, down(gap3 * factor))

def discrete_charge_lower(Qlo, Shi, terms=12, mode='full'):
    assert mode in ('full', 'cubic', 'gap8')
    Qlo, Shi = np.broadcast_arrays(np.asarray(Qlo, dtype=float), np.asarray(Shi, dtype=float))
    continuous = CONTINUOUS(Qlo, Shi, terms)
    m, a, b, D, valid = _extremal_coordinates(Qlo, Shi)
    if not np.any(valid):
        return continuous
    if mode == 'full':
        pair = np.stack((a, b), axis=-1)
        square = np.maximum(0.0, down(pair * pair))
        argument = np.maximum(0.0, down(2 * pair))
        phi_extra = np.maximum(0.0, down(weighted._phi_lower(argument, terms) - 1))
        extra = np.maximum(0.0, down(square * phi_extra))
        extra = np.maximum(0.0, down(down(m * extra[..., 0]) + extra[..., 1]))
        candidate = np.maximum(0.0, down(Qlo + extra))
    else:
        order = 3 if mode == 'cubic' else 8
        extra = _gap_lower(Qlo, Shi, m, a, b, D, order)
        candidate = np.maximum(0.0, down(continuous + extra))
    return np.maximum(continuous, np.where(valid, candidate, continuous))

@contextmanager
def using(mode='full'):
    assert mode in ('full', 'cubic', 'gap8')
    previous = weighted.weighted_charge_lower

    def selected(Qlo, Shi, terms=12):
        return discrete_charge_lower(Qlo, Shi, terms=terms, mode=mode)
    weighted.weighted_charge_lower = selected
    try:
        yield
    finally:
        weighted.weighted_charge_lower = previous
