from pathlib import Path
import importlib.util
import sys
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
for path in (HERE.parent / 'odd_refinement', HERE.parent.parent / '26'):
    if path.exists():
        sys.path.insert(0, str(path))
from interval_bounds import up, down, au, ia, sqpos, mul, add, sub, divpos
from paired_bounds import rational_poly, alpha_upper, CosineWeights
ODD_SOURCE = HERE.parent / 'odd_refinement' / 'odd_bounds.py'
if not ODD_SOURCE.is_file():
    ODD_SOURCE = HERE / 'odd_prefactor.py'
if not ODD_SOURCE.is_file():
    raise ImportError('The odd-prefactor module is missing.')
_odd_spec = importlib.util.spec_from_file_location('_round27_concave_odd_input', ODD_SOURCE)
_odd_input = importlib.util.module_from_spec(_odd_spec)
_odd_spec.loader.exec_module(_odd_input)
if not hasattr(_odd_input, 'exact_prefactor'):
    raise ImportError('The supplied odd-prefactor module lacks exact_prefactor.')
odd_cells = _odd_input.odd_cells
old_exact_prefactor = _odd_input.exact_prefactor

def sqrt_nonnegative(interval):
    return (down(np.sqrt(np.maximum(0.0, interval[0]))), up(np.sqrt(np.maximum(0.0, interval[1]))))

def concave_cells(thi, Bhi, dhi, taulo, tauhi):
    baseline = odd_cells(thi, Bhi, dhi, taulo, tauhi)
    tau = min(Bhi, max(0.0, taulo))
    if tau <= 0 or Bhi <= 0:
        return baseline
    x = up(np.asarray(thi) * au(arb(dhi).sqrt()))
    safe_x = np.minimum(x, 2.0)
    z = sqpos(ia(safe_x))
    p = rational_poly(z, [arb(1), arb(1) / 6, -arb(1) / 40, arb(1) / 1008])[1]
    a = alpha_upper(safe_x)
    q = up(1 + up(up(safe_x * safe_x) / 3))
    b = up(up(8 * up(p * p)) / 3)
    excess = max(0.0, au(arb(Bhi) - arb(tau)))
    w = mul(ia(q), ia(a))
    w2, b2 = (sqpos(w), sqpos(ia(b)))
    curvature_gap = sub(b2, mul(ia(8.0), w2))
    eligible = (x <= 2) & (a > 1e-100) & (curvature_gap[0] > 0)
    d2 = sub(b2, mul(ia(4.0), w2))
    d = sqrt_nonnegative((np.maximum(d2[0], 1e-200), np.maximum(d2[1], 1e-200)))
    rho = divpos(sub(d, mul(ia(2.0), w)), ia(b))
    eligible &= rho[0] > 0
    safe_rho = (np.maximum(rho[0], 1e-100), np.maximum(rho[1], 1e-100))
    scale = mul(ia(tau), divpos(ia(a), ia(q)))
    contact_low = mul(scale, safe_rho)
    contact_high = divpos(scale, safe_rho)
    denominator = mul(ia(2.0), sqrt_nonnegative(mul(w, d)))
    safe_denominator = (np.maximum(denominator[0], 1e-200), np.maximum(denominator[1], 1e-200))
    intercept = divpos(mul(ia(a), ia(b)), safe_denominator)[1]
    slope = divpos(mul(ia(q), add(d, mul(ia(2.0), w))), safe_denominator)[1]
    line = up(up(tau * intercept) + up(excess * slope))
    base_linear = up(up(b * tau) * excess)
    first_squared = up(up(up(a * tau) ** 2) + base_linear)
    quadratic = up(up(q * q) * up(excess * excess))
    branch1 = up(np.sqrt(np.maximum(0.0, up(first_squared + quadratic))))
    branch2 = up(np.sqrt(np.maximum(0.0, up(base_linear + up(2 * quadratic)))))
    candidate = np.where(excess <= contact_low[0], branch1, np.where(excess >= contact_high[1], branch2, line))
    candidate = np.minimum(Bhi, candidate)
    return np.minimum(baseline, np.where(eligible, candidate, baseline))

class ConcaveWeights(CosineWeights):
    prefactor = staticmethod(concave_cells)

def exact_majorant(a, b, q, e):
    a, b, q, e = map(arb, (a, b, q, e))
    f2 = (b * e + 2 * q * q * e * e).sqrt()
    if a == 0:
        return min(1 + e, f2)
    if b <= arb(8).sqrt() * q * a:
        return min(1 + e, a + arb(2).sqrt() * q * e)
    r = b / (q * a)
    delta = (r * r - 4).sqrt()
    rho = (delta - 2) / r
    left, right = (a * rho / q, a / (q * rho))
    if e <= left:
        value = (a * a + b * e + q * q * e * e).sqrt()
    elif e >= right:
        value = f2
    else:
        value = (a * r + q * (delta + 2) * e) / (2 * delta.sqrt())
    return min(1 + e, value)

def exact_prefactor(t, B, d, tau):
    t, B, d, tau = map(arb, (t, B, d, tau))
    x = t * d.sqrt()
    if x > 2 or tau == 0:
        return B
    a = 2 * x / 3 - x ** 3 / 15 + x ** 5 / 420
    p = 1 + x * x / 6 - x ** 4 / 40 + x ** 6 / 1008
    q = 1 + x * x / 3
    return min(B, tau * exact_majorant(a, 8 * p * p / 3, q, (B - tau) / tau))
