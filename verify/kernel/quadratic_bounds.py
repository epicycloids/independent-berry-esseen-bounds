from pathlib import Path
import sys
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
if (HERE.parent / 'concave_bounds.py').is_file():
    sys.path.insert(0, str(HERE.parent))
from concave_bounds import concave_cells as previous_cells, CosineWeights, sqrt_nonnegative
from interval_bounds import up, au, ia, sqpos, mul, add, sub, divpos
from paired_bounds import rational_poly, alpha_upper

def quadratic_cells(thi, Bhi, dhi, taulo, tauhi):
    previous = previous_cells(thi, Bhi, dhi, taulo, tauhi)
    tau = min(Bhi, max(0.0, taulo))
    if tau <= 0 or Bhi <= 0:
        return previous
    x = up(np.asarray(thi) * au(arb(dhi).sqrt()))
    safe_x = np.minimum(x, 2.0)
    z = sqpos(ia(safe_x))
    p = rational_poly(z, [arb(1), arb(1) / 6, -arb(1) / 40, arb(1) / 1008])[1]
    a = alpha_upper(safe_x)
    q = up(1 + up(up(safe_x * safe_x) / 3))
    b = up(up(8 * up(p * p)) / 3)
    excess = max(0.0, au(arb(Bhi) - arb(tau)))
    a2, b2, q2 = (sqpos(ia(a)), sqpos(ia(b)), sqpos(ia(q)))
    v = add(ia(1.0), q2)
    excess2 = up(excess * excess)
    base_linear = up(up(b * tau) * excess)
    constant = up(up(a * tau) ** 2)
    quadratic = up(v[1] * excess2)
    D = np.maximum(b, mul(mul(ia(2.0), ia(a)), sqrt_nonnegative(v))[1])
    polynomial_squared = up(up(constant + up(up(D * tau) * excess)) + quadratic)
    polynomial = np.minimum(Bhi, up(np.sqrt(np.maximum(0.0, polynomial_squared))))
    branch_gap = sub(b2, mul(ia(4.0), mul(a2, v)))
    eligible = (x <= 2) & (a > 1e-100) & (branch_gap[0] > 0)
    d2 = sub(b2, mul(ia(4.0), a2))
    d = sqrt_nonnegative((np.maximum(d2[0], 1e-200), np.maximum(d2[1], 1e-200)))
    aq = mul(ia(a), ia(q))
    rho = divpos(sub(d, mul(ia(2.0), aq)), ia(b))
    eligible &= rho[0] > 0
    safe_rho = (np.maximum(rho[0], 1e-100), np.maximum(rho[1], 1e-100))
    scale = mul(ia(tau), divpos(ia(a), ia(q)))
    contact_low, contact_high = (mul(scale, safe_rho), divpos(scale, safe_rho))
    Z = sub(mul(ia(q), d), mul(ia(a), sub(q2, ia(1.0))))
    eligible &= Z[0] > 0
    safe_Z = (np.maximum(Z[0], 1e-200), np.maximum(Z[1], 1e-200))
    denominator = mul(ia(2.0), sqrt_nonnegative(mul(ia(a), safe_Z)))
    safe_denominator = (np.maximum(denominator[0], 1e-200), np.maximum(denominator[1], 1e-200))
    intercept = divpos(mul(ia(a), ia(b)), safe_denominator)[1]
    slope = divpos(add(mul(ia(q), d), mul(ia(2.0), ia(a))), safe_denominator)[1]
    line = up(up(tau * intercept) + up(excess * slope))
    first = up(np.sqrt(np.maximum(0.0, up(up(constant + base_linear) + excess2))))
    second = up(np.sqrt(np.maximum(0.0, up(base_linear + quadratic))))
    candidate = np.where(excess <= contact_low[0], first, np.where(excess >= contact_high[1], second, line))
    candidate = np.minimum(polynomial, np.where(eligible, candidate, polynomial))
    return np.minimum(previous, np.where(x <= 2, candidate, Bhi))

class QuadraticWeights(CosineWeights):
    prefactor = staticmethod(quadratic_cells)

def exact_majorant(a, b, u, v, e, capped=True):
    a, b, u, v, e = map(arb, (a, b, u, v, e))
    if not (a >= 0 and b >= 0 and (u >= 0) and (v > u) and (e >= 0)):
        raise ValueError('Require a,b,e>=0 and 0<=u<v.')
    second = (b * e + v * e * e).sqrt()
    if a == 0:
        return min(1 + e, second) if capped else second
    if b <= 2 * a * v.sqrt():
        value = a + v.sqrt() * e
        return min(1 + e, value) if capped else value
    w = (v - u).sqrt()
    s, r = (u / (v - u), b / (a * w))
    delta = (r * r - 4 * s).sqrt()
    rho = (delta - 2) / r
    left, right = (a * rho / w, a / (w * rho))
    if e <= left:
        value = (a * a + b * e + u * e * e).sqrt()
    elif e >= right:
        value = second
    else:
        value = (a * r + w * (delta + 2 * s) * e) / (2 * (delta + s - 1).sqrt())
    return min(1 + e, value) if capped else value

def exact_prefactor(t, B, variance_cap, tau):
    t, B, variance_cap, tau = map(arb, (t, B, variance_cap, tau))
    x = t * variance_cap.sqrt()
    if x > 2 or tau == 0:
        return B
    a = 2 * x / 3 - x ** 3 / 15 + x ** 5 / 420
    p = 1 + x * x / 6 - x ** 4 / 40 + x ** 6 / 1008
    q = 1 + x * x / 3
    return min(B, tau * exact_majorant(a, 8 * p * p / 3, 1, 1 + q * q, (B - tau) / tau))
