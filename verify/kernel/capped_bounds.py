from pathlib import Path
import sys
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'concave_majorant' / 'quadratic_variant'
if (SOURCE / 'quadratic_bounds.py').is_file():
    sys.path.insert(0, str(SOURCE))
from quadratic_bounds import quadratic_cells as previous_cells, CosineWeights, exact_majorant as uncapped_reference, sqrt_nonnegative
from interval_bounds import up, au, ia, add, sub, mul, divpos, sqpos
from paired_bounds import rational_poly, alpha_upper

def _positive(interval):
    return (np.maximum(interval[0], 1e-200), np.maximum(interval[1], 1e-200))

def capped_cells(thi, Bhi, dhi, taulo, tauhi):
    previous = previous_cells(thi, Bhi, dhi, taulo, tauhi)
    tau = min(Bhi, max(0.0, taulo))
    if tau <= 0 or Bhi <= 0:
        return previous
    x = up(np.asarray(thi) * au(arb(dhi).sqrt()))
    sx = np.minimum(x, 2.0)
    z = sqpos(ia(sx))
    p = rational_poly(z, [arb(1), arb(1) / 6, -arb(1) / 40, arb(1) / 1008])[1]
    a = alpha_upper(sx)
    b = up(up(8 * up(p * p)) / 3)
    q = up(1 + up(up(sx * sx) / 3))
    aa, bb, qq = (ia(a), ia(b), ia(q))
    a2, b2, q2 = (sqpos(aa), sqpos(bb), sqpos(qq))
    v = add(ia(1.0), q2)
    valid = (x <= 2) & (a >= 0) & (a < 1)
    excess = max(0.0, au(arb(Bhi) - arb(tau)))
    linear = up(up(b * tau) * excess)
    constant = up(up(a * tau) ** 2)
    first = up(np.sqrt(np.maximum(0.0, up(up(constant + linear) + up(excess * excess)))))
    S = sub(mul(aa, sub(bb, ia(2.0))), mul(qq, sub(ia(1.0), a2)))
    result = np.where(valid & (S[0] >= 0), first, Bhi)
    bm2 = sub(bb, ia(2.0))
    radicand = add(mul(bm2, bm2), mul(ia(4.0), q2))
    kden = add(bm2, sqrt_nonnegative(radicand))
    k = divpos(ia(2.0), _positive(kden))
    after_crossing = valid & (S[1] < 0) & (kden[0] > 1e-100) & (k[0] > 1e-100)
    T = add(ia(1.0), k)
    J = sub(mul(mul(ia(2.0), aa), sub(ia(1.0), aa)), mul(sub(bb, mul(ia(2.0), aa)), k))
    chord_slope = add(ia(1.0), divpos(sub(ia(1.0), aa), _positive(k)))
    chord = up(up(tau * a) + up(excess * chord_slope[1]))
    chord_valid = after_crossing & (J[0] >= 0)
    result = np.minimum(result, np.where(chord_valid, chord, Bhi))
    delta = sub(b2, mul(ia(4.0), a2))
    gap = sub(mul(q2, sqpos(k)), a2)
    U = add(a2, add(mul(bb, k), sqpos(k)))
    slope_den = mul(ia(2.0), U)
    numerator = add(mul(T, add(bb, mul(ia(2.0), k))), sqrt_nonnegative(mul(delta, gap)))
    slope = divpos(numerator, _positive(slope_den))
    intercept = sub(T, mul(k, slope))
    corner_slope = divpos(add(bb, mul(ia(2.0), mul(v, k))), _positive(mul(ia(2.0), T)))
    contact_num = sub(bb, mul(ia(2.0), mul(intercept, slope)))
    contact_den = mul(ia(2.0), sub(sqpos(slope), ia(1.0)))
    contact = divpos(contact_num, _positive(contact_den))
    crossing = divpos(aa, qq)
    tangent_valid = after_crossing & (J[1] < 0) & (delta[0] > 0) & (gap[0] > 0) & (slope_den[0] > 1e-100) & (slope[0] > 1) & (intercept[0] >= 0) & (slope[1] <= corner_slope[0]) & (contact_num[0] > 0) & (contact_den[0] > 1e-100) & (contact[1] <= crossing[0])
    contact_scaled = mul(ia(tau), contact)
    tangent = up(up(tau * intercept[1]) + up(excess * slope[1]))
    candidate = np.where(excess <= contact_scaled[0], first, tangent)
    result = np.minimum(result, np.where(tangent_valid, candidate, Bhi))
    return np.minimum(previous, np.minimum(Bhi, result))

class CappedWeights(CosineWeights):
    prefactor = staticmethod(capped_cells)

def exact_majorant(a, b, q, e):
    a, b, q, e = map(arb, (a, b, q, e))
    if not (a >= 0 and a < 1 and (b >= 0) and (q > 0) and (e >= 0)):
        raise ValueError('Require 0<=a<1, b,e>=0, and q>0.')
    cap = 1 + e
    first = (a * a + b * e + e * e).sqrt()
    second = (b * e + (1 + q * q) * e * e).sqrt()
    if a == 0:
        return min(cap, second)
    S = a * (b - 2) - q * (1 - a * a)
    if S >= 0:
        return min(cap, first)
    if not S < 0:
        raise ValueError('Resolve the cap/crossing comparison at higher precision.')
    k = 2 / (b - 2 + ((b - 2) ** 2 + 4 * q * q).sqrt())
    J = 2 * a * (1 - a) - (b - 2 * a) * k
    if J >= 0:
        return min(cap, a + (1 + (1 - a) / k) * e)
    if not J < 0:
        raise ValueError('Resolve the origin-contact comparison at higher precision.')
    T = 1 + k
    U = a * a + b * k + k * k
    delta = b * b - 4 * a * a
    gap = q * q * k * k - a * a
    if not (delta > 0 and gap > 0):
        raise ValueError('Resolve the positive tangent gaps at higher precision.')
    M = (T * (b + 2 * k) + (delta * gap).sqrt()) / (2 * U)
    C = T - k * M
    n = (b + 2 * (1 + q * q) * k) / (2 * T)
    if M <= n:
        left = (b - 2 * C * M) / (2 * (M * M - 1))
        return min(cap, first if e <= left else C + M * e)
    if not M > n:
        raise ValueError('Resolve the corner/tangent comparison at higher precision.')
    return uncapped_reference(a, b, 1, 1 + q * q, e, capped=True)

def exact_prefactor(t, B, variance_cap, tau):
    t, B, variance_cap, tau = map(arb, (t, B, variance_cap, tau))
    x = t * variance_cap.sqrt()
    if x > 2 or tau == 0:
        return B
    a = 2 * x / 3 - x ** 3 / 15 + x ** 5 / 420
    p = 1 + x * x / 6 - x ** 4 / 40 + x ** 6 / 1008
    q = 1 + x * x / 3
    return min(B, tau * exact_majorant(a, 8 * p * p / 3, q, (B - tau) / tau))
