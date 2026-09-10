from pathlib import Path
import sys
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
QUADRATIC = HERE.parent / 'concave_majorant' / 'quadratic_variant'
if QUADRATIC.is_dir():
    sys.path.insert(0, str(QUADRATIC))
from quadratic_bounds import quadratic_cells, QuadraticWeights
from concave_bounds import sqrt_nonnegative
from interval_bounds import al, au, up, down, ia, add, sub, mul, divpos, sqpos
from paired_bounds import rational_poly, alpha_upper

def packing_data(tau, d):
    if tau <= 0 or d <= 0:
        return None
    mass = arb(d) * arb(d).sqrt()
    quotient = arb(tau) / mass
    lo, hi = (al(quotient), au(quotient))
    if not (np.isfinite(lo) and np.isfinite(hi) and (0 <= lo <= hi < 2.0 ** 52)):
        return None
    nlo, nhi = (int(np.floor(lo)), int(np.floor(hi)))
    if nlo != nhi:
        return None
    remainder = arb(tau) - nlo * mass
    if not (remainder >= 0 and remainder <= mass):
        return None
    return (nlo, au(mass), max(0.0, au(remainder)), max(0.0, au(remainder ** (arb(1) / 3))))

def support_upper(a, b, q, price):
    a, b, q, price = np.broadcast_arrays(a, b, q, price)
    a2, b2, q2, l2 = map(lambda x: sqpos(ia(x)), (a, b, q, price))
    gap1 = sub(l2, ia(1.0))
    gap2 = sub(l2, add(ia(1.0), q2))
    discriminant = sub(b2, mul(ia(4.0), a2))
    valid = (a >= 0) & (b > 0) & (q >= 1) & (price > 1) & (gap2[0] >= 0) & (discriminant[0] > 0) & (up(up(2 * a) * price) <= b)
    safe_gap = (np.maximum(0.0, gap1[0]), np.maximum(0.0, gap1[1]))
    safe_discriminant = (np.maximum(0.0, discriminant[0]), np.maximum(0.0, discriminant[1]))
    radical = sqrt_nonnegative(mul(safe_gap, safe_discriminant))
    numerator = add(b2, mul(ia(4.0), mul(a2, safe_gap)))
    denominator = mul(ia(2.0), add(mul(ia(price), ia(b)), radical))
    safe_denominator = (np.maximum(1e-200, denominator[0]), np.maximum(1e-200, denominator[1]))
    first = divpos(numerator, safe_denominator)[1]
    second_root = sqrt_nonnegative((np.maximum(0.0, gap2[0]), np.maximum(0.0, gap2[1])))
    second_denominator = mul(ia(2.0), add(ia(price), second_root))
    second = divpos(ia(b), second_denominator)[1]
    valid &= (denominator[0] > 1e-200) & (second_denominator[0] > 0)
    return (np.maximum(first, second), valid)

def variance_cells(thi, Bhi, dhi, taulo, tauhi, return_details=False):
    thi = np.asarray(thi, dtype=float)
    assert np.all(np.isfinite(thi) & (thi >= 0))
    baseline = quadratic_cells(thi, Bhi, dhi, taulo, tauhi)
    tau = min(Bhi, max(0.0, taulo))
    data = packing_data(tau, dhi)
    if data is None:
        details = dict(enabled=np.zeros_like(thi, bool), reason='uncertified packing quotient')
        return (baseline, details) if return_details else baseline
    count, mass_hi, remainder_hi, remainder_root_hi = data
    x = up(thi * au(arb(dhi).sqrt()))
    safe_x = np.minimum(x, 1.5)
    z = sqpos(ia(safe_x))
    a = alpha_upper(safe_x)
    p = rational_poly(z, [arb(1), arb(1) / 6, -arb(1) / 40, arb(1) / 1008])[1]
    b = up(up(8 * up(p * p)) / 3.0)
    q = up(1.0 + up(up(safe_x * safe_x) / 3.0))
    excess = max(0.0, au(arb(Bhi) - arb(tau)))
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        mean = excess / tau
        guessed = (b + 2 * mean) / (2 * np.sqrt(a * a + b * mean + mean * mean))
        upper_price = b / (2 * a)
        lower_price = np.sqrt(1 + q * q)
        price = np.minimum(upper_price * (1 - 2.0 ** (-20)), np.maximum(lower_price * (1 + 2.0 ** (-20)), guessed))
    price = np.where(np.isfinite(price), price, 2.0)
    full_intercept, valid_full = support_upper(a, b, q, price)
    remainder_x = np.minimum(safe_x, up(thi * remainder_root_hi))
    remainder_a = alpha_upper(remainder_x)
    remainder_intercept, valid_remainder = support_upper(remainder_a, b, q, price)
    derivative_upper = up(full_intercept + up(up(2 * safe_x) / 9.0))
    eligible = (x <= 1.5) & (a > 1e-100) & valid_full & valid_remainder & (derivative_upper <= price)
    packed_full = up(up(float(count) * mass_hi) * full_intercept)
    packed_remainder = up(remainder_hi * remainder_intercept)
    candidate = up(up(packed_full + packed_remainder) + up(price * excess))
    result = np.minimum(baseline, np.where(eligible, candidate, Bhi))
    if return_details:
        return (result, dict(enabled=eligible, count=count, mass_upper=mass_hi, remainder_upper=remainder_hi, price=price, positive=a, linear=b, q=q, full_intercept=full_intercept, remainder_intercept=remainder_intercept, candidate=candidate, baseline=baseline))
    return result

class VarianceResolvedWeights(QuadraticWeights):
    prefactor = staticmethod(variance_cells)
