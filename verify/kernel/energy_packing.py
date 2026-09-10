from fractions import Fraction
from functools import lru_cache
import math
from flint import arb, ctx
PRICE = Fraction(9, 8)
PADDED_ARGUMENT = Fraction(10001, 10000)
TAU_DECREASING_THRESHOLD = Fraction(15, 23)

def exact(value):
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)

def upper(value):
    if not value.is_finite():
        raise ArithmeticError('A finite outward enclosure is required')
    return 0.0 if value == 0 else math.nextafter(float(value.upper()), math.inf)

def _box(values):
    left, right = map(Fraction, values)
    if left > right:
        raise ValueError('A box interval is reversed')
    return (left, right)

def q_ball(x):
    return arb(8) * x.sinc() ** 4 / (3 * x.cos() ** 2)

@lru_cache(maxsize=1)
def derivative_cap_certificate():
    z = PADDED_ARGUMENT
    sin_lower = z - z ** 3 / 6 + z ** 5 / 120 - z ** 7 / 5040
    sin_upper = sin_lower + z ** 9 / math.factorial(9)
    cos_lower = 1 - z * z / 2 + z ** 4 / 24 - z ** 6 / 720 + z ** 8 / 40320 - z ** 10 / math.factorial(10)
    cos_upper = cos_lower + z ** 12 / math.factorial(12)
    assert 0 < sin_lower < sin_upper and 0 < cos_lower < cos_upper
    bracket = 2 * z * (1 + cos_upper * cos_upper) - sin_lower * cos_lower
    cap = 8 * sin_upper ** 3 * bracket / (9 * z ** 4 * cos_lower ** 3)
    assert 0 < cap < Fraction(29, 4)
    return dict(argument=str(z), rational_upper=str(cap), strict_upper='29/4', passed=True)

def _packed_sum(frequency, weight_cap, total_weight):
    frequency, weight_cap, total_weight = map(Fraction, (frequency, weight_cap, total_weight))
    if frequency < 0 or weight_cap <= 0 or total_weight < 0:
        raise ValueError('Invalid packing data')
    count = total_weight // weight_cap
    residual = total_weight - count * weight_cap
    assert 0 <= residual < weight_cap
    cap_argument = exact(frequency) * exact(weight_cap).root(3)
    value = exact(count * weight_cap) * q_ball(cap_argument)
    if residual:
        value += exact(residual) * q_ball(exact(frequency) * exact(residual).root(3))
    return (value, int(count), residual)

def packing_moment_upper(frequency, *, L, d, tau, precision_bits=128):
    frequency = Fraction(frequency)
    Llo, Lhi = _box(L)
    dlo, dhi = _box(d)
    taulo, tauhi = _box(tau)
    taulo, tauhi = (max(Fraction(0), taulo), min(tauhi, Lhi))
    if not (frequency >= 0 and 0 < Llo <= Lhi and (0 <= dlo <= dhi <= 1) and (dhi > 0) and (taulo <= tauhi)):
        raise ValueError('Invalid frequency or moment box')
    derivative_cap_certificate()
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        actual_argument = exact(frequency) * exact(dhi).sqrt()
        if not actual_argument <= 1:
            return None
        cap = Fraction(upper(exact(dhi) * exact(dhi).sqrt()))
        padded_argument = exact(frequency) * exact(cap).root(3)
        if not padded_argument <= exact(PADDED_ARGUMENT):
            return None
        decreasing = taulo >= TAU_DECREASING_THRESHOLD * Lhi
        total = taulo if decreasing else tauhi
        packed, count, residual = _packed_sum(frequency, cap, total)
        excess = Lhi - taulo
        square = exact(excess) * (packed + exact(PRICE * excess))
        result = 0.0 if excess == 0 else upper(square.sqrt())
        return dict(upper=result, coefficient_sum_upper=upper(packed), branch='tau-decreasing' if decreasing else 'separate-factor-maxima', L_used=str(Lhi), tau_coefficient_used=str(total), tau_excess_used=str(taulo), variance_cap_used=str(dhi), weight_cap_rational=str(cap), full_weights=count, residual_weight=str(residual), actual_argument_upper=upper(actual_argument), padded_argument_upper=upper(padded_argument), price=str(PRICE), normalizer='none', source_multiplier='s^2 D(s)/2', precision_bits=ctx.prec)
    finally:
        ctx.prec = previous

def packing_source_upper(frequency, deletion_cap, *, L, d, tau, precision_bits=128):
    deletion_cap = Fraction(deletion_cap)
    if deletion_cap < 0:
        raise ValueError('The external common deletion cap must be nonnegative')
    result = packing_moment_upper(frequency, L=L, d=d, tau=tau, precision_bits=precision_bits)
    if result is None:
        return None
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        return upper(exact(frequency) ** 2 * exact(deletion_cap) * arb(result['upper']) / 2)
    finally:
        ctx.prec = previous
