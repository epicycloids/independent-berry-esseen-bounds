from fractions import Fraction
from functools import lru_cache
import math
from flint import arb, ctx
PRICE = Fraction(9, 8)
TAIL_START = 48

def rational(v):
    return v if isinstance(v, Fraction) else Fraction(v)

def exact(v):
    v = rational(v)
    return arb(v.numerator) / arb(v.denominator)

def interval(lo, hi):
    return exact(lo).union(exact(hi))

def upper(v):
    if not v.is_finite():
        raise ArithmeticError('Nonfinite interval')
    return 0.0 if v == 0 else math.nextafter(float(v.upper()), math.inf)

def lower(v):
    if not v.is_finite():
        raise ArithmeticError('Nonfinite interval')
    return 0.0 if v == 0 else math.nextafter(float(v.lower()), -math.inf)

def normalized_gap(x, y, price=PRICE):
    import energy_majorant as module
    c, S, u = (x.cos(), x.sinc(), (2 * x).sinc())
    a = (1 + u) / 2
    R = arb(2) * S * S * (1 + 2 * u) / 3 - arb(4) * S ** 4 / 9
    small_E3 = module.exponential_remainders(2 * x)[3].real
    K = -u * S * S - 4 * x * x * small_E3 * small_E3
    z = y - 1
    E = module.exponential_remainders(x * z)
    return exact(price) * c * c * (y + arb(1) / 2) ** 2 / 4 + R * (y / 2 + arb(5) / 4) + K - 2 * a * u * (z + 2) * E[2].real - 2 * (a * c * c * z * z + (a * (1 + 2 * c * c) - u) * z + 1) * E[3].real + 2 * (c * c * z * z + (1 + c * c) * z + 1) * E[4].real

def tail_coefficients():
    F = Fraction
    polynomial = (F(4, 21), F(8, 7), F(20, 21), F(1))
    square = [sum((polynomial[i] * polynomial[j] for i in range(4) for j in range(4) if i + j == k), F(0)) for k in range(7)]
    h_squared = (F(1, 4), F(0), F(-3, 2), F(1), F(9, 4), F(-3), F(1))
    coefficients = [PRICE * h_squared[k] - square[k] for k in range(7)]
    return tuple((sum((coefficients[k] * math.comb(k, j) * TAIL_START ** (k - j) for k in range(j, 7)), F(0)) for j in range(7)))

def q_upper(x):
    x = exact(x)
    return upper(arb(8) * x.sinc() ** 4 / (3 * x.cos() ** 2))
