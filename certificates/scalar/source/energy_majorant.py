from fractions import Fraction
import math
from flint import arb, acb, ctx

def rational(v):
    return v if isinstance(v, Fraction) else Fraction(v)

def exact(v):
    v = rational(v)
    return arb(v.numerator) / arb(v.denominator)

def interval(a, b):
    return exact(a).union(exact(b))

def upper(v):
    return math.nextafter(float(v.upper()), math.inf)

def lower(v):
    return math.nextafter(float(v.lower()), -math.inf)

def exponential_remainders(v):
    iv = acb(0, v)
    small = max(abs(lower(v)), abs(upper(v))) <= 2
    if small or v.contains(0):
        K = 20 if small else 64
        value = acb(1) / math.factorial(K + 4)
        for k in range(K - 1, -1, -1):
            value = value * iv + acb(1) / math.factorial(k + 4)
        error = upper(abs(v) ** (K + 1) / math.factorial(K + 5))
        result = {4: value + acb(arb(0, error), arb(0, error))}
        for n in (3, 2, 1):
            result[n] = acb(1) / math.factorial(n) + iv * result[n + 1]
        return result
    result = {1: acb(v.sin() / v, 2 * (v / 2).sin() ** 2 / v)}
    for n in (1, 2, 3):
        result[n + 1] = (result[n] - acb(1) / math.factorial(n)) / iv
    return result

def gap_quotient(x, y, price=1):
    s, c = (x.sin(), x.cos())
    A = (x + s * c) / 2
    B, D = (s * s, x * x + 2 * x * s * c)
    slope = arb(2) * B * D / 3 - arb(4) * B * B / 9
    scaled_price = exact(price) * x ** 4 * c * c / 4
    rem = exponential_remainders(x * (y - 1))
    coefficients = {4: acb(-2 * x ** 4), 3: acb(2 * x ** 4, 2 * x ** 3 * (1 + c * c)), 2: acb(x * x * (4 * A * s * c + 2 * c * c), -x * x * (2 * A * (1 + 2 * c * c) - 2 * s * c)), 1: acb(-2 * x * A * c * c, -2 * x * A * s * c)}
    oscillation = sum((coefficients[n] * rem[n] for n in (1, 2, 3, 4)), acb(0))
    return scaled_price * (y + arb(1) / 2) ** 2 + slope * (y / 2 + arb(5) / 4) - A * A - oscillation.real

def tail_coefficients():
    F = Fraction
    coefficients = [F(1, 400) - 1, -4, F(-3, 200) - F(11, 2), F(1, 100) - 3, F(9, 400) - F(9, 16), F(-3, 100), F(1, 100)]
    return tuple((sum((coefficients[k] * math.comb(k, j) * 16 ** (k - j) for k in range(j, 7))) for j in range(7)))
