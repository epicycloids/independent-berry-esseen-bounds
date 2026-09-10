from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial, isfinite
from typing import NamedTuple
from flint import arb
Q = Fraction
_DERIVATIVES = ((0, 0), (0, 1), (0, 2), (1, 0))
_NAMES = ('value', 'dw', 'dww', 'dx')

def exact(value):
    if isinstance(value, Q):
        return value
    if isinstance(value, int):
        return Q(value)
    if isinstance(value, str):
        return Q(value)
    if isinstance(value, float) and isfinite(value):
        return Q.from_float(value)
    raise TypeError('An exact input must be an integer, rational string, Fraction, or finite float.')

def _point(value):
    value = exact(value)
    return arb(value.numerator) / value.denominator

def _ball(value):
    if isinstance(value, arb):
        result = value
    elif isinstance(value, (tuple, list)):
        if len(value) != 2:
            raise ValueError('An interval requires two endpoints.')
        lo, hi = map(exact, value)
        if lo > hi:
            raise ValueError('Interval endpoints are reversed.')
        result = _point(lo).union(_point(hi))
    else:
        result = _point(value)
    if not result.is_finite():
        raise ValueError('Require finite input enclosures.')
    return result

def _lower(value):
    return Q(str(value.lower().fmpq()))

def _upper(value):
    return Q(str(value.upper().fmpq()))

def _absolute_upper(value):
    return max(abs(_lower(value)), abs(_upper(value)))

def _powers(value, degree):
    result = [arb(1)]
    for _ in range(degree):
        result.append(result[-1] * value)
    return result

def _falling(degree, order):
    if degree < order:
        return 0
    result = 1
    for j in range(order):
        result *= degree - j
    return result

@lru_cache(maxsize=None)
def _families(n):
    sign = -1 if n % 2 else 1
    return ((0, Q(2 * sign, factorial(2 * n)), 2 * n - 1, 2 * n), (0, Q(-2 * sign, factorial(2 * n + 1)), 2 * n - 1, 2 * n + 2), (1, Q(2 * sign, factorial(2 * n + 1)), 2 * n, 2 * n + 1), (1, Q(2 * sign, factorial(2 * n)), 2 * n - 2, 2 * n + 1))

class TargetBounds(NamedTuple):
    value: arb
    dw: arb
    dww: arb
    dx: arb

@dataclass(frozen=True)
class SeriesComponents:
    real: TargetBounds
    imaginary: TargetBounds
    tail_bounds: dict
    terms: int
    ratio_bound: Q
    x_abs_upper: Q
    w_abs_upper: Q

def series_components(w, x, *, terms=12):
    if not isinstance(terms, int) or terms < 1:
        raise ValueError('Use a positive integer number of series terms.')
    W, X = (_ball(w), _ball(x))
    Xabs, Wabs = (_absolute_upper(X), _absolute_upper(W))
    ratio = (Xabs * Wabs / Q(terms + 1)) ** 2
    if ratio >= 1:
        raise ValueError('The geometric ratio is not below one; increase terms.')
    xp, wp = (_powers(X, 2 * terms), _powers(W, 2 * terms + 2))
    sums = [[arb(0) for _ in _DERIVATIVES] for _ in range(2)]
    for n in range(1, terms + 1):
        for component, coefficient, a, b in _families(n):
            for j, (dx, dw) in enumerate(_DERIVATIVES):
                factor = _falling(a, dx) * _falling(b, dw)
                if factor:
                    sums[component][j] += _point(coefficient * factor) * xp[a - dx] * wp[b - dw]
    tails = [[Q(0) for _ in _DERIVATIVES] for _ in range(2)]
    for component, coefficient, a, b in _families(terms + 1):
        for j, (dx, dw) in enumerate(_DERIVATIVES):
            factor = _falling(a, dx) * _falling(b, dw)
            if factor:
                first = abs(coefficient) * factor * Xabs ** (a - dx) * Wabs ** (b - dw)
                tails[component][j] += first / (1 - ratio)
    for component in range(2):
        for j in range(len(_DERIVATIVES)):
            tail = tails[component][j]
            if tail:
                error = _point(-tail).union(_point(tail))
                sums[component][j] += error
    return SeriesComponents(TargetBounds(*sums[0]), TargetBounds(*sums[1]), {name: (tails[0][j], tails[1][j]) for j, name in enumerate(_NAMES)}, terms, ratio, Xabs, Wabs)

def centered_target(w, x, C, S, *, terms=12):
    components = series_components(w, x, terms=terms)
    C, S = (_ball(C), _ball(S))
    return TargetBounds(*(C * real + S * imaginary for real, imaginary in zip(components.real, components.imaginary)))

def centered_tail_coefficients(coefficients, direction, x_abs_hi, radius):
    if len(coefficients) != 4 or len(direction) != 2:
        raise ValueError('Require (c0,c1,c2,Lambda) and (C,S).')
    coefficients = tuple(map(exact, coefficients))
    direction = tuple(map(exact, direction))
    x_abs_hi, radius = (exact(x_abs_hi), exact(radius))
    if x_abs_hi < 0 or radius <= 0 or coefficients[3] <= 0:
        raise ValueError('Require x_abs_hi>=0, radius>0, and Lambda>0.')
    c0, c1, c2, lam = map(_point, coefficients)
    C, S = map(_point, direction)
    norm = (C * C + S * S).sqrt()
    a2, a3 = (c2 - norm * _point(x_abs_hi), lam - norm)
    R = _point(radius)
    rows = []
    for sign in (-1, 1):
        rows.append(tuple(map(_lower, (c0 + sign * c1 * R + a2 * R * R + a3 * R * R * R, sign * c1 + 2 * a2 * R + 3 * a3 * R * R, a2 + 3 * a3 * R, a3))))
    return {'certified': all((value >= 0 for row in rows for value in row)), 'rows': tuple(rows), 'cubic_margin_lower': _lower(a3), 'radius': radius, 'x_abs_hi': x_abs_hi}
