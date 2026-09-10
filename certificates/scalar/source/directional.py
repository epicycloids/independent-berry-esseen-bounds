from dataclasses import dataclass
from fractions import Fraction
import math
from flint import arb, ctx
Q = Fraction

def exact(value):
    if isinstance(value, Q):
        return value
    if isinstance(value, int):
        return Q(value)
    if isinstance(value, str):
        return Q(value)
    if isinstance(value, float) and math.isfinite(value):
        return Q.from_float(value)
    raise TypeError('Use integers, finite binary64 floats, rational strings, or Fraction.')

def _arb(value):
    value = exact(value)
    return arb(value.numerator) / value.denominator

def _lower(value):
    return Q(str(value.lower().fmpq()))

def _upper(value):
    return Q(str(value.upper().fmpq()))

def _interval(lo, hi):
    return _arb(lo).union(_arb(hi))

@dataclass(frozen=True)
class _X:
    value: arb
    inverse: arb | None

def _xinterval(lo, hi):
    return _X(_interval(lo, hi), None if lo == 0 else _interval(1 / hi, 1 / lo))

@dataclass(frozen=True)
class Problem:
    coefficients: tuple
    x_interval: tuple
    direction: tuple
    center_moments: bool = False

    @classmethod
    def make(cls, coefficients, x, direction, *, center_moments=False):
        coefficients = tuple(map(exact, coefficients))
        direction = tuple(map(exact, direction))
        if len(coefficients) != 4 or len(direction) != 2:
            raise ValueError('Require (c0, c1, c2, Lambda) and (C, S).')
        if isinstance(x, (tuple, list)):
            if len(x) != 2:
                raise ValueError('An x interval requires two endpoints.')
            x_interval = tuple(map(exact, x))
        else:
            x_interval = (exact(x), exact(x))
        if type(center_moments) is not bool:
            raise TypeError('center_moments must be an explicit boolean.')
        if not (0 <= x_interval[0] <= x_interval[1] and (x_interval[0] > 0 or center_moments)):
            raise ValueError('Zero frequency requires the centered target; frequencies must be nonnegative.')
        if coefficients[3] <= 0:
            raise ValueError('Require Lambda > 0.')
        return cls(coefficients, x_interval, direction, center_moments)

class _Evaluator:

    def __init__(self, problem):
        self.problem = problem
        self.c0, self.c1, self.c2, self.lam = map(_arb, problem.coefficients)
        self.C, self.S = map(_arb, problem.direction)

    def terms(self, w, x):
        angle = x.value * w
        si, co = (angle.sin(), angle.cos())
        return (self.C * co + self.S * si, self.S * co - self.C * si)

    def entire(self, w, x):
        from entire_target import centered_target
        size = _upper(abs(w)) * _upper(abs(x.value))
        terms = max(12, math.ceil(size))
        if terms > 128:
            raise ValueError('Split the frequency interval: zero-frequency evaluation exceeds 128 series terms.')
        return centered_target(w, x.value, self.C, self.S, terms=terms)

    def use_entire(self, w, x):
        return x.inverse is None or (self.problem.x_interval[0] == 0 and _upper(abs(w)) * _upper(abs(x.value)) <= 8)

    def difference(self, w, x):
        if self.use_entire(w, x):
            absolute = abs(w)
            polynomial = self.c0 + self.c1 * w + self.c2 * w * w + self.lam * absolute * absolute * absolute
            return polynomial - self.entire(w, x).value
        t0, t1 = self.terms(w, x)
        inv = x.inverse
        if self.problem.center_moments:
            z = x.value * w
            si, co = (z.sin(), z.cos())
            real = 2 * (co - 1) * inv + 2 * w * (z - si) * inv * inv
            imag = 2 * (si - z) * inv + 2 * w * (co - 1) * inv * inv
            target = self.C * real + self.S * imag
        else:
            target = 2 * t0 * inv + 2 * w * (t1 - self.S) * inv * inv
        absolute = abs(w)
        polynomial = self.c0 + self.c1 * w + self.c2 * w * w + self.lam * absolute * absolute * absolute
        return polynomial - target

    def dw(self, w, x):
        if self.use_entire(w, x):
            return self.c1 + 2 * self.c2 * w + 3 * self.lam * w * abs(w) - self.entire(w, x).dw
        t0, t1 = self.terms(w, x)
        inv = x.inverse
        if self.problem.center_moments:
            z = x.value * w
            si, co = (z.sin(), z.cos())
            real = -2 * si + 2 * (2 * z - si - z * co) * inv * inv
            imag = 2 * (co - 1) + 2 * (co - 1 - z * si) * inv * inv
            target = self.C * real + self.S * imag
        else:
            target = 2 * (1 + inv * inv) * t1 - 2 * self.S * inv * inv - 2 * w * t0 * inv
        return self.c1 + 2 * self.c2 * w + 3 * self.lam * w * abs(w) - target

    def dww(self, w, x):
        if self.use_entire(w, x):
            return 2 * self.c2 + 6 * self.lam * abs(w) - self.entire(w, x).dww
        t0, t1 = self.terms(w, x)
        if self.problem.center_moments:
            z = x.value * w
            si, co = (z.sin(), z.cos())
            real = -2 * x.value * co + 4 * (1 - co) * x.inverse + 2 * w * si
            imag = -2 * (x.value + 2 * x.inverse) * si - 2 * w * co
            target = self.C * real + self.S * imag
        else:
            target = -2 * (x.value + 2 * x.inverse) * t0 - 2 * w * t1
        return 2 * self.c2 + 6 * self.lam * abs(w) - target

    def target_dx(self, w, x):
        if self.use_entire(w, x):
            return self.entire(w, x).dx
        t0, t1 = self.terms(w, x)
        inv = x.inverse
        if self.problem.center_moments:
            z = x.value * w
            si, co = (z.sin(), z.cos())
            real = -2 * w * si * inv - 2 * (co - 1) * inv * inv + 2 * w * w * (1 - co) * inv * inv - 4 * w * (z - si) * inv * inv * inv
            imag = 2 * w * (co - 1) * inv - 2 * (si - z) * inv * inv - 2 * w * w * si * inv * inv - 4 * w * (co - 1) * inv * inv * inv
            return self.C * real + self.S * imag
        return 2 * w * t1 * inv - 2 * (1 + w * w) * t0 * inv * inv - 4 * w * (t1 - self.S) * inv * inv * inv

    def point_w_lower(self, w, x0, x1):
        X = _xinterval(x0, x1)
        W = _arb(w)
        direct = _lower(self.difference(W, X))
        if x0 == x1:
            return direct
        xm, hx = ((x0 + x1) / 2, (x1 - x0) / 2)
        centered = _lower(self.difference(W, _xinterval(xm, xm))) - hx * _upper(abs(self.target_dx(W, X)))
        return max(direct, centered)

    def tail(self, radius):
        if self.problem.x_interval[0] == 0:
            from entire_target import centered_tail_coefficients
            if radius < 1:
                raise ValueError('The centered tail requires radius >= 1.')
            return centered_tail_coefficients(self.problem.coefficients, self.problem.direction, self.problem.x_interval[1], radius)['rows']
        x0 = _arb(self.problem.x_interval[0])
        R = _arb(radius)
        norm = (self.C * self.C + self.S * self.S).sqrt()
        constant = 2 * norm / x0
        quadratic = self.c2
        if self.problem.center_moments:
            if radius < 1:
                raise ValueError('The centered tail requires radius >= 1.')
            q = _arb(_upper(2 * self.C * _xinterval(*self.problem.x_interval).inverse))
            constant -= q
            quadratic -= q
        rows = []
        for sign in (-1, 1):
            linear = 2 * (norm - sign * self.S) / (x0 * x0)
            if self.problem.center_moments:
                linear -= 2 * sign * self.S
            b1 = sign * self.c1 - linear
            rows.append(tuple(map(_lower, (self.c0 - constant + b1 * R + quadratic * R * R + self.lam * R ** 3, b1 + 2 * quadratic * R + 3 * self.lam * R * R, quadratic + 3 * self.lam * R, self.lam))))
        return tuple(rows)

@dataclass
class _Node:
    x0: Q
    x1: Q
    w0: Q
    w1: Q
    lower: Q
    sample_upper: Q
    x_score: Q
    w_score: Q
    children: tuple | None = None
    split_axis: str | None = None
    split_at: Q | None = None

    def as_dict(self):
        return {'x': list(map(str, (self.x0, self.x1))), 'w': list(map(str, (self.w0, self.w1))), 'lower': str(self.lower), 'sample_upper': str(self.sample_upper), 'children': self.children, 'split_axis': self.split_axis, 'split_at': None if self.split_at is None else str(self.split_at)}

def _node(evaluator, x0, x1, w0, w1):
    X, W = (_xinterval(x0, x1), _interval(w0, w1))
    xm, wm = ((x0 + x1) / 2, (w0 + w1) / 2)
    hx, hw = ((x1 - x0) / 2, (w1 - w0) / 2)
    f0 = evaluator.point_w_lower(wm, x0, x1)
    D = _upper(abs(evaluator.dw(_arb(wm), X)))
    second = evaluator.dww(W, X)
    L = _lower(second)
    if L > 0 and D <= L * hw:
        quadratic_minimum = -D * D / (2 * L)
    else:
        quadratic_minimum = -D * hw + L * hw * hw / 2
    lower = max(_lower(evaluator.difference(W, X)), f0 + quadratic_minimum)
    first = evaluator.dw(W, X)
    if first >= 0:
        lower = max(lower, evaluator.point_w_lower(w0, x0, x1))
    if first <= 0:
        lower = max(lower, evaluator.point_w_lower(w1, x0, x1))
    sample_upper = min((_upper(evaluator.difference(_arb(wm), _xinterval(xx, xx))) for xx in sorted(set((x0, xm, x1)))))
    if lower > sample_upper:
        raise ArithmeticError('A cell bound contradicts its enclosed sample.')
    x_score = hx * _upper(abs(evaluator.target_dx(W, X))) if hx else Q(0)
    w_score = hw * D + hw * hw * _upper(abs(second)) / 2
    return _Node(x0, x1, w0, w1, lower, sample_upper, x_score, w_score)
