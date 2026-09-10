from dataclasses import dataclass
from fractions import Fraction
import math
from flint import arb, ctx
Q = Fraction
_CENTERED_MAX_TERMS = 128

def exact(value):
    if isinstance(value, Q):
        return value
    if isinstance(value, int):
        return Q(value)
    if isinstance(value, str):
        return Q(value)
    if isinstance(value, float) and math.isfinite(value):
        return Q.from_float(value)
    raise TypeError('Use integers, finite floats, rational strings, or Fraction')

def _a(value):
    value = exact(value)
    return arb(value.numerator) / value.denominator

def _lo(value):
    return Q(str(value.lower().fmpq()))

def _hi(value):
    return Q(str(value.upper().fmpq()))

def _iv(lo, hi):
    return _a(lo).union(_a(hi))

@dataclass(frozen=True)
class DiskProblem:
    coefficients: tuple
    x_interval: tuple
    center_moments: bool = False

    def __post_init__(self):
        if type(self.center_moments) is not bool:
            raise ValueError('center_moments must be an actual boolean')

    @classmethod
    def make(cls, coefficients, x, *, center_moments=False):
        if type(center_moments) is not bool:
            raise ValueError('center_moments must be an actual boolean')
        coefficients = tuple(map(exact, coefficients))
        if len(coefficients) != 5:
            raise ValueError('Require (beta,gamma,c0,c2,Lambda)')
        if coefficients[4] <= 0:
            raise ValueError('Require Lambda > 0')
        if isinstance(x, (tuple, list)):
            if len(x) != 2:
                raise ValueError('An x interval needs two endpoints')
            endpoints = tuple(map(exact, x))
        else:
            endpoints = (exact(x), exact(x))
        if center_moments:
            if not 0 <= endpoints[0] <= endpoints[1]:
                raise ValueError('Centered mode requires a closed nonnegative x interval')
            if endpoints[0] == 0 and coefficients[4] <= 1:
                raise ValueError('A zero-touching centered cell requires Lambda > 1')
        elif not 0 < endpoints[0] <= endpoints[1]:
            raise ValueError('Require a closed positive x interval; x=0 needs centered mode')
        return cls(coefficients, endpoints, center_moments)

@dataclass(frozen=True)
class _CenteredBounds:
    real: tuple
    imaginary: tuple

def _centered_components(y, x0, x1):
    from entire_target import series_components, TargetBounds
    x = _iv(x0, x1)
    argument = max(abs(_lo(x)), abs(_hi(x))) * max(abs(_lo(y)), abs(_hi(y)))
    if x0 > 0 and argument > 4:
        inv = _iv(1 / x1, 1 / x0)
        inv2 = inv * inv
        si, co = ((x * y).sin(), (x * y).cos())
        u = 2 * inv * (co + y * y - 1) - 2 * y * inv2 * si
        v = 2 * inv * si + 2 * y * inv2 * (co - 1) - 2 * y
        uy = -2 * (1 + inv2) * si - 2 * y * inv * co + 4 * y * inv
        vy = 2 * (1 + inv2) * co - 2 * y * inv * si - 2 * inv2 - 2
        uyy = -2 * (x + 2 * inv) * co + 2 * y * si + 4 * inv
        vyy = -2 * (x + 2 * inv) * si - 2 * y * co
        ux = -2 * (1 + y * y) * inv2 * co + (4 * y * inv2 * inv - 2 * y * inv) * si - 2 * (y * y - 1) * inv2
        vx = -2 * (1 + y * y) * inv2 * si + 2 * y * inv * co - 4 * y * inv2 * inv * (co - 1)
        return _CenteredBounds(TargetBounds(u, uy, uyy, ux), TargetBounds(v, vy, vyy, vx))
    if argument >= _CENTERED_MAX_TERMS + 1:
        raise ValueError('Centered entire target exceeds the 128-term argument cap')
    rounded_up = -(-argument.numerator // argument.denominator)
    terms = min(_CENTERED_MAX_TERMS, max(12, rounded_up + 8))
    return series_components(y, x, terms=terms)

def _p(problem, correction, y):
    _, _, c0, c2, lam = problem.coefficients
    return c0 + correction + c2 * y * y + lam * y * y * y

def _radius_minimum(problem, correction, lo=Q(0), hi=None):
    _, _, _, c2, lam = problem.coefficients
    points = [lo]
    if hi is not None:
        points.append(hi)
    critical = -2 * c2 / (3 * lam)
    if critical >= lo and (hi is None or critical <= hi):
        points.append(critical)
    return min((_p(problem, correction, y) for y in points))

def _tail(problem, correction, radius):
    beta, gamma, c0, c2, lam = problem.coefficients
    x0 = problem.x_interval[0]
    if problem.center_moments:
        x1 = problem.x_interval[1]
        if lam > 1:
            cubic = lam - 1
            quadratic = c2 - abs(gamma)
            linear = -(x1 * x1 / 2 + abs(beta))
            constant = c0 + correction + abs(gamma)
        else:
            bound_gamma = max(abs(gamma - 2 / x0), abs(gamma - 2 / x1))
            cubic = lam
            quadratic = c2 - bound_gamma
            linear = -(4 / (x0 * x0) + abs(beta + 2))
            constant = c0 + correction + bound_gamma - 2 / x0
        r = radius
        return (constant + linear * r + quadratic * r * r + cubic * r * r * r, linear + 2 * quadratic * r + 3 * cubic * r * r, quadratic + 3 * cubic * r, cubic)
    quadratic = c2 - abs(gamma)
    linear = -(4 / (x0 * x0) + abs(beta))
    constant = c0 + correction + abs(gamma) - 2 / x0
    r = radius
    return (constant + linear * r + quadratic * r * r + lam * r * r * r, linear + 2 * quadratic * r + 3 * lam * r * r, quadratic + 3 * lam * r, lam)

class _Evaluator:

    def __init__(self, problem, correction):
        self.problem, self.correction = (problem, correction)
        self.beta, self.gamma, c0, self.c2, self.lam = map(_a, problem.coefficients)
        self.c0 = c0 + _a(correction)

    def quantities(self, y, x0, x1):
        if self.problem.center_moments:
            components = _centered_components(y, x0, x1)
            u, uy, uyy, ux = components.real
            v, vy, vyy, vx = components.imaginary
        else:
            x, inv = (_iv(x0, x1), _iv(1 / x1, 1 / x0))
            z = x * y
            si, co = (z.sin(), z.cos())
            inv2 = inv * inv
            u = 2 * inv * co - 2 * y * inv2 * si
            v = 2 * inv * si + 2 * y * inv2 * (co - 1)
            uy = -2 * (1 + inv2) * si - 2 * y * inv * co
            vy = 2 * (1 + inv2) * co - 2 * y * inv * si - 2 * inv2
            uyy = -2 * (x + 2 * inv) * co + 2 * y * si
            vyy = -2 * (x + 2 * inv) * si - 2 * y * co
            ux = -2 * (1 + y * y) * inv2 * co + (4 * y * inv2 * inv - 2 * y * inv) * si
            vx = -2 * (1 + y * y) * inv2 * si + 2 * y * inv * co - 4 * y * inv2 * inv * (co - 1)
        r, s = (u - self.gamma * (y * y - 1), v - self.beta * y)
        ry, sy = (uy - 2 * self.gamma * y, vy - self.beta)
        ryy, syy = (uyy - 2 * self.gamma, vyy)
        p = self.c0 + self.c2 * y * y + self.lam * y * y * y
        py, pyy = (2 * self.c2 * y + 3 * self.lam * y * y, 2 * self.c2 + 6 * self.lam * y)
        h = p * p - r * r - s * s
        hy = 2 * (p * py - r * ry - s * sy)
        hyy = 2 * (py * py + p * pyy - ry * ry - r * ryy - sy * sy - s * syy)
        hx = -2 * (r * ux + s * vx)
        return (h, hy, hyy, hx)

    def point_lower(self, y, x0, x1):
        values = self.quantities(_a(y), x0, x1)
        result = _lo(values[0])
        if x0 != x1:
            xm, hx = ((x0 + x1) / 2, (x1 - x0) / 2)
            middle = self.quantities(_a(y), xm, xm)[0]
            result = max(result, _lo(middle) - hx * _hi(abs(values[3])))
        return result

    def correction_witness(self, x, y):
        yy = _a(y)
        if self.problem.center_moments:
            components = _centered_components(yy, x, x)
            r = components.real.value - self.gamma * (yy * yy - 1)
            s = components.imaginary.value - self.beta * yy
        else:
            xx, inv = (_a(x), _a(1 / x))
            si, co = ((xx * yy).sin(), (xx * yy).cos())
            r = 2 * inv * co - 2 * yy * inv * inv * si - self.gamma * (yy * yy - 1)
            s = 2 * inv * si + 2 * yy * inv * inv * (co - 1) - self.beta * yy
        rlo, slo = (max(Q(0), _lo(abs(r))), max(Q(0), _lo(abs(s))))
        modulus_lower = _lo(_a(rlo * rlo + slo * slo).sqrt())
        lower = max(Q(0), modulus_lower - _p(self.problem, Q(0), y))
        return dict(x=x, y=y, lower=lower)

def _additional_correction(lower_h, lower_p):
    if lower_h >= 0 and lower_p >= 0:
        return Q(0)
    radicand = max(Q(0), lower_p * lower_p - lower_h)
    correction = max(Q(0), -lower_p + _hi(_a(radicand).sqrt()))
    if lower_p + correction < 0 or lower_h + 2 * correction * lower_p + correction ** 2 < 0:
        raise ArithmeticError('The outward correction did not prove the quadratic')
    return correction

def _cell(evaluator, x0, x1, y0, y1):
    y = _iv(y0, y1)
    ym, hy = ((y0 + y1) / 2, (y1 - y0) / 2)
    hx = (x1 - x0) / 2
    h, first, second, frequency_derivative = evaluator.quantities(y, x0, x1)
    midpoint_first = evaluator.quantities(_a(ym), x0, x1)[1]
    d, curvature = (_hi(abs(midpoint_first)), _lo(second))
    if curvature > 0 and d <= curvature * hy:
        quadratic = -d * d / (2 * curvature)
    else:
        quadratic = -d * hy + curvature * hy * hy / 2
    lower_h = max(_lo(h), evaluator.point_lower(ym, x0, x1) + quadratic)
    if first >= 0:
        lower_h = max(lower_h, evaluator.point_lower(y0, x0, x1))
    if first <= 0:
        lower_h = max(lower_h, evaluator.point_lower(y1, x0, x1))
    lower_p = _radius_minimum(evaluator.problem, evaluator.correction, y0, y1)
    return dict(x=(x0, x1), y=(y0, y1), lower_h=lower_h, lower_p=lower_p, additional=_additional_correction(lower_h, lower_p), x_score=hx * _hi(abs(frequency_derivative)), y_score=hy * d + hy * hy * _hi(abs(second)) / 2, children=None, split_axis=None, split_at=None)
