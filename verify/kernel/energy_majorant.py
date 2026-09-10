from fractions import Fraction
import math
import time
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

def direct_gap(x, y, price=1):
    s, c = (x.sin(), x.cos())
    A, B, D = ((x + s * c) / 2, s * s, x * x + 2 * x * s * c)
    q = arb(2) * B * B / 3
    slope = arb(2) * B * D / 3 - arb(4) * B * B / 9
    h = (y - 1) ** 2 * (y + arb(1) / 2)
    V = acb(s, y * c) * acb(0, x * y).exp() - acb(0, y) + A * (y * y - 1)
    return h * (q + slope * (y * y - 1) / 2 + exact(price) * x ** 4 * c * c * h / 4) - abs(V) ** 2

def tail_coefficients():
    F = Fraction
    coefficients = [F(1, 400) - 1, -4, F(-3, 200) - F(11, 2), F(1, 100) - 3, F(9, 400) - F(9, 16), F(-3, 100), F(1, 100)]
    return tuple((sum((coefficients[k] * math.comb(k, j) * 16 ** (k - j) for k in range(j, 7))) for j in range(7)))

def certify(xlo=1, xhi=1, *, price=1, cpu_seconds=3.0, max_cells=20000, precision_bits=128, retain_boxes=True):
    xlo, xhi, price = map(rational, (xlo, xhi, price))
    if not (Fraction(1, 2) <= xlo <= xhi <= 1 and price >= 1):
        raise ValueError('The retained tail proof requires 1/2<=x<=1 and price>=1')
    previous = ctx.prec
    began = time.process_time()
    accepted, pending, evaluated = ([], [], 0)
    cuts = tuple(map(Fraction, (0, Fraction(3, 4), Fraction(5, 4), 2, 4, 8, 16)))
    queue = [(xlo, xhi, l, r, 0) for l, r in zip(cuts[:-1], cuts[1:])]
    try:
        ctx.prec = max(previous, precision_bits, 128)
        while queue and evaluated < max_cells and (time.process_time() - began < cpu_seconds):
            xl, xr, yl, yr, depth = queue.pop()
            value = gap_quotient(interval(xl, xr), interval(yl, yr), price)
            evaluated += 1
            if value >= 0:
                accepted.append((xl, xr, yl, yr, lower(value)))
                continue
            if depth >= 48:
                pending.append((xl, xr, yl, yr, depth))
                continue
            if xl < xr and 32 * (xr - xl) > yr - yl:
                middle = (xl + xr) / 2
                queue.extend([(xl, middle, yl, yr, depth + 1), (middle, xr, yl, yr, depth + 1)])
            else:
                middle = (yl + yr) / 2
                queue.extend([(xl, xr, yl, middle, depth + 1), (xl, xr, middle, yr, depth + 1)])
        pending.extend(queue)
        coefficients = tail_coefficients()
        assert all((c > 0 for c in coefficients))
        complete = not pending
        result = dict(status='AUTHOR interval energy-majorant certificate' if complete else 'AUTHOR incomplete interval proposal; no scalar majorant established', complete=complete, xlo=str(xlo), xhi=str(xhi), price=str(price), compact_y_range=[0, 16], tail_shift_coefficients=list(map(str, coefficients)), evaluated_cells=evaluated, accepted_cells=len(accepted), pending_cells=len(pending), precision_bits=ctx.prec, cpu_seconds=time.process_time() - began, minimum_accepted_lower=min((box[4] for box in accepted), default=None))
        if retain_boxes:
            result['accepted_boxes'] = [dict(xlo=str(xl), xhi=str(xr), ylo=str(yl), yhi=str(yr), quotient_lower=lo) for xl, xr, yl, yr, lo in accepted]
            result['pending_boxes'] = [list(map(str, box[:4])) for box in pending]
        return result
    finally:
        ctx.prec = previous

def quick_checks():
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128)
        checks = 0
        for x in (Fraction(1, 2), Fraction(3, 4), Fraction(1)):
            for y in (Fraction(0), Fraction(3, 4), Fraction(1), Fraction(5, 4), Fraction(4), Fraction(16)):
                xx, yy = (exact(x), exact(y))
                difference = gap_quotient(xx, yy) * (yy - 1) ** 4 - direct_gap(xx, yy)
                assert difference.contains(0), (x, y, difference)
                checks += 1
        assert all((c > 0 for c in tail_coefficients()))
        return dict(all_passed=True, quotient_identity_points=checks, rational_tail_coefficients_positive=True)
    finally:
        ctx.prec = previous
