from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib
import importlib.util
import math
from pathlib import Path
import time
from flint import arb, ctx
REMAINDER_SHA256 = 'cb18a4904895d7bae7f0a5fed5a71895c179f6192197e8b69b04762a2e5d508f'
PRICE = Fraction(9, 8)
TAIL_START = 48

@lru_cache(maxsize=1)
def remainders():
    try:
        module = importlib.import_module('energy_majorant')
    except ModuleNotFoundError as error:
        if error.name != 'energy_majorant':
            raise
        path = Path(__file__).resolve().parent.parent / 'majorant.py'
        spec = importlib.util.spec_from_file_location('small_frequency_reviewed_remainders', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() != REMAINDER_SHA256:
        raise AssertionError('The exponential-remainder kernel is not the reviewed source')
    return module

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
    module = remainders()
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

def certify(xlo=0, xhi=Fraction(1, 128), *, cpu_seconds=15.0, max_cells=400000, precision_bits=128, retain_boxes=True):
    xlo, xhi = map(rational, (xlo, xhi))
    if not 0 <= xlo <= xhi <= Fraction(1, 2) or cpu_seconds <= 0:
        raise ValueError('Require 0<=xlo<=xhi<=1/2 and a positive CPU budget')
    remainders()
    previous = ctx.prec
    began = time.process_time()
    cuts = tuple(map(Fraction, (0, Fraction(3, 4), Fraction(5, 4), 2, 4, 8, 16, 32, 48)))
    queue = [(xlo, xhi, l, r, 0) for l, r in zip(cuts[:-1], cuts[1:])]
    accepted, pending, evaluated = ([], [], 0)
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        while queue and evaluated < int(max_cells) and (time.process_time() - began < cpu_seconds):
            xl, xr, yl, yr, depth = queue.pop()
            value = normalized_gap(interval(xl, xr), interval(yl, yr))
            evaluated += 1
            if value >= 0:
                accepted.append((xl, xr, yl, yr, lower(value)))
                continue
            if depth >= 48:
                pending.append((xl, xr, yl, yr, depth))
                continue
            x_score = (xr - xl) * (8 + max(abs(yl - 1), abs(yr - 1)))
            y_score = (yr - yl) * (1 + xr)
            if xl < xr and x_score > y_score:
                m = (xl + xr) / 2
                queue.extend(((xl, m, yl, yr, depth + 1), (m, xr, yl, yr, depth + 1)))
            else:
                m = (yl + yr) / 2
                queue.extend(((xl, xr, yl, m, depth + 1), (xl, xr, m, yr, depth + 1)))
        pending.extend(queue)
        coefficients = tail_coefficients()
        assert all((v > 0 for v in coefficients))
        result = dict(status='AUTHOR small-frequency energy certificate' if not pending else 'AUTHOR incomplete small-frequency proposal; requested cell not proved', certificate_family='small_frequency_normalized_gap_uniform_tail_v1', complete=not pending, xlo=str(xlo), xhi=str(xhi), price=str(PRICE), compact_y_range=[0, TAIL_START], tail_shift_coefficients=list(map(str, coefficients)), evaluated_cells=evaluated, accepted_cells=len(accepted), pending_cells=len(pending), minimum_accepted_lower=min((row[4] for row in accepted), default=None), q_right_upper=q_upper(xhi), precision_bits=ctx.prec, cpu_seconds=time.process_time() - began, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), remainder_source_sha256=REMAINDER_SHA256, dispatch_performed=False)
        if retain_boxes:
            result['accepted_boxes'] = [dict(xlo=str(a), xhi=str(b), ylo=str(l), yhi=str(r), normalized_gap_lower=v) for a, b, l, r, v in accepted]
            result['pending_boxes'] = [list(map(str, row[:4])) for row in pending]
        return result
    finally:
        ctx.prec = previous

def run_cells(indices, denominator=128, *, cpu_seconds=110.0, per_cell_seconds=15.0, max_cells_per_cell=400000, retain_boxes=True, precision_bits=128):
    denominator = int(denominator)
    indices = tuple((int(i) for i in indices))
    if denominator < 2 or denominator % 2 or any((not 0 <= i < denominator // 2 for i in indices)):
        raise ValueError('Require an even grid denominator and indices in [0,denominator/2)')
    began = time.process_time()
    rows = []
    for i in indices:
        remaining = cpu_seconds - (time.process_time() - began)
        if remaining <= 0:
            break
        result = certify(Fraction(i, denominator), Fraction(i + 1, denominator), cpu_seconds=min(per_cell_seconds, remaining), max_cells=max_cells_per_cell, retain_boxes=retain_boxes, precision_bits=precision_bits)
        result.update(index=i, denominator=denominator)
        rows.append(result)
    return dict(status='AUTHOR bounded small-frequency worker; no integrated bound', rows=rows, requested_indices=list(indices), attempted_cells=len(rows), completed_cells=sum((row['complete'] for row in rows)), unattempted_indices=list(indices[len(rows):]), cpu_seconds=time.process_time() - began, dispatch_performed=False)

def quick_checks():
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128)
        module = remainders()
        positive_checks = zero_checks = 0
        for x in (Fraction(1, 128), Fraction(1, 4), Fraction(1, 2)):
            for y in (0, Fraction(3, 4), 1, Fraction(5, 4), 4, 48):
                xx, yy = (exact(x), exact(y))
                difference = normalized_gap(xx, yy) * xx ** 4 - module.gap_quotient(xx, yy, PRICE)
                assert difference.contains(0), (x, y, difference)
                positive_checks += 1
        for y in (0, Fraction(3, 4), 1, Fraction(5, 4), 4, 48):
            yy = exact(y)
            endpoint = (4 * yy + 1) / 144 + (exact(PRICE) - 1) * (yy + arb(1) / 2) ** 2 / 4
            assert (normalized_gap(arb(0), yy) - endpoint).contains(0)
            zero_checks += 1
        expected = ('51794446417/288', '354649842/7', '28964396957/7056', '181934075/1176', '43072585/14112', '5161/168', '1/8')
        assert tuple(map(str, tail_coefficients())) == expected
        return dict(all_passed=True, positive_x_identity_points=positive_checks, nontrivial_zero_endpoint_points=zero_checks, exact_uniform_tail_checked=True)
    finally:
        ctx.prec = previous
