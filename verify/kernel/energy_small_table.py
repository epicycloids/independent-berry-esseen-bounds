from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib
import importlib.util
import math
from pathlib import Path
from flint import arb, ctx
WORKER_SHA256 = '43b265e75e63be5d154d4a41897ee1bde670cd9fa1568838e513f2cc084115f6'
REMAINDER_SHA256 = 'f1725128c15ba70b0e51029293e22b1a1e77c589fd907adea196a3fc41821c46'
FAMILY = 'small_frequency_normalized_gap_uniform_tail_v1'
COMPILED_FAMILY = 'replayed_small_frequency_prefix_v1'
PRICE = Fraction(9, 8)

@lru_cache(maxsize=1)
def worker():
    try:
        module = importlib.import_module('energy_small_worker')
    except ModuleNotFoundError as error:
        if error.name != 'energy_small_worker':
            raise
        path = Path(__file__).resolve().parent / 'worker.py'
        spec = importlib.util.spec_from_file_location('small_frequency_author_worker', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() != WORKER_SHA256:
        raise AssertionError('The small-frequency worker differs from the recorded source')
    if module.REMAINDER_SHA256 != REMAINDER_SHA256:
        raise AssertionError('Unexpected exponential-remainder source')
    module.remainders()
    return module

def _cover_y(intervals):
    end = Fraction(0)
    for left, right in sorted(intervals):
        if left > end:
            return False
        end = max(end, right)
    return end == 48

def verify_record(record, *, precision_bits=128):
    module = worker()
    if record.get('certificate_family') != FAMILY or record.get('source_sha256') != WORKER_SHA256 or record.get('remainder_source_sha256') != REMAINDER_SHA256 or (record.get('complete') is not True) or (record.get('pending_cells') != 0) or record.get('pending_boxes'):
        raise ValueError('Require a complete retained certificate from the frozen sources')
    xl, xr, price = map(Fraction, (record['xlo'], record['xhi'], record['price']))
    if not (0 <= xl <= xr <= Fraction(1, 2) and price == PRICE):
        raise ValueError('Unexpected frequency interval or energy price')
    if record.get('compact_y_range') != [0, 48]:
        raise ValueError('Unexpected compact y domain')
    coefficients = module.tail_coefficients()
    if not all((v > 0 for v in coefficients)) or record.get('tail_shift_coefficients') != list(map(str, coefficients)):
        raise ValueError('The exact uniform-tail certificate differs')
    raw = record.get('accepted_boxes')
    if not raw or len(raw) != record.get('accepted_cells'):
        raise ValueError('The complete accepted rectangle partition must be retained')
    boxes = []
    for row in raw:
        a, b, lo, hi = map(Fraction, (row['xlo'], row['xhi'], row['ylo'], row['yhi']))
        bound = Fraction(row['normalized_gap_lower'])
        if not (xl <= a <= b <= xr and 0 <= lo < hi <= 48 and (bound >= 0)):
            raise ValueError('A rectangle or its lower bound is outside the claimed domain')
        boxes.append((a, b, lo, hi, bound))
    xnodes = sorted({xl, xr, *(v for box in boxes for v in box[:2])})
    checks = xnodes + [(a + b) / 2 for a, b in zip(xnodes[:-1], xnodes[1:])]
    for x in checks:
        if not _cover_y([(lo, hi) for a, b, lo, hi, _ in boxes if a <= x <= b]):
            raise ValueError('The accepted rectangles leave an exact coverage gap')
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        for a, b, lo, hi, bound in boxes:
            value = module.normalized_gap(module.interval(a, b), module.interval(lo, hi), PRICE)
            if not value >= 0:
                raise ArithmeticError('A stored rectangle fails the energy inequality on replay')
            if not module.exact(bound) <= value:
                raise ArithmeticError('A stored outward lower bound is overstated')
        return dict(verified=True, rectangles=len(boxes), xlo=str(xl), xhi=str(xr), price=str(PRICE), precision_bits=ctx.prec)
    finally:
        ctx.prec = previous

def compile_small_cover(cover, *, endpoint=Fraction(1, 2), precision_bits=128):
    endpoint = Fraction(endpoint)
    if not 0 < endpoint <= Fraction(1, 2):
        raise ValueError('Require a positive prefix endpoint at most 1/2')
    records = cover['rows'] if isinstance(cover, dict) else list(cover)
    records = sorted(records, key=lambda row: (Fraction(row['xlo']), Fraction(row['xhi'])))
    end = Fraction(0)
    for row in records:
        left, right = map(Fraction, (row['xlo'], row['xhi']))
        if not (left == end and left < right <= endpoint):
            raise ValueError('Frequency records must form a contiguous nonoverlapping prefix')
        end = right
    if end != endpoint:
        raise ValueError('The requested prefix endpoint is not covered')
    proofs = [verify_record(row, precision_bits=precision_bits) for row in records]
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        q_right = worker().q_upper(endpoint)
    finally:
        ctx.prec = previous
    return dict(status='small-frequency bounds with all rectangles verified', compiled_family=COMPILED_FAMILY, verified=True, verification='all_rectangles_replayed', left='0', right=str(endpoint), price=str(PRICE), lambda_upper=float(PRICE), q_formula='(8/3)*sinc(x)^4/cos(x)^2', q_right_upper=q_right, verified_scalar_records=len(proofs), verified_rectangles=sum((p['rectangles'] for p in proofs)), precision_bits=max((p['precision_bits'] for p in proofs)), source_sha256=WORKER_SHA256, remainder_source_sha256=REMAINDER_SHA256, replay_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), scope='Scalar energy-sum bounds on the small-frequency interval')

def _compiled_endpoint(compiled):
    if compiled.get('compiled_family') != COMPILED_FAMILY or compiled.get('verified') is not True or compiled.get('verification') != 'all_rectangles_replayed' or (compiled.get('source_sha256') != WORKER_SHA256) or (compiled.get('remainder_source_sha256') != REMAINDER_SHA256) or (Fraction(compiled['left']) != 0) or (Fraction(compiled['price']) != PRICE):
        raise ValueError('Require metadata returned by complete small-prefix replay')
    endpoint = Fraction(compiled['right'])
    if not 0 < endpoint <= Fraction(1, 2):
        raise ValueError('Invalid compiled prefix endpoint')
    return endpoint

def small_prefix_pair(argument, compiled, *, precision_bits=128):
    endpoint = _compiled_endpoint(compiled)
    module = worker()
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        if isinstance(argument, arb):
            if not argument.is_finite():
                raise ValueError('The frequency enclosure must be finite')
            if not (argument >= 0 and argument <= module.exact(endpoint)):
                return None
            x = arb(argument.upper())
        else:
            X = Fraction(argument)
            if not 0 <= X <= endpoint:
                return None
            x = module.exact(X)
        q = arb(8) * x.sinc() ** 4 / (3 * x.cos() ** 2)
        return (module.upper(q), float(PRICE))
    finally:
        ctx.prec = previous

def replace_small_prefix(old_table, compiled, *, precision_bits=128):
    if _compiled_endpoint(compiled) != Fraction(1, 2):
        raise ValueError('The replacement requires the entire [0,1/2] prefix')
    if old_table.get('source_sha256') != REMAINDER_SHA256:
        raise ValueError('Unexpected upper-frequency table source')
    result = deepcopy(old_table)
    rows = result['rows']
    if not rows or Fraction(rows[0]['left']) != 0 or Fraction(rows[0]['right']) != Fraction(1, 2):
        raise ValueError('The supplied table must begin with its [0,1/2] prefix row')
    end = Fraction(1, 2)
    for row in rows[1:]:
        left, right = map(Fraction, (row['left'], row['right']))
        if not (left == end and left < right <= 1):
            raise ValueError('The upper-frequency table has a gap or overlap')
        for key in ('Q_upper', 'lambda_upper'):
            if not (math.isfinite(row[key]) and row[key] >= 0):
                raise ValueError('Nonfinite or negative upper-frequency coefficient')
        end = right
    if end != 1:
        raise ValueError('The supplied upper-frequency table must end at one')
    Q, price = small_prefix_pair(Fraction(1, 2), compiled, precision_bits=precision_bits)
    rows[0] = dict(left='0', right='1/2', branch='small-energy', Q_upper=Q, lambda_upper=price, prefix_Q_upper=Q, prefix_lambda_upper=price)
    for row in rows[1:]:
        Q, price = (max(Q, row['Q_upper']), max(price, row['lambda_upper']))
        row['prefix_Q_upper'], row['prefix_lambda_upper'] = (Q, price)
    result['small_frequency_prefix'] = deepcopy(compiled)
    result['status'] = 'energy table with verified small-frequency bounds'
    return result
