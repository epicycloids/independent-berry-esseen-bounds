from fractions import Fraction
import math
from flint import arb, ctx
try:
    from energy_cover_worker import reviewed, REVIEWED_SHA256
except ModuleNotFoundError as error:
    if error.name != 'energy_cover_worker':
        raise
    from cover_worker import reviewed, REVIEWED_SHA256

def rational(v):
    return v if isinstance(v, Fraction) else Fraction(v)

def exact(v):
    v = rational(v)
    return arb(v.numerator) / arb(v.denominator)

def upper(v):
    if not v.is_finite():
        raise ArithmeticError('Nonfinite enclosure')
    return 0.0 if v == 0 else math.nextafter(float(v.upper()), math.inf)

def _old_pair(X):
    x = exact(X)
    p = arb(1) / 2 if x == 0 else (x * x.sin() - 2 * (x / 2).sin() ** 2) / (x * x * x.cos())
    price = 4 * p * p
    return (upper(arb(8) * price / 3), upper(price))

def _q_upper(X):
    x = exact(X)
    return upper(arb(8) * x.sin() ** 4 / (3 * x ** 4 * x.cos() ** 2))

def _cover_y(intervals):
    end = Fraction(0)
    for lo, hi in sorted(intervals):
        if lo > end:
            return False
        end = max(end, hi)
    return end == 16

def verify_scalar_record(record, precision_bits=128):
    module = reviewed()
    if record.get('source_sha256') != REVIEWED_SHA256 or not record.get('complete') or record.get('pending_cells') != 0 or record.get('pending_boxes'):
        raise ValueError('A completed certificate from the reviewed source is required')
    xl, xr, price = map(rational, (record['xlo'], record['xhi'], record['price']))
    if not (Fraction(1, 2) <= xl <= xr <= 1 and price in (1, Fraction(9, 8))):
        raise ValueError('Unexpected frequency range or energy price')
    if list(map(str, module.tail_coefficients())) != record['tail_shift_coefficients']:
        raise ValueError('Tail certificate mismatch')
    raw = record.get('accepted_boxes')
    if not raw or len(raw) != record['accepted_cells']:
        raise ValueError('Retain the accepted partition for independent replay')
    boxes = []
    for row in raw:
        a, b, l, r = map(rational, (row['xlo'], row['xhi'], row['ylo'], row['yhi']))
        if not (xl <= a <= b <= xr and 0 <= l < r <= 16):
            raise ValueError('A stored rectangle lies outside the requested domain')
        boxes.append((a, b, l, r))
    xnodes = sorted({xl, xr, *(v for box in boxes for v in box[:2])})
    checks = xnodes + [(a + b) / 2 for a, b in zip(xnodes[:-1], xnodes[1:])]
    for x in checks:
        if not _cover_y([(l, r) for a, b, l, r in boxes if a <= x <= b]):
            raise ValueError('The accepted partition has a gap')
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        for (a, b, l, r), stored in zip(boxes, raw):
            value = module.gap_quotient(module.interval(a, b), module.interval(l, r), price)
            if not value >= 0:
                raise ArithmeticError('A stored rectangle does not prove the majorant')
            if not arb(stored['quotient_lower']) <= value:
                raise ArithmeticError('A reported lower bound is overstated')
    finally:
        ctx.prec = previous
    return dict(verified=True, rectangles=len(boxes), xlo=str(xl), xhi=str(xr), price=str(price))

def compile_envelope(records, *, denominator=128, precision_bits=128):
    denominator = int(denominator)
    if denominator < 2 or denominator % 2:
        raise ValueError('Require an even denominator')
    usable, ignored = ([], 0)
    for record in records:
        if not record.get('complete'):
            ignored += 1
            continue
        proof = verify_scalar_record(record, precision_bits)
        usable.append((rational(proof['xlo']), rational(proof['xhi']), rational(proof['price'])))
    nodes = sorted({Fraction(1, 2), Fraction(1), *(Fraction(j, denominator) for j in range(denominator // 2, denominator + 1)), *(v for row in usable for v in row[:2])})
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        Q, price = _old_pair(Fraction(1, 2))
        rows = [dict(left='0', right='1/2', branch='old', Q_upper=Q, lambda_upper=price, prefix_Q_upper=Q, prefix_lambda_upper=price)]
        for left, right in zip(nodes[:-1], nodes[1:]):
            choices = [p for l, r, p in usable if l <= left and right <= r]
            if choices:
                qj, pj, branch = (_q_upper(right), upper(exact(min(choices))), 'new')
            else:
                qj, pj = _old_pair(right)
                branch = 'old'
            Q, price = (max(Q, qj), max(price, pj))
            rows.append(dict(left=str(left), right=str(right), branch=branch, Q_upper=qj, lambda_upper=pj, prefix_Q_upper=Q, prefix_lambda_upper=price))
        return dict(status='whole-frequency energy bounds', rows=rows, new_upper_interval_complete=all((row['branch'] == 'new' for row in rows[1:])), verified_scalar_records=len(usable), ignored_incomplete_records=ignored, source_sha256=REVIEWED_SHA256, precision_bits=ctx.prec, scope='Scalar energy-sum bounds over the frequency range')
    finally:
        ctx.prec = previous

def lookup_prefix(table, X):
    X = rational(X)
    if X <= Fraction(1, 2) or X > 1:
        return None
    for row in table['rows']:
        if X <= rational(row['right']):
            return (row['prefix_Q_upper'], row['prefix_lambda_upper'])
    raise ValueError('The table does not cover X')

def compressed_box_upper(Q, price, *, L, tau, precision_bits=128):
    Q, price = (rational(Q), rational(price))
    Llo, Lhi = map(rational, L)
    taulo, tauhi = map(rational, tau)
    taulo, tauhi = (max(Fraction(0), taulo), min(Lhi, tauhi))
    if not (Q >= 0 and price >= 0 and (0 < Llo <= Lhi) and (taulo <= tauhi)):
        raise ValueError('Invalid coefficients or moment box')
    candidates = [taulo, tauhi]
    if Q > price:
        vertex = (Q - 2 * price) * Lhi / (2 * (Q - price))
        if taulo <= vertex <= tauhi:
            candidates.append(vertex)
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        squares = []
        for t in candidates:
            if t == Lhi:
                squares.append(0.0)
            else:
                e = exact(Lhi - t)
                squares.append(upper(e * (exact(Q * t) + exact(price) * e)))
        return upper(arb(max(squares)).sqrt())
    finally:
        ctx.prec = previous

def moment_forcing_upper(X, table, *, L, tau, precision_bits=128):
    X = rational(X)
    if X < 0:
        raise ValueError('Require a nonnegative standardized frequency')
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        if not exact(X) < arb.pi() / 2:
            raise ValueError('The reference source cannot cross its pole')
        old_Q, old_price = _old_pair(X)
        old = compressed_box_upper(old_Q, old_price, L=L, tau=tau, precision_bits=ctx.prec)
        pair = lookup_prefix(table, X)
        new = None if pair is None else compressed_box_upper(*pair, L=L, tau=tau, precision_bits=ctx.prec)
        return dict(upper=old if new is None else min(old, new), old_upper=old, energy_upper=new, normalizer='none', source_multiplier='s^2 D(s)/2')
    finally:
        ctx.prec = previous
