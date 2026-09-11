from fractions import Fraction
import math
from flint import arb, ctx

def _rational(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('All input endpoints and caps must be finite')
    return Fraction.from_float(value)

def _arb(value):
    value = _rational(value)
    return arb(value.numerator) / arb(value.denominator)

def _upper(value):
    value = arb(value)
    if not value.is_finite():
        raise ArithmeticError('Nonfinite Arb enclosure')
    return 0.0 if value == 0 else math.nextafter(float(value.upper()), math.inf)

def _lower(value):
    value = arb(value)
    if not value.is_finite():
        raise ArithmeticError('Nonfinite Arb enclosure')
    return 0.0 if value == 0 else math.nextafter(float(value.lower()), -math.inf)

def _interval(value):
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise ValueError('Require a pair of ordered moment endpoints')
    lo, hi = map(_rational, value)
    if lo > hi:
        raise ValueError('Reversed interval')
    return (lo, hi)

def _nodes(values):
    result = []
    for value in values:
        pair = _interval(value) if isinstance(value, (tuple, list)) else (_rational(value),) * 2
        result.append(pair)
    if len(result) < 2 or result[0] != (0, 0):
        raise ValueError('At least two nodes, starting at exact zero, are required')
    for (ll, lh), (rl, rh) in zip(result[:-1], result[1:]):
        if not (0 <= ll <= lh <= rl <= rh and ll < rh):
            raise ValueError('Node enclosures must be ordered and nonoverlapping')
    return tuple(result)

def _p(x):
    if x == 0:
        return arb(1) / 2
    return (x * x.sin() - 2 * (x / 2).sin() ** 2) / (x * x * x.cos())

def _correction(z):
    return arb(0) if z == 0 else -z.cos().log() - z * z / 2

def _moment_upper(Lhi, taulo, tauhi):
    if Lhi == taulo == tauhi:
        return 0.0
    L, lo, hi = map(_arb, (Lhi, taulo, tauhi))

    def square(t):
        return (L - t) * (L + arb(5) * t / 3)
    candidates = [_upper(square(lo)), _upper(square(hi))]
    vertex = L / 5
    if not (vertex < lo or vertex > hi):
        candidates.append(_upper(arb(16) * L * L / 15))
    return _upper(arb(max(0.0, *candidates)).sqrt())

def _ratio_upper(left_hi, right_lo, tau):
    l, r, u = map(_arb, (left_hi, right_lo, tau))
    exponent = -(r * r - l * l) / 2 if tau == 0 else ((r * u).cos().log() - (l * u).cos().log()) / (u * u)
    return min(1.0, _upper(exponent.exp()))

def _center_offsets(left_lo, right_hi, taulo, tauhi, dhi):
    l, r, a, b, d = map(_arb, (left_lo, right_hi, taulo, tauhi, dhi))
    abar = arb(0) if taulo == 0 else _correction(l * a) / (a * a)
    alower = b * _correction(r * d.sqrt()) / (d * d.sqrt())
    lo = (-l * l / 2).exp() * (-alower).expm1()
    hi = (-r * r / 2).exp() * (-abar).expm1()
    return (_lower(lo), min(0.0, _upper(hi)))

def reference_bounds(nodes, deletion_caps, *, L, d, tau, max_argument=Fraction(3, 2), precision_bits=128):
    nodes = _nodes(nodes)
    deletion_caps = tuple(deletion_caps)
    if len(deletion_caps) + 1 != len(nodes):
        raise ValueError('Exactly one deletion cap per cell is required')
    Llo, Lhi = _interval(L)
    dlo, dhi = _interval(d)
    taulo, tauhi = _interval(tau)
    taulo, tauhi = (max(Fraction(0), taulo), min(Lhi, tauhi))
    if not (Llo > 0 and dhi > 0 and (taulo <= tauhi)):
        raise ValueError('Invalid moment box or empty feasible tau projection')
    if taulo * taulo > dhi:
        raise ValueError('Empty box: every feasible array has tau²<=d')
    cap_argument = None if max_argument is None else _rational(max_argument)
    caps = tuple((None if cap is None else _rational(cap) for cap in deletion_caps))
    if any((cap is not None and cap < 0 for cap in caps)):
        raise ValueError('Deletion caps must be nonnegative')
    previous_precision = ctx.prec
    try:
        ctx.prec = max(previous_precision, int(precision_bits), 128)
        if cap_argument is not None and (not 0 < _arb(cap_argument) < arb.pi() / 2):
            raise ValueError('Require 0<max_argument<pi/2')
        M = _moment_upper(Lhi, taulo, tauhi)
        E_upper, lower_offset, upper_offset, eligible = ([], [], [], [])
        endpoints, records, active = ([0.0], [], True)
        for index, ((ll, lh), (rl, rh), deletion) in enumerate(zip(nodes[:-1], nodes[1:], caps)):
            argument = _arb(rh) * _arb(dhi).sqrt()
            active = bool(active and argument < arb.pi() / 2 and (cap_argument is None or argument <= _arb(cap_argument)))
            if not active:
                E_upper.append(math.inf)
                lower_offset.append(-math.inf)
                upper_offset.append(math.inf)
                eligible.append(False)
                endpoints.append(math.inf)
                records.append(dict(index=index, eligible=False))
                continue
            if deletion is None:
                raise ValueError('An eligible reference cell needs a certified deletion cap')
            l, r, width = (_arb(ll), _arb(rh), _arb(rh - ll))
            source = _upper(arb(M) * r * r * _p(argument) * _arb(deletion))
            decay = _ratio_upper(lh, rl, taulo)
            u = _arb(taulo)
            k = l if taulo == 0 else (l * u).tan() / u
            exp_decay = (-k * width).exp()
            psi = width if k == 0 else -(-k * width).expm1() / k
            previous = endpoints[-1]
            whole = max(previous, _upper(exp_decay * arb(previous) + arb(source) * psi))
            endpoint = min(whole, _upper(arb(decay) * arb(previous) + arb(source) * psi))
            dl, du = _center_offsets(ll, rh, taulo, tauhi, dhi)
            E_upper.append(whole)
            lower_offset.append(dl)
            upper_offset.append(du)
            eligible.append(True)
            endpoints.append(endpoint)
            records.append(dict(index=index, eligible=True, source_upper=source, propagation_ratio_upper=decay, E_upper=whole, endpoint_upper=endpoint, lower_offset=dl, upper_offset=du))
        return dict(status='outward reference components for the supplied deletion caps', E_upper=E_upper, lower_offset=lower_offset, upper_offset=upper_offset, eligible=eligible, endpoint_upper=endpoints, M_upper=M, taulo_used=str(taulo), tauhi_used=str(tauhi), exact_grid=all((lo == hi for lo, hi in nodes)), cells=records, precision_bits=ctx.prec, full_smoothing_certificate=False)
    finally:
        ctx.prec = previous_precision

def quick_checks():
    import cmath
    import time
    began = time.process_time()
    old_precision = ctx.prec
    exact = (Fraction(0), Fraction(1, 8), Fraction(1, 4))
    box = dict(L=(0.59, 0.6), d=(0.35, 0.354), tau=(0.51, 0.53))
    result = reference_bounds(exact, (1.0, 1.0), **box)
    assert all(result['eligible']) and result['exact_grid']
    assert all((a >= b >= 0 for a, b in zip(result['E_upper'], result['endpoint_upper'][1:])))
    assert all((l <= u <= 0 for l, u in zip(result['lower_offset'], result['upper_offset'])))
    wide = reference_bounds(exact, (1.0, 1.0), L=(0.1, 0.6), d=(0.01, 1.0), tau=(0.0, 0.5))
    expected = 4 * 0.6 / math.sqrt(15)
    assert expected <= wide['M_upper'] < expected + 1e-14
    assert wide['upper_offset'] == [0.0, 0.0]
    boundary = reference_bounds(exact, (1.0, 1.0), L=(Fraction(1, 2),) * 2, d=(Fraction(1, 4),) * 2, tau=(Fraction(1, 2),) * 2)
    assert boundary['M_upper'] == 0 and boundary['E_upper'] == [0.0, 0.0]
    nondyadic = reference_bounds(exact, (1.0, 1.0), L=(Fraction(1, 3),) * 2, d=(Fraction(1, 9),) * 2, tau=(Fraction(1, 3),) * 2)
    assert nondyadic['M_upper'] == 0 and nondyadic['E_upper'] == [0.0, 0.0]
    e = Fraction(1, 10 ** 12)
    paired = ((0, 0), (exact[1] - e, exact[1] + e), (exact[2] - e, exact[2] + e))
    broad = reference_bounds(paired, (1.0, 1.0), **box)
    assert not broad['exact_grid']
    assert all((b >= a for a, b in zip(result['E_upper'], broad['E_upper'])))
    assert all((b <= a for a, b in zip(result['lower_offset'], broad['lower_offset'])))
    assert all((b >= a for a, b in zip(result['upper_offset'], broad['upper_offset'])))
    stopped = reference_bounds((0, 1, 2, 3), (1, None, None), L=(0.5, 0.6), d=(1.0, 1.0), tau=(0.5, 0.5))
    assert stopped['eligible'] == [True, False, False]
    assert math.isinf(stopped['E_upper'][1])
    zero = reference_bounds(exact, (0, 0), **box)
    assert zero['E_upper'] == [0.0, 0.0]
    law_edges = tuple((Fraction(j, 4) for j in range(5)))
    law_bounds = reference_bounds(law_edges, (1,) * 4, L=(0.6, 0.6), d=(0.25, 0.25), tau=(0.5, 0.5))
    rho = 1.2
    hypotenuse = (rho + math.sqrt(rho * rho + 8)) / 2
    skew = math.sqrt(hypotenuse * hypotenuse - 4)
    atom = (skew + hypotenuse) / 2
    probability = 1 / (1 + atom * atom)
    law_points = 0
    for j in range(4):
        for k in range(9):
            t = (j + k / 8) / 4
            f = (probability * cmath.exp(0.5j * t * atom) + (1 - probability) * cmath.exp(-0.5j * t / atom)) ** 4
            center, gaussian = (math.cos(t / 2) ** 4, math.exp(-t * t / 2))
            assert abs(f - center) <= law_bounds['E_upper'][j] + 2e-14
            assert law_bounds['lower_offset'][j] - 2e-14 <= center - gaussian
            assert center - gaussian <= law_bounds['upper_offset'][j] + 2e-14
            law_points += 1
    assert ctx.prec == old_precision
    return dict(all_passed=True, scope='Fixed checks of the reference recurrence components', exact_fraction_grid=True, paired_outward_grid=True, interior_M_maximum=True, exact_Rademacher_zero=True, zero_source=True, pole_prefix_fallback=True, M_upper=result['M_upper'], direct_law_points=law_points, precision_restored=True, cpu_seconds=time.process_time() - began)
