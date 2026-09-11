from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib
import importlib.util
import math
from pathlib import Path
from flint import arb, ctx
REFERENCE_SHA256 = '7c9c4d9d58383780969ace001ac62112578225e457ce9466fea1bd2ec5e7c4a6'

@lru_cache(maxsize=1)
def reviewed_reference():
    try:
        module = importlib.import_module('reference_bounds')
    except ModuleNotFoundError as error:
        if error.name != 'reference_bounds':
            raise
        path = Path(__file__).resolve().parent.parent / 'reference_bounds.py'
        spec = importlib.util.spec_from_file_location('attained_maximum_frozen_reference', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() != REFERENCE_SHA256:
        raise AssertionError('The reference parser is not the reviewed frozen source')
    return module

def exact(value):
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)

def upper(value):
    if not value.is_finite():
        raise ArithmeticError('A finite outward enclosure is required')
    return 0.0 if value == 0 else math.nextafter(float(value.upper()), math.inf)

def lower(value):
    if not value.is_finite():
        raise ArithmeticError('A finite outward enclosure is required')
    return 0.0 if value == 0 else math.nextafter(float(value.lower()), -math.inf)

def _box(value):
    lo, hi = map(Fraction, value)
    if lo > hi:
        raise ValueError('A box interval is reversed')
    return (lo, hi)

def _geometry(d, tau):
    dlo, dhi = _box(d)
    taulo, tauhi = _box(tau)
    taulo, tauhi = (max(Fraction(0), taulo), min(Fraction(1), tauhi))
    if not 0 <= dlo <= dhi <= 1:
        raise ValueError('Require 0<=d_lo<=d_hi<=1')
    if taulo > tauhi or taulo * taulo > dhi or tauhi * tauhi < dlo ** 3:
        return None
    D = max(dlo, taulo * taulo)
    if D > dhi or tauhi * tauhi < D ** 3:
        return None
    if D == 1:
        M, m, branch = (arb(1), arb(0), 'single-variance-one')
    elif D == taulo * taulo:
        M, m, branch = (exact(taulo), exact(taulo), 'ordinary-jensen')
    elif taulo * taulo <= D ** 3:
        M, m, branch = (exact(D).sqrt(), arb(0), 'zero-residual-moment')
    else:
        M = exact(D).sqrt()
        m = exact(taulo * taulo - D ** 3) / (exact(1 - D) * (exact(taulo) + exact(D) * M))
        branch = 'attained-maximum-plus-residual'
    return dict(D=D, V=1 - D, M=M, m=m, dhi=dhi, taulo=taulo, tauhi=tauhi, metadata=dict(D0=str(D), tau_lower_used=str(taulo), tau_upper_used=str(tauhi), branch=branch, M_upper=upper(M), m_upper=upper(m)))

def _j(z):
    return arb(0) if z == 0 else -z.cos().log() - z * z / 2

def _jdelta_lower(left, right, x):
    if left == right or x.contains(0):
        return 0.0
    value = (_j(exact(right) * x) - _j(exact(left) * x)) / (x * x)
    return max(0.0, lower(value))

def _extra_lower(left, right, geometry):
    a = _jdelta_lower(left, right, geometry['M'])
    b = _jdelta_lower(left, right, geometry['m']) if geometry['V'] else 0.0
    return max(0.0, lower(exact(geometry['D']) * arb(a) + exact(geometry['V']) * arb(b)))

def _drift_lower(t, geometry):
    if t == 0:
        return 0.0
    time = exact(t)

    def factor(x):
        z = time * x
        return z.sinc() / z.cos()
    value = time * (exact(geometry['D']) * factor(geometry['M']) + exact(geometry['V']) * factor(geometry['m']))
    return max(0.0, lower(time), lower(value))

def _eligible(t, geometry):
    return bool(exact(t) * exact(geometry['dhi']).sqrt() < arb.pi() / 2)

def point_components(t, *, d, tau, precision_bits=128):
    t = Fraction(t)
    if t < 0:
        raise ValueError('A nonnegative frequency is required')
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        geometry = _geometry(d, tau)
        if geometry is None or not _eligible(t, geometry):
            return None
        correction = _extra_lower(0, t, geometry)
        center = upper((-exact(t * t) / 2 - arb(correction)).exp())
        return dict(center_upper=min(1.0, center), drift_lower=_drift_lower(t, geometry), correction_lower=correction, geometry=geometry['metadata'], precision_bits=ctx.prec)
    finally:
        ctx.prec = previous

def reference_components(nodes, *, d, tau, precision_bits=128):
    nodes = reviewed_reference()._nodes(tuple(nodes))
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128, int(precision_bits))
        geometry = _geometry(d, tau)
        if geometry is None:
            return dict(usable=False, reason='empty necessary moment projection', cells=[])
        cells = []
        for index, ((ll, lh), (rl, rh)) in enumerate(zip(nodes[:-1], nodes[1:])):
            if not _eligible(rh, geometry):
                cells.append(dict(index=index, eligible=False))
                continue
            k = _drift_lower(ll, geometry)
            extra = _extra_lower(lh, rl, geometry)
            ratio = min(1.0, upper((-exact(rl * rl - lh * lh) / 2 - arb(extra)).exp()))
            correction = _extra_lower(0, ll, geometry)
            offset = min(0.0, upper((-exact(rh * rh) / 2).exp() * (-arb(correction)).expm1()))
            cells.append(dict(index=index, eligible=True, drift_lower=k, propagation_ratio_upper=ratio, upper_offset=offset, left_correction_lower=correction, ratio_extra_lower=extra))
        return dict(usable=True, geometry=geometry['metadata'], cells=cells, precision_bits=ctx.prec)
    finally:
        ctx.prec = previous

def refine_reference(inherited, nodes, *, d, tau, precision_bits=128):
    nodes = tuple(nodes)
    base = reviewed_reference()
    parsed = base._nodes(nodes)
    parts = reference_components(nodes, d=d, tau=tau, precision_bits=precision_bits)
    result = dict(inherited)
    if not parts['usable']:
        result['maximum_variance_refinement'] = parts
        return result
    count = len(parsed) - 1
    if not (len(inherited['cells']) == len(inherited['E_upper']) == len(inherited['upper_offset']) == count and len(inherited['endpoint_upper']) == count + 1 and (inherited['endpoint_upper'][0] == 0)):
        raise ValueError('The inherited reference result has a different grid shape')
    radii, endpoints, offsets, cells = ([], [0.0], [], [])
    previous_precision = ctx.prec
    try:
        ctx.prec = max(previous_precision, 128, int(precision_bits))
        for j, (((ll, lh), (rl, rh)), improvement) in enumerate(zip(zip(parsed[:-1], parsed[1:]), parts['cells'])):
            original = inherited['cells'][j]
            row = dict(original)
            if not original['eligible'] or not improvement['eligible']:
                radii.append(inherited['E_upper'][j])
                endpoints.append(inherited['endpoint_upper'][j + 1])
                offsets.append(inherited['upper_offset'][j])
                cells.append(row)
                continue
            source = Fraction(original['source_upper'])
            if source < 0:
                raise ValueError('A nonnegative inherited source is required')
            width, k = (exact(rh - ll), arb(improvement['drift_lower']))
            decay = (-k * width).exp()
            psi = width if k == 0 else -(-k * width).expm1() / k
            previous = endpoints[-1]
            whole = max(previous, upper(decay * arb(previous) + exact(source) * psi))
            whole = min(whole, inherited['E_upper'][j])
            ratio = min(original['propagation_ratio_upper'], improvement['propagation_ratio_upper'])
            endpoint = min(whole, inherited['endpoint_upper'][j + 1], upper(arb(ratio) * arb(previous) + exact(source) * psi))
            offset = min(inherited['upper_offset'][j], improvement['upper_offset'])
            radii.append(whole)
            endpoints.append(endpoint)
            offsets.append(offset)
            row.update(E_upper=whole, endpoint_upper=endpoint, upper_offset=offset, propagation_ratio_upper=ratio, maximum_variance_drift_lower=improvement['drift_lower'])
            cells.append(row)
    finally:
        ctx.prec = previous_precision
    result.update(status='attained-maximum reference bound', E_upper=radii, endpoint_upper=endpoints, upper_offset=offsets, cells=cells, maximum_variance_refinement=parts, maximum_variance_improved_radius_cells=sum((a < b for a, b in zip(radii, inherited['E_upper']))), maximum_variance_improved_center_cells=sum((a < b for a, b in zip(offsets, inherited['upper_offset']))), source_upper_values_preserved=True, center_lower_preserved=True, reference_source_sha256=REFERENCE_SHA256)
    return result
