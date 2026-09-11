from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib
import importlib.util
import math
from pathlib import Path
from flint import arb, ctx
try:
    from energy_envelope import compressed_box_upper, REVIEWED_SHA256
except ModuleNotFoundError as error:
    if error.name != 'energy_envelope':
        raise
    from envelope import compressed_box_upper, REVIEWED_SHA256
REFERENCE_SHA256 = '7c9c4d9d58383780969ace001ac62112578225e457ce9466fea1bd2ec5e7c4a6'

@lru_cache(maxsize=1)
def reviewed_reference():
    try:
        module = importlib.import_module('reference_bounds')
    except ModuleNotFoundError as error:
        if error.name != 'reference_bounds':
            raise
        path = Path(__file__).resolve().parent.parent / 'reference_bounds.py'
        spec = importlib.util.spec_from_file_location('author_frozen_reference_bounds', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() != REFERENCE_SHA256:
        raise AssertionError('The reference helper differs from its reviewed source')
    return module

def _prefix_for_ball(argument, table, base):
    if argument <= arb(1) / 2 or not argument <= 1:
        return None
    for row in table['rows']:
        if argument <= base._arb(Fraction(row['right'])):
            return (row['prefix_Q_upper'], row['prefix_lambda_upper'])
    raise ValueError('Compiled energy table is missing a required prefix')

def energy_reference_bounds(nodes, deletion_caps, *, energy_table, L, d, tau, max_argument=Fraction(3, 2), precision_bits=128):
    if energy_table.get('source_sha256') != REVIEWED_SHA256:
        raise ValueError('Require a compiled envelope for the reviewed transformed target')
    base = reviewed_reference()
    nodes, deletion_caps = (tuple(nodes), tuple(deletion_caps))
    old = base.reference_bounds(nodes, deletion_caps, L=L, d=d, tau=tau, max_argument=max_argument, precision_bits=precision_bits)
    parsed = base._nodes(nodes)
    _, dhi = base._interval(d)
    taulo = Fraction(old['taulo_used'])
    radii, endpoints, cells = ([], [0.0], [])
    source_changes = 0
    previous_precision = ctx.prec
    try:
        ctx.prec = max(previous_precision, 128, int(precision_bits))
        for j, (((ll, lh), (rl, rh)), deletion) in enumerate(zip(zip(parsed[:-1], parsed[1:]), deletion_caps)):
            original = old['cells'][j]
            if not original['eligible']:
                radii.append(math.inf)
                endpoints.append(math.inf)
                cells.append(dict(original))
                continue
            old_source = original['source_upper']
            source, energy_source = (old_source, None)
            argument = base._arb(rh) * base._arb(dhi).sqrt()
            pair = _prefix_for_ball(argument, energy_table, base)
            if pair is not None:
                moment = compressed_box_upper(*pair, L=L, tau=tau, precision_bits=ctx.prec)
                energy_source = base._upper(base._arb(rh) ** 2 * base._arb(deletion) * arb(moment) / 2)
                source = min(old_source, energy_source)
            source_changes += source < old_source
            width, left, u = (base._arb(rh - ll), base._arb(ll), base._arb(taulo))
            k = left if taulo == 0 else (left * u).tan() / u
            exp_decay = (-k * width).exp()
            psi = width if k == 0 else -(-k * width).expm1() / k
            previous = endpoints[-1]
            whole = max(previous, base._upper(exp_decay * arb(previous) + arb(source) * psi))
            whole = min(whole, old['E_upper'][j])
            endpoint = min(whole, old['endpoint_upper'][j + 1], base._upper(arb(original['propagation_ratio_upper']) * arb(previous) + arb(source) * psi))
            radii.append(whole)
            endpoints.append(endpoint)
            row = dict(original)
            row.update(source_upper=source, old_source_upper=old_source, energy_source_upper=energy_source, source_branch='energy' if source < old_source else 'old', E_upper=whole, endpoint_upper=endpoint, old_E_upper=old['E_upper'][j])
            cells.append(row)
    finally:
        ctx.prec = previous_precision
    result = dict(old)
    result.update(status='energy and reference component bounds', E_upper=radii, endpoint_upper=endpoints, cells=cells, old_E_upper=old['E_upper'], old_endpoint_upper=old['endpoint_upper'], energy_improved_source_cells=int(source_changes), energy_improved_radius_cells=sum((eligible and new < prior for eligible, new, prior in zip(old['eligible'], radii, old['E_upper']))), M_upper_role='Moment factor shared by the cellwise energy bounds', reference_source_sha256=REFERENCE_SHA256, energy_source_sha256=REVIEWED_SHA256, energy_table_replayed_per_box=False)
    return result
