from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
PACKED_SHA = '5e2de769e7b13f28370974ad0ff77cf2672e322fa80dcd8386549a4441c1526e'
REFINE_SHA = '0fbd1217d20b2447f50039ccd7af8ba554a981ed27277ac77431bd0abe67054f'
MODES = ('baseline', 'identity', 'contracted')
REFINED_KEYS = ('E_upper', 'endpoint_upper', 'lower_offset', 'upper_offset', 'eligible', 'maximum_variance_refinement', 'precision_bits', 'source_upper_values_preserved', 'center_lower_preserved')

def require(value, message):
    if not value:
        raise ValueError(message)

def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def encode(value):
    if isinstance(value, float) and (not math.isfinite(value)):
        require(not math.isnan(value), 'NaN is not reference evidence')
        return {'nonfinite': 'positive_infinity' if value > 0 else 'negative_infinity'}
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    if hasattr(value, 'tolist'):
        return encode(value.tolist())
    require(value is None or isinstance(value, (str, int, float, bool)), 'Unexpected reference evidence type')
    return value

def certificate(Llo, Lhi, box):
    dl, dh, bl, bh, a, z = map(Q, box)
    lo, M = (Q(Llo), Q(Lhi))
    require(0 < lo <= M and 0 < dl <= dh <= 1 and (0 <= bl <= bh) and (0 <= a <= z), 'Invalid original moment rectangle')
    scale = 2 ** 64
    kroot = math.isqrt(dh.numerator * scale * scale // dh.denominator)
    require(Q(kroot * kroot, scale * scale) <= dh < Q((kroot + 1) ** 2, scale * scale), 'Invalid square-root enclosure')
    W0 = dh * Q(kroot, scale)
    Wupper = dh * Q(kroot + 1, scale)
    balance = bh * a / M
    c = (balance + W0) / 2
    E = M - a
    rest = a - Wupper - Q(15, 23) * (M - bl)
    gamma = Q(1, 4)
    require(0 < balance < c < W0 and c * c > dl ** 3 and (Wupper <= bl) and (rest > 0), 'The complete actual b rectangle does not support this analytic cut')
    count = a // W0
    residual = a - count * c
    require(count == a // c == 2 and 0 <= residual < c and (gamma ** 3 >= residual ** 2) and (dl > gamma), 'The complete weight interval is not in the fixed packing branch')
    C = Q(40 * count, 27) * (dl - gamma)
    A = Q(29, 4) * a + Q(9, 8) * E
    K = C * (W0 - c) / A
    delta = M * c - bh * a
    D = E * (a + Q(27, 64) * E)
    require(0 < delta <= D and K > 0, 'No packed-source defect')
    H = (delta / D) ** 2
    return dict(kind='packed_weight_defect.v1', L=[str(lo), str(M)], actual_box=list(map(str, map(Q, box))), weight_endpoint_lower=str(W0), weight_endpoint_upper=str(Wupper), balance_weight=str(balance), analytic_weight_cut=str(c), endpoint_transport_margin=str(rest), packing_count=count, residual_upper=str(residual), residual_two_thirds_upper=str(gamma), derivative_gap_coefficient=str(C), weight_gap=str(W0 - c), denominator_upper=str(A), frequency_squared_loss_coefficient=str(K), positive_defect=str(delta), Cauchy_denominator=str(D), high_weight_loss=str(H), packed_source_only=True, final_radius_scaling=False, reference_center_unchanged=True)

def source_factor(proof, left, mode):
    require(mode in MODES and Q(left) >= 0, 'Invalid source mode or left frequency')
    if mode != 'contracted':
        return Q(1)
    loss = min(Q(proof['frequency_squared_loss_coefficient']) * Q(left) ** 2, Q(proof['high_weight_loss']))
    factor = 1 - loss / 2
    require(0 < factor <= 1 and factor * factor + loss >= 1, 'Invalid rational packed-source factor')
    return factor

def refined_fields(result):
    return encode({key: result[key] for key in REFINED_KEYS if key in result})

def check_receipt(receipt, Llo, Lhi, box, reference_box, mode):
    proof = certificate(Llo, Lhi, box)
    require(receipt['certificate'] == proof and receipt['mode'] == mode and (mode in MODES), 'Changed packed-source proof or operation')
    require(receipt['actual_box'] == box and receipt['reference_box'] == reference_box and (box[:2] == reference_box[:2]) and (box[4:] == reference_box[4:]) and (reference_box[2] <= box[2] <= box[3] <= reference_box[3]), 'Changed original scoring/reference domain')
    require(receipt['packed_function_sha256'] == PACKED_SHA and receipt['refine_function_sha256'] == REFINE_SHA and receipt['all_three_caches_cleared_before_and_after'] and receipt['original_bindings_restored'], 'Changed recurrence source or cache isolation')
    columns = receipt['columns']
    require(len(columns) == 3, 'Missing source column')
    counts = dict(eligible_cells=0, packed_available_cells=0, original_packed_minimizers=0, strictly_reduced_packed_candidates=0, strictly_reduced_selected_sources=0)
    for column in columns:
        pre = column['packed_reference']
        post = column['refined_reference']
        nodes = column['node_enclosures']
        original = column['original_packed_candidates']
        factors = column['source_factors']
        require(len(nodes) == 513 and len(original) == len(factors) == len(column['deletion_caps']) == len(pre['cells']) == 512, 'Changed complete source grid')
        require(len(pre['endpoint_upper']) == len(post['endpoint_upper']) == 513 and pre['endpoint_upper'][0] == post['endpoint_upper'][0] == 0, 'Changed zero initial radius or endpoint cover')
        require(pre['precision_bits'] >= 128 and post['precision_bits'] == pre['precision_bits'], 'Changed reference precision')
        if post['maximum_variance_refinement']['usable']:
            require(post['source_upper_values_preserved'] and post['center_lower_preserved'], 'Changed reference refinement scope')
        else:
            require(all((post[k] == pre[k] for k in ('E_upper', 'endpoint_upper', 'lower_offset', 'upper_offset', 'eligible'))), 'Unavailable maximum-variance refinement changed its inherited bounds')
        require(post['eligible'] == pre['eligible'] and post['lower_offset'] == pre['lower_offset'] and (len(post['E_upper']) == len(post['upper_offset']) == 512), 'Changed eligibility or lower center')
        for index, cell in enumerate(pre['cells']):
            require(cell['index'] == index and type(cell['eligible']) is bool, 'Changed source cell identity')
            require(Q(nodes[index][0]) <= Q(nodes[index][1]) <= Q(nodes[index + 1][0]) <= Q(nodes[index + 1][1]), 'Source node enclosures overlap or reverse')
            if not cell['eligible']:
                require(original[index] is None and factors[index] == '1', 'Changed ineligible packed source')
                continue
            counts['eligible_cells'] += 1
            old_source = Q(cell['old_source_upper'])
            energy = cell['energy_source_upper']
            prior = original[index]
            fresh = cell['packed_source_upper']
            require(old_source >= 0 and (energy is None or Q(energy) >= 0), 'Invalid old complete source')
            before = [old_source] + ([] if energy is None else [Q(energy)])
            after = list(before)
            if prior is None:
                require(fresh is None and factors[index] == '1', 'Bypassed the original packed-frequency guard')
            else:
                eta = source_factor(proof, Q(nodes[index][0]), mode)
                require(factors[index] == str(eta) and Q(0) <= Q(fresh) <= Q(prior) and (Q(fresh) >= eta * Q(prior)), 'Packed source reduced beyond its proved factor')
                if mode != 'contracted':
                    require(fresh == prior, 'Identity changed packed arithmetic')
                counts['packed_available_cells'] += 1
                counts['original_packed_minimizers'] += Q(prior) <= min(before)
                counts['strictly_reduced_packed_candidates'] += Q(fresh) < Q(prior)
                before.append(Q(prior))
                after.append(Q(fresh))
            require(Q(cell['source_upper']) == min(after), 'Changed the original source-minimum operation')
            counts['strictly_reduced_selected_sources'] += min(after) < min(before)
            require(Q(0) <= Q(post['E_upper'][index]) <= Q(pre['E_upper'][index]) and Q(0) <= Q(post['endpoint_upper'][index + 1]) <= Q(post['E_upper'][index]) and (Q(post['upper_offset'][index]) <= Q(pre['upper_offset'][index])), 'Invalid frozen attained-maximum recurrence/refinement output')
    return counts

def compare_receipts(baseline, identity, contracted):
    require({k: v for k, v in baseline.items() if k != 'mode'} == {k: v for k, v in identity.items() if k != 'mode'}, 'Identity changed an original source candidate or recurrence field')
    for old, new in zip(baseline['columns'], contracted['columns'], strict=True):
        require(old['node_enclosures'] == new['node_enclosures'] and old['deletion_caps'] == new['deletion_caps'] and (old['original_packed_candidates'] == new['original_packed_candidates']), 'Changed original frequency cells, deletion caps, or packed candidates')
        a, b = (old['packed_reference'], new['packed_reference'])
        for key in ('old_E_upper', 'old_endpoint_upper', 'lower_offset', 'upper_offset', 'eligible', 'precision_bits'):
            require(a[key] == b[key], 'Changed inherited elementary bound or source center')
        for ca, cb in zip(a['cells'], b['cells'], strict=True):
            require(ca['eligible'] == cb['eligible'], 'Changed source eligibility')
            if ca['eligible']:
                for key in ('old_source_upper', 'energy_source_upper', 'old_E_upper', 'propagation_ratio_upper', 'lower_offset', 'upper_offset'):
                    require(ca[key] == cb[key], 'Changed an unmodified candidate, center, or propagator')
        for key in ('lower_offset', 'upper_offset', 'eligible', 'maximum_variance_refinement', 'precision_bits'):
            require(old['refined_reference'][key] == new['refined_reference'][key], 'Changed final signed center, eligibility, or attained-maximum geometry')
    return True

def score(consumer, Llo, Lhi, actual_box, reference_box, envelopes, mode):
    require(mode in MODES, 'Invalid complete portfolio mode')
    proof = certificate(Llo, Lhi, actual_box)
    require(actual_box[:2] == reference_box[:2] and actual_box[4:] == reference_box[4:] and (reference_box[2] <= actual_box[2] <= actual_box[3] <= reference_box[3]), 'Require the unchanged original D/tau reference and actual b rectangle')
    import fast_maximum_variance_signed as binding
    from fast_packed_reference import upper, exact
    from flint import arb, ctx
    original_bounds = binding.fast_packed_energy_reference_bounds
    original_refine = binding.refine_reference
    require(sha(original_bounds.__code__.co_filename) == PACKED_SHA and sha(original_refine.__code__.co_filename) == REFINE_SHA, 'Changed active numerical recurrence functions')
    require(consumer._reference.__func__.__globals__ is vars(binding), 'The selected consumer does not use the inspected reference binding')

    def clear():
        for name in ('reference_cache', 'prepared_reference', 'score_cache'):
            getattr(consumer, name).clear()
            require(not getattr(consumer, name), 'Reference cache did not clear')
    pending = {}
    columns = []

    class PackedProxy:

        def __init__(self, original, geometry):
            self.original = original
            self.cells = {id(row[0]): (index, row, Q(geometry.parsed[index][0])) for index, row in enumerate(geometry.rows)}
            self.candidates = {}
            self.factors = {}

        def __getattr__(self, name):
            return getattr(self.original, name)

        def packing_source(self, frequency, frequency_square, deletion):
            index, row, left = self.cells[id(frequency)]
            require(frequency is row[0] and frequency_square is row[1] and (index not in self.candidates), 'The packed source bypassed its original complete frequency cell')
            old = self.original.packing_source(frequency, frequency_square, deletion)
            self.candidates[index] = old
            eta = Q(1) if old is None else source_factor(proof, left, mode)
            self.factors[index] = str(eta)
            if old is None or eta == 1:
                return old
            require(ctx.prec >= 128, 'Packed contraction used insufficient precision')
            value = min(old, upper(arb(old) * exact(eta)))
            require(Q(value) >= Q(old) * eta, 'Packed factor multiplication rounded inward')
            return value

    def bounds(nodes, deletion_caps, **kwargs):
        geometry = kwargs['prepared_geometry']
        box = kwargs['prepared_box']
        require(tuple(geometry.parsed) == tuple(box.base._nodes(nodes)), 'Changed prepared source geometry')
        proxy = None
        if mode != 'baseline':
            proxy = PackedProxy(box, geometry)
            kwargs = dict(kwargs, prepared_box=proxy)
        result = original_bounds(nodes, deletion_caps, **kwargs)
        originals = []
        factors = []
        for index, cell in enumerate(result['cells']):
            if not cell['eligible']:
                require(proxy is None or index not in proxy.candidates, 'An ineligible packed cell was evaluated')
                originals.append(None)
                factors.append('1')
            elif proxy is None:
                originals.append(cell['packed_source_upper'])
                factors.append('1')
            else:
                require(index in proxy.candidates, 'Omitted an eligible original packed candidate')
                originals.append(proxy.candidates[index])
                factors.append(proxy.factors[index])
        require(len(result['cells']) == 512, 'Changed fixed source mesh')
        pending[id(result)] = dict(node_enclosures=[[str(a), str(b)] for a, b in geometry.parsed], deletion_caps=encode(list(deletion_caps)), original_packed_candidates=encode(originals), source_factors=factors, packed_reference=encode(result))
        return result

    def refine(result, nodes, **kwargs):
        require(id(result) in pending, 'Reference refinement skipped its recorded original source')
        record = pending.pop(id(result))
        answer = original_refine(result, nodes, **kwargs)
        record['refined_reference'] = refined_fields(answer)
        columns.append(record)
        return answer
    clear()
    try:
        binding.fast_packed_energy_reference_bounds = bounds
        binding.refine_reference = refine
        answer = consumer.score(Llo, Lhi, reference_box, envelopes, detail=True)
        require(len(columns) == 3 and (not pending), 'Incomplete source-column or recurrence capture')
    finally:
        binding.fast_packed_energy_reference_bounds = original_bounds
        binding.refine_reference = original_refine
        clear()
    receipt = dict(certificate=proof, mode=mode, actual_box=actual_box, reference_box=reference_box, packed_function_sha256=PACKED_SHA, refine_function_sha256=REFINE_SHA, columns=columns, all_three_caches_cleared_before_and_after=True, original_bindings_restored=True)
    check_receipt(receipt, Llo, Lhi, actual_box, reference_box, mode)
    return (answer, receipt)
