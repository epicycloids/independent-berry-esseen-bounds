from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
PACKED_SHA = '5e2de769e7b13f28370974ad0ff77cf2672e322fa80dcd8386549a4441c1526e'
REFINE_SHA = '0fbd1217d20b2447f50039ccd7af8ba554a981ed27277ac77431bd0abe67054f'
MODES = ('baseline', 'identity', 'polynomial', 'scalar_price1', 'scalar_price1_individual')
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

def cube_root(n):
    assert type(n) is int and n >= 0
    lo = 0
    hi = 1 << (n.bit_length() + 2) // 3
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** 3 <= n:
            lo = mid
        else:
            hi = mid
    if hi ** 3 <= n:
        lo = hi
    assert lo ** 3 <= n < (lo + 1) ** 3
    return lo

def two_thirds_bounds(value):
    value = Q(value)
    assert value >= 0
    scale = 2 ** 64
    square = value * value
    k = cube_root(square.numerator * scale ** 3 // square.denominator)
    lower, upper = (Q(k, scale), Q(k + 1, scale))
    assert lower ** 3 <= square < upper ** 3
    return (lower, upper)

def general_certificate(Llo, Lhi, box):
    dl, dh, bl, bh, a, z = map(Q, box)
    lo, M = (Q(Llo), Q(Lhi))
    assert 0 < lo <= M and 0 < dl <= dh <= 1 and (0 <= bl <= bh) and (0 < a <= z)
    E = M - a
    assert E > 0
    scale = 2 ** 64
    kroot = math.isqrt(dh.numerator * scale ** 2 // dh.denominator)
    assert Q(kroot ** 2, scale ** 2) <= dh < Q((kroot + 1) ** 2, scale ** 2)
    W0 = dh * Q(kroot, scale)
    Wup = dh * Q(kroot + 1, scale)
    balance = bh * a / M
    margins = dict(balance_gap=str(W0 - balance), transported_distinguished_excess_lower=str(bl - Wup), rest_transport_margin=str(a - Wup - Q(15, 23) * (M - bl)))
    failure = []
    if not 0 < balance < W0:
        failure.append('balance_not_strictly_below_weight_endpoint')
    if bl < Wup:
        failure.append('transported_distinguished_excess_can_be_negative')
    if Q(margins['rest_transport_margin']) < 0:
        failure.append('packed_rest_transport_condition_fails')
    if failure:
        return dict(covers=False, failures=failure, exact_margins=margins)
    k = a // W0
    if k < 1:
        return dict(covers=False, failures=['no_positive_packing_count'], exact_margins=margins)
    lower = max(balance, a / (k + 1))
    c = (lower + W0) / 2
    assert 0 < balance < c < W0 and a // c == a // W0 == k
    residual = a - k * c
    assert 0 <= residual < c
    d0, _ = two_thirds_bounds(c)
    _, gamma = two_thirds_bounds(residual)
    assert d0 > gamma and d0 ** 3 <= c * c and (gamma ** 3 >= residual * residual)
    C = Q(40 * k, 27) * (d0 - gamma)
    A = Q(29, 4) * a + Q(9, 8) * E
    K = C * (W0 - c) / A
    delta = M * c - bh * a
    D = E * (a + Q(27, 64) * E)
    assert 0 < delta <= D and K > 0
    H = (delta / D) ** 2
    eta = 1 - min(K, H) / 2
    assert 0 < eta < 1 and eta * eta + min(K, H) >= 1
    return dict(covers=True, failures=[], exact_margins=margins, certificate_kind='packed_weight_general_witness.v1', L=[str(lo), str(M)], actual_box=list(map(str, map(Q, box))), weight_endpoint_lower=str(W0), weight_endpoint_upper=str(Wup), balance_weight=str(balance), analytic_cut=str(c), lower_packing_join=str(a / (k + 1)), packing_count=int(k), residual_upper=str(residual), power_lower_bound=str(d0), residual_power_upper_bound=str(gamma), low_weight_power_margin=str(c * c - d0 ** 3), residual_power_margin=str(gamma ** 3 - residual ** 2), derivative_gap_coefficient=str(C), denominator_upper=str(A), frequency_squared_loss_coefficient=str(K), positive_defect=str(delta), Cauchy_denominator=str(D), high_weight_loss=str(H), rational_source_factor_at_one=str(eta), original_packed_guard_required=True, original_radius_recurrence_required=True, signed_center_unchanged=True)

def certificate(Llo, Lhi, box):
    value = general_certificate(Llo, Lhi, box)
    require(value['covers'], 'The complete regional rectangle fails the original packed prerequisites')
    M = Q(Lhi)
    a = Q(box[4])
    bh = Q(box[3])
    W0 = Q(value['weight_endpoint_lower'])
    k = value['packing_count']
    lower = max(bh * a / M, a / (k + 1))
    A = Q(value['denominator_upper'])
    D = Q(value['Cauchy_denominator'])
    cuts = []
    for part in (Q(1, 4), Q(1, 2), Q(3, 4)):
        c = lower + part * (W0 - lower)
        r = a - k * c
        d0, _ = two_thirds_bounds(c)
        _, gamma = two_thirds_bounds(r)
        require(a // c == a // W0 == k and d0 > gamma >= 0 and (d0 ** 3 <= c * c) and (gamma ** 3 >= r * r), 'Invalid complete fixed packing interval')
        delta = M * c - bh * a
        require(0 < delta <= D, 'Invalid Cauchy loss')
        numerator = [k * (W0 - c) * q * (1 + Q(2 * n, 3)) * (d0 ** n - gamma ** n) for n, q in ((1, Q(8, 9)), (2, Q(8, 15)), (3, Q(736, 2835)))]
        cuts.append(dict(fraction=str(part), cut=str(c), packing_count=k, d0=str(d0), gamma=str(gamma), polynomial_numerators=list(map(str, numerator)), high_weight_loss=str((delta / D) ** 2)))
    require(Q(cuts[1]['polynomial_numerators'][0]) / A == Q(value['frequency_squared_loss_coefficient']) and cuts[1]['high_weight_loss'] == value['high_weight_loss'], 'Changed midpoint first-term certificate')
    e = max(Q(0), Q(box[2]) - Q(value['weight_endpoint_upper']))
    excess = max(Q(0), Q(Llo) - Q(box[5]))
    dlo = Q(box[0])
    scale = 2 ** 64
    k = math.isqrt(dlo.numerator * scale ** 2 // dlo.denominator)
    root_lower = Q(k, scale)
    require(root_lower ** 2 <= dlo, 'Invalid lower distinguished square root')
    weight_lower = dlo * root_lower
    excess_upper = max(Q(0), Q(box[3]) - weight_lower)
    return dict(value, kind='packed_polynomial_scalar_price1_individual.v1', packed_source_only=True, final_radius_scaling=False, reference_center_unchanged=True, fixed_cuts=cuts, denominator_rule='min(global_A,a*q_upper_at_outward_cell_argument+9E/8)', original_actual_and_padded_frequency_guards_required=True, scalar_price1=dict(normalized_frequency=['3/4', '1'], distinguished_excess_lower=str(e), total_excess_lower=str(excess), unmodified_square_saving=str(excess * e / 8), distinguished_weight_lower=str(weight_lower), distinguished_excess_upper=str(excess_upper), individual_rule='R-K*e^2/(8*(U+U_prime))', linear_rule='R-K^2*e^2/(16*R)', square_rule='sqrt(R_original^2-K^2*E*e/8)', distinguished_coordinate_preserved=True, actual_scoring_box_used=True))

def upward_rational(value):
    value = Q(value)
    require(value >= 0, 'Require a nonnegative exact source')
    result = 0.0 if value == 0 else math.nextafter(float(value), math.inf)
    require(math.isfinite(result) and Q(result) >= value, 'Inward rational source conversion')
    return result

def upward_sqrt(value):
    value = Q(value)
    require(value >= 0, 'Require a nonnegative exact source square')
    if not value:
        return 0.0
    scale = 2 ** 128
    k = math.isqrt(value.numerator * scale ** 2 // value.denominator)
    result = upward_rational(Q(k + 1, scale))
    require(Q(result) ** 2 >= value, 'Inward source square root')
    return result

def price1_candidates(proof, box, left, right, deletion, original, polynomial, mode, q_upper=None):
    require(mode in MODES and Q(0) <= Q(left) <= Q(right) and (Q(deletion) >= 0), 'Invalid complete frequency cell or deletion multiplier')
    K = Q(right) ** 2 * Q(deletion) / 2
    e = Q(proof['scalar_price1']['distinguished_excess_lower'])
    active = mode in ('scalar_price1', 'scalar_price1_individual') and original is not None and (e > 0) and (Q(left) ** 2 * Q(box[0]) >= Q(9, 16)) and (Q(right) ** 2 * Q(box[1]) <= 1)
    linear = square = individual = old_coordinate = new_coordinate = denominator = None
    if active:
        require(Q(0) <= Q(polynomial) <= Q(original), 'Invalid packed candidate ordering')
        if not K:
            linear = square = 0.0
        elif polynomial > 0:
            target = Q(polynomial) - K * K * e * e / (16 * Q(polynomial))
            if target >= 0:
                linear = min(polynomial, upward_rational(target))
        radicand = Q(original) ** 2 - K * K * Q(proof['scalar_price1']['unmodified_square_saving'])
        if radicand >= 0:
            square = min(original, upward_sqrt(radicand))
        if mode == 'scalar_price1_individual' and q_upper is not None:
            q = Q(q_upper)
            require(q >= Q(8, 3), 'Invalid individual q upper')
            emax = Q(proof['scalar_price1']['distinguished_excess_upper'])
            wmax = Q(proof['weight_endpoint_upper'])
            require(emax >= e > 0 and wmax > 0, 'Invalid actual distinguished coordinate bounds')
            old_coordinate = upward_sqrt(emax * (q * wmax + Q(9, 8) * emax))
            new_coordinate = upward_sqrt(emax * (q * wmax + emax))
            denominator = 8 * (Q(old_coordinate) + Q(new_coordinate))
            require(denominator > 0, 'Nonpositive individual rationalization denominator')
            target = Q(polynomial) - K * e * e / denominator
            if target >= 0:
                individual = min(polynomial, upward_rational(target))
    fresh = polynomial
    if fresh is not None:
        fresh = min([fresh] + [v for v in (linear, square, individual) if v is not None])
    return dict(active=bool(active), multiplier=str(K), polynomial=polynomial, linear=linear, unmodified_square=square, individual=individual, individual_old_upper=old_coordinate, individual_new_upper=new_coordinate, individual_denominator=None if denominator is None else str(denominator), selected=fresh)

def source_factor(proof, left, mode, denominator=None):
    require(mode in MODES and Q(left) >= 0, 'Invalid source mode or left frequency')
    if mode in ('baseline', 'identity'):
        return Q(1)
    else:
        A = Q(proof['denominator_upper']) if denominator is None else Q(denominator)
        require(0 < A <= Q(proof['denominator_upper']), 'Invalid cellwise denominator upper bound')
        loss = max((min(sum((Q(v) * Q(left) ** (2 * n) for n, v in enumerate(cut['polynomial_numerators'], 1))) / A, Q(cut['high_weight_loss'])) for cut in proof['fixed_cuts']))
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
    counts = dict(eligible_cells=0, packed_available_cells=0, original_packed_minimizers=0, strictly_reduced_packed_candidates=0, strictly_reduced_selected_sources=0, price1_guarded_cells=0, price1_strictly_reduced_polynomial_candidates=0, price1_strictly_reduced_selected_sources=0)
    for column in columns:
        pre = column['packed_reference']
        post = column['refined_reference']
        nodes = column['node_enclosures']
        original = column['original_packed_candidates']
        factors = column['source_factors']
        require(len(nodes) == 513 and len(original) == len(factors) == len(column['deletion_caps']) == len(pre['cells']) == 512, 'Changed complete source grid')
        require(all((len(column[key]) == 512 for key in ('cell_denominators', 'q_argument_uppers', 'q_value_uppers', 'cap_root_uppers'))), 'Missing complete cell denominator evidence')
        require(len(column['price1_transfers']) == 512, 'Missing complete scalar transfer evidence')
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
                require(column['price1_transfers'][index] is None, 'Applied price one to an ineligible cell')
                require(column['cell_denominators'][index] == proof['denominator_upper'] and column['q_argument_uppers'][index] is None and (column['q_value_uppers'][index] is None) and (column['cap_root_uppers'][index] is None), 'Changed ineligible denominator evidence')
                continue
            counts['eligible_cells'] += 1
            old_source = Q(cell['old_source_upper'])
            energy = cell['energy_source_upper']
            prior = original[index]
            fresh = cell['packed_source_upper']
            denominator = Q(column['cell_denominators'][index])
            argument = column['q_argument_uppers'][index]
            qvalue = column['q_value_uppers'][index]
            caproot = column['cap_root_uppers'][index]
            if mode not in ('polynomial', 'scalar_price1', 'scalar_price1_individual') or prior is None or argument is None:
                require(argument is None and qvalue is None and (caproot is None) and (denominator == Q(proof['denominator_upper'])), 'Changed unmodified or guarded denominator')
            else:
                require(Q(0) <= Q(argument) <= Q(10001, 10000) and Q(qvalue) >= Q(8, 3), 'Invalid outward q argument or lower q constraint')
                require(Q(caproot) > 0 and Q(caproot) ** 2 >= Q(box[1]) and (Q(argument) >= Q(nodes[index + 1][1]) * Q(caproot)), 'Recorded q argument is below the packed cap argument')
                a = Q(box[4])
                E = Q(Lhi) - a
                require(denominator == min(Q(proof['denominator_upper']), a * Q(qvalue) + Q(9, 8) * E), 'Changed exact cellwise denominator')
            require(old_source >= 0 and (energy is None or Q(energy) >= 0), 'Invalid old complete source')
            before = [old_source] + ([] if energy is None else [Q(energy)])
            after = list(before)
            transfer = column['price1_transfers'][index]
            polynomial = transfer['polynomial']
            expected_transfer = price1_candidates(proof, box, Q(nodes[index][0]), Q(nodes[index + 1][1]), Q(column['deletion_caps'][index]), prior, polynomial, mode, qvalue)
            require(transfer == expected_transfer and fresh == transfer['selected'], 'Changed exact scalar guard, forcing multiplier, or compatible source subtraction')
            if prior is None:
                require(fresh is polynomial is None and factors[index] == '1', 'Bypassed the original packed-frequency guard')
            else:
                eta = source_factor(proof, Q(nodes[index][0]), mode, denominator)
                require(factors[index] == str(eta) and Q(0) <= Q(fresh) <= Q(polynomial) <= Q(prior) and (Q(polynomial) >= eta * Q(prior)), 'Packed polynomial source reduced beyond its proved factor')
                if mode in ('baseline', 'identity'):
                    require(fresh == prior, 'Identity changed packed arithmetic')
                counts['packed_available_cells'] += 1
                counts['original_packed_minimizers'] += Q(prior) <= min(before)
                counts['strictly_reduced_packed_candidates'] += Q(fresh) < Q(prior)
                before.append(Q(prior))
                after.append(Q(fresh))
                counts['price1_guarded_cells'] += transfer['active']
                counts['price1_strictly_reduced_polynomial_candidates'] += Q(fresh) < Q(polynomial)
                counts['price1_strictly_reduced_selected_sources'] += min(after) < min(before[:-1] + [Q(polynomial)])
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

        def __init__(self, original, geometry, deletion_caps):
            self.original = original
            self.cells = {id(row[0]): (index, row, Q(geometry.parsed[index][0])) for index, row in enumerate(geometry.rows)}
            self.candidates = {}
            self.factors = {}
            self.denominators = {}
            self.arguments = {}
            self.qvalues = {}
            self.caproots = {}
            self.geometry = geometry
            self.deletion_caps = deletion_caps
            self.transfers = {}

        def __getattr__(self, name):
            return getattr(self.original, name)

        def packing_source(self, frequency, frequency_square, deletion):
            index, row, left = self.cells[id(frequency)]
            require(frequency is row[0] and frequency_square is row[1] and (index not in self.candidates), 'The packed source bypassed its original complete frequency cell')
            old = self.original.packing_source(frequency, frequency_square, deletion)
            self.candidates[index] = old
            A = Q(proof['denominator_upper'])
            argument = qvalue = caproot = None
            if old is not None and mode in ('polynomial', 'scalar_price1', 'scalar_price1_individual'):
                packing = self.original.packing
                require(packing is not None, 'Missing the original prepared packed moment')
                caproot = upper(packing.cap_root)
                argument = upper(frequency * exact(Q(caproot)))
                if Q(argument) > Q(10001, 10000):
                    argument = caproot = None
                else:
                    qvalue = upper(packing.q_ball(exact(Q(argument))))
                    require(Q(qvalue) >= Q(8, 3), 'Invalid outward packed q value')
                    a = Q(actual_box[4])
                    E = Q(Lhi) - a
                    A = min(A, a * Q(qvalue) + Q(9, 8) * E)
            self.denominators[index] = str(A)
            self.arguments[index] = argument
            self.qvalues[index] = qvalue
            self.caproots[index] = caproot
            eta = Q(1) if old is None else source_factor(proof, left, mode, A)
            self.factors[index] = str(eta)
            value = old
            if old is not None and eta != 1:
                require(ctx.prec >= 128, 'Packed contraction used insufficient precision')
                value = min(old, upper(arb(old) * exact(eta)))
                require(Q(value) >= Q(old) * eta, 'Packed factor multiplication rounded inward')
            transfer = price1_candidates(proof, actual_box, left, Q(self.geometry.parsed[index + 1][1]), Q(self.deletion_caps[index]), old, value, mode, qvalue)
            self.transfers[index] = transfer
            return transfer['selected']

    def bounds(nodes, deletion_caps, **kwargs):
        geometry = kwargs['prepared_geometry']
        box = kwargs['prepared_box']
        require(tuple(geometry.parsed) == tuple(box.base._nodes(nodes)), 'Changed prepared source geometry')
        proxy = None
        if mode != 'baseline':
            proxy = PackedProxy(box, geometry, deletion_caps)
            kwargs = dict(kwargs, prepared_box=proxy)
        result = original_bounds(nodes, deletion_caps, **kwargs)
        originals = []
        factors = []
        denominators = []
        arguments = []
        qvalues = []
        caproots = []
        transfers = []
        for index, cell in enumerate(result['cells']):
            if not cell['eligible']:
                require(proxy is None or index not in proxy.candidates, 'An ineligible packed cell was evaluated')
                originals.append(None)
                factors.append('1')
                transfers.append(None)
            elif proxy is None:
                originals.append(cell['packed_source_upper'])
                factors.append('1')
                transfers.append(price1_candidates(proof, actual_box, Q(geometry.parsed[index][0]), Q(geometry.parsed[index + 1][1]), Q(deletion_caps[index]), cell['packed_source_upper'], cell['packed_source_upper'], mode))
            else:
                require(index in proxy.candidates, 'Omitted an eligible original packed candidate')
                originals.append(proxy.candidates[index])
                factors.append(proxy.factors[index])
                transfers.append(proxy.transfers[index])
            denominators.append(str(proof['denominator_upper']) if proxy is None or index not in proxy.denominators else proxy.denominators[index])
            arguments.append(None if proxy is None or index not in proxy.arguments else proxy.arguments[index])
            qvalues.append(None if proxy is None or index not in proxy.qvalues else proxy.qvalues[index])
            caproots.append(None if proxy is None or index not in proxy.caproots else proxy.caproots[index])
        require(len(result['cells']) == 512, 'Changed fixed source mesh')
        pending[id(result)] = dict(node_enclosures=[[str(a), str(b)] for a, b in geometry.parsed], deletion_caps=encode(list(deletion_caps)), original_packed_candidates=encode(originals), source_factors=factors, cell_denominators=denominators, q_argument_uppers=arguments, q_value_uppers=qvalues, cap_root_uppers=caproots, price1_transfers=encode(transfers), packed_reference=encode(result))
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
