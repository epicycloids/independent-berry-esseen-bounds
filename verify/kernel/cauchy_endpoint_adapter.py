from fractions import Fraction as Q
import hashlib
import json
from math import isqrt
from types import MethodType
RATIO = Q(27, 64)

def require(condition, message):
    if not condition:
        raise ValueError(message)

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()

def certificate(Llo, Lhi, box, eta):
    dlo, dhi, blo, bhi, taulo, tauhi = map(Q, box)
    lo, hi, eta = (Q(Llo), Q(Lhi), Q(eta))
    require(0 < lo <= hi and 0 <= dlo <= dhi <= 1 and (0 <= blo <= bhi) and (0 <= taulo <= tauhi), 'Invalid original closed moment rectangle')
    scale = 2 ** 64
    root_lo = isqrt((dlo.numerator << 128) // dlo.denominator)
    root_hi = isqrt((dhi.numerator << 128) // dhi.denominator)
    require(Q(root_lo ** 2, scale ** 2) <= dlo < Q((root_lo + 1) ** 2, scale ** 2) and Q(root_hi ** 2, scale ** 2) <= dhi < Q((root_hi + 1) ** 2, scale ** 2), 'Invalid two-sided rational square-root enclosures')
    wlo = dlo * Q(root_lo, scale)
    whi = dhi * Q(root_hi + 1, scale)
    excess = hi - taulo
    rest_margin = taulo - whi - Q(15, 23) * (hi - blo)
    delta = hi * wlo - bhi * taulo
    denominator = excess * (taulo + RATIO * excess)
    require(0 < eta <= 1 and denominator > 0, 'Invalid reference contraction')
    if eta < 1:
        require(rest_margin >= 0 and 0 < delta <= denominator and (eta ** 2 + (delta / denominator) ** 2 >= 1), 'No proved monotone-endpoint contraction on this actual rectangle')
    return dict(kind='monotone_endpoint_transport.v1', L=[str(lo), str(hi)], actual_box=list(map(str, map(Q, box))), eta=str(eta), sqrt_d_lower=str(Q(root_lo, scale)), sqrt_d_upper=str(Q(root_hi + 1, scale)), d_three_halves_lower=str(wlo), d_three_halves_upper=str(whi), transported_L=str(hi), transported_tau=str(taulo), packed_transport_ratio='15/23', rest_transport_margin=str(rest_margin), original_whole_box_defect_lower=str(lo * wlo - bhi * tauhi), excess_upper=str(excess), defect_lower=str(delta), denominator_upper=str(denominator), price_over_Q=str(RATIO), contraction_margin=str(eta ** 2 + (delta / denominator) ** 2 - 1) if eta < 1 else None)

def contracted_radii(E, eligible, proof):
    eta = _upper(Q(proof['eta']))
    require(Q(eta) >= Q(proof['eta']) and eta <= 1, 'Contraction coefficient was rounded inward')
    contracted = np.minimum(E, up(E * eta))
    require(np.all(np.isfinite(E)) and np.all(E >= 0) and np.all(contracted <= E), 'Invalid finite eligible reference radius')
    record = dict(certificate=proof, eta_float=eta, eligible=eligible.tolist(), old_radius=E.tolist(), contracted_radius=contracted.tolist(), strict_radius_improvements=int(np.count_nonzero(contracted < E)), first_native_cell_unchanged=True, reference_center_unchanged=True)
    return (contracted, record)

def proposal(self, Llo, Lhi, dlo, dhi, taulo, tauhi):
    require((Llo, Lhi, dlo, dhi, taulo, tauhi) == self._cauchy_geometry, 'Changed containing reference geometry')
    w = self.weights
    E, offset_lo, offset_hi, eligible = self._reference(Llo, Lhi, dlo, dhi, taulo, tauhi)
    eligible = eligible.copy()
    eligible[0, :] = False
    E = np.where(eligible, E, 0.0)
    offset_lo = np.where(eligible, offset_lo, 0.0)
    offset_hi = np.where(eligible, offset_hi, 0.0)
    E, receipt = contracted_radii(E, eligible, self._cauchy_certificate)
    product = np.maximum.reduce([up(self.qlo * offset_lo), up(self.qlo * offset_hi), up(self.qhi * offset_lo), up(self.qhi * offset_hi)])
    center = up(product * np.where(product >= 0, self.width_twice, self.width_low))
    numerator = up(up(E * w.low) + center)
    proposed = normalize(numerator, Llo, Lhi)
    proposed.flags.writeable = False
    eligible.flags.writeable = False
    self._cauchy_receipts.append(receipt)
    return (proposed, eligible)

def score(consumer, Llo, Lhi, actual_box, reference_box, envelopes, eta):
    global np, up, normalize, _upper
    import numpy as np
    from interval_signed import up, normalize, _upper
    import stable_bounds
    require(all((reference_box[2 * j] <= actual_box[2 * j] <= actual_box[2 * j + 1] <= reference_box[2 * j + 1] for j in range(3))), 'The actual scoring rectangle is outside the original reference box')
    proof = certificate(Llo, Lhi, actual_box, eta)
    clipped = stable_bounds.feasible_clip(reference_box, Lhi)
    require(clipped is not None, 'No valid containing reference rectangle')
    dlo, dhi, _, _, taulo, tauhi = clipped
    names = ('_proposal', '_cauchy_certificate', '_cauchy_geometry', '_cauchy_receipts')
    sentinel = object()
    previous = {name: consumer.__dict__.get(name, sentinel) for name in names}
    try:
        consumer._cauchy_certificate = proof
        consumer._cauchy_geometry = (Llo, Lhi, dlo, dhi, taulo, tauhi)
        consumer._cauchy_receipts = []
        consumer._proposal = MethodType(proposal, consumer)
        answer = consumer.score(Llo, Lhi, reference_box, envelopes, detail=True)
        require(len(consumer._cauchy_receipts) == 1, 'The complete signed score bypassed the reference contraction')
        return (answer, consumer._cauchy_receipts[0])
    finally:
        for name, value in previous.items():
            if value is sentinel:
                consumer.__dict__.pop(name, None)
            else:
                setattr(consumer, name, value)

def check_receipt(record, Llo, Lhi, box, eta):
    proof = certificate(Llo, Lhi, box, eta)
    require(record['certificate'] == proof and Q(record['eta_float']) >= Q(eta) and (record['eta_float'] <= 1) and record['first_native_cell_unchanged'] and record['reference_center_unchanged'], 'Changed uniform contraction proof or scope')
    old, new, eligible = (record[k] for k in ('old_radius', 'contracted_radius', 'eligible'))
    require(len(old) == len(new) == len(eligible) == 512 and all((len(v) == 3 for a in (old, new, eligible) for v in a)), 'Changed native reference grid')
    count = 0
    for j, (a, b, flags) in enumerate(zip(old, new, eligible, strict=True)):
        for x, y, flag in zip(a, b, flags, strict=True):
            require(type(flag) is bool and Q(0) <= Q(y) <= Q(x) and (Q(y) >= Q(eta) * Q(x)), 'A radius was reduced beyond the exact contraction')
            require(flag or x == y == 0, 'An ineligible radius entered the modified signed score')
            if j == 0:
                require(not flag and x == y == 0, 'The excluded first native cell changed')
            count += y < x
    require(count == record['strict_radius_improvements'], 'Changed radius-improvement count')
    if Q(eta) == 1:
        require(old == new, 'Identity contraction changed arithmetic')
    return proof
