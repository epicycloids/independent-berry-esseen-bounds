from functools import lru_cache
from pathlib import Path
import math
import sys
import numpy as np
from flint import arb, ctx
HERE = Path(__file__).resolve().parent
R27 = HERE.parent
for path in (R27.parent / '26', R27, R27 / 'dual_logarithm'):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
from interval_bounds import K, al, au
from paired_bounds import cos_abs_upper
from stable_bounds import feasible_clip
from dual_log_bounds import support_table
PRICES = (0.75, 1.0, 1.5, 2.0, 3.0)
MODES = ('simple', 'dual', 'dual_cubic', 'dual_tyurin')

def _lo(value):
    if not value.is_finite():
        raise ArithmeticError('A finite Arb enclosure is required')
    return 0.0 if value == 0 else al(value)

def _hi(value):
    if not value.is_finite():
        raise ArithmeticError('A finite Arb enclosure is required')
    return 0.0 if value == 0 else au(value)

@lru_cache(maxsize=1)
def tyurin_supports():
    four = arb(4)
    n4 = (four * four - 6) * four.cos() - 4 * four * four.sin() + 6
    assert 4 < 3 * arb.pi() / 2 < 2 * arb.pi() and n4 > 0
    rows = [(0.0, 0.5, _hi(K))]
    for j in range(129, 201):
        z = arb(j) / 32
        assert 4 < z < 2 * arb.pi()
        f = (1 - z.cos()) / (z * z)
        derivative = z.sin() / (z * z) - 2 * (1 - z.cos()) / z ** 3
        alpha, beta = (f - z * derivative, -derivative)
        assert beta >= 0 and beta <= K
        assert alpha - 4 * beta <= arb(1) / 2 - 4 * K
        a, b = (max(0.0, _lo(alpha)), max(0.0, _hi(beta)))
        assert arb(a) <= alpha and arb(b) >= beta
        rows.append((j / 32, a, b))
    return tuple(rows)

def _roots(a, b, c):
    if a == 0:
        if b == 0:
            return []
        if b.contains(0):
            return [None]
        return [-c / b]
    if a.contains(0):
        return [None]
    discriminant = b * b - 4 * a * c
    if discriminant < 0:
        return []
    radical = discriminant.sqrt() if discriminant > 0 else arb(0).union(discriminant.upper().sqrt())
    return [(-b - radical) / (2 * a), (-b + radical) / (2 * a)]

def _interval(root, lo, hi):
    if root is None:
        return (lo, hi)
    if not root.is_finite():
        return (lo, hi)
    a, b = (max(lo, _lo(root)), min(hi, _hi(root)))
    return None if a > b else (a, b)

def _active_lines(lines, elo, ehi):
    active = []
    for line in lines:
        h, r = map(arb, line)
        redundant = False
        for old in active:
            a, b = map(arb, old)
            if h - r * arb(elo) >= a - b * arb(elo) and h - r * arb(ehi) >= a - b * arb(ehi):
                redundant = True
                break
        if redundant:
            continue
        active = [old for old in active if not (arb(old[0]) - arb(old[1]) * arb(elo) >= h - r * arb(elo) and arb(old[0]) - arb(old[1]) * arb(ehi) >= h - r * arb(ehi))]
        active.append(line)
    return active

def quadratic_exponential_max(polynomials, lines, elo, ehi, *, details=False):
    if not 0 <= elo <= ehi or not polynomials or (not lines):
        raise ValueError('A nonempty nonnegative excess interval is required')
    if not all((math.isfinite(x) and x >= 0 for p in polynomials for x in p)):
        raise ValueError('Finite nonnegative quadratic coefficients required')
    if not all((math.isfinite(h) and math.isfinite(r) and (r >= 0) for h, r in lines)):
        raise ValueError('Finite affine caps with nonnegative rates required')
    lines = _active_lines(lines, elo, ehi)
    polys = [tuple(map(arb, p)) for p in polynomials]
    exponentials = [tuple(map(arb, line)) for line in lines]
    candidates = [(elo, elo), (ehi, ehi)]

    def add_roots(a, b, c):
        for root in _roots(a, b, c):
            candidate = _interval(root, elo, ehi)
            if candidate is not None:
                candidates.append(candidate)
    for i, (a, b, c) in enumerate(polys):
        add_roots(c, b, a - 1)
        for aa, bb, cc in polys[:i]:
            add_roots(c - cc, b - bb, a - aa)
        for _, r in exponentials:
            add_roots(-2 * r * c, 2 * c - 2 * r * b, b - 2 * r * a)
    for j, (h, r) in enumerate(exponentials):
        add_roots(arb(0), -r, h)
        for a, b in exponentials[:j]:
            add_roots(arb(0), b - r, h - a)
    best = 0.0
    for left, right in candidates:
        x, y = (arb(left), arb(right))
        single = min(1.0, *(_hi((a + y * (b + y * c)).sqrt()) for a, b, c in polys))
        rest = min(1.0, *(_hi((h - r * x).exp()) for h, r in exponentials))
        best = max(best, min(1.0, _hi(arb(single) * arb(rest))))
    if details:
        return (best, dict(candidate_count=len(candidates), active_lines=len(lines), polynomials=len(polys)))
    return best

def _projection(Llo, Lhi, box):
    clipped = feasible_clip(box, Lhi)
    if clipped is None:
        return None
    dl, dh, bl, bh, tl, th = clipped
    if not 0 < Llo <= Lhi or not 0 <= dl <= dh <= 1:
        raise ValueError('Require positive L and a normalized variance box')
    wlo = max(0.0, _lo(arb(dl) * arb(dl).sqrt()))
    whi = max(0.0, _hi(arb(dh) * arb(dh).sqrt()))
    elo = max(0.0, _lo(arb(bl) - arb(whi)))
    ehi = min(_hi(arb(bh) - arb(wlo)), _hi(arb(Lhi) - arb(tl)))
    if ehi < elo:
        return None
    vlo, vhi = (max(0.0, _lo(1 - arb(dh))), max(0.0, _hi(1 - arb(dl))))
    srlo, srhi = (max(0.0, _lo(arb(tl) - arb(whi))), max(0.0, _hi(arb(th) - arb(wlo))))
    d2 = min(dh, vhi, _hi(arb(srhi) ** (arb(2) / 3)))
    return dict(box=clipped, dl=dl, dh=dh, tl=tl, th=th, wlo=wlo, whi=whi, elo=elo, ehi=ehi, vlo=vlo, vhi=vhi, srlo=srlo, srhi=srhi, d2=d2, Lhi=Lhi)

def _selected_tyurin(arguments):
    rows = tyurin_supports()
    selected = {0}
    for z in arguments:
        if math.isfinite(z) and z >= 4:
            selected.add(min(len(rows) - 1, max(1, round(32 * z) - 128)))
    return [rows[j] for j in sorted(selected)]

def _cell_portfolio(tlo, thi, C, g, mode):
    lo, hi = (arb(tlo), arb(thi))
    t2lo, t2hi, t3lo, t3hi = (lo * lo, hi * hi, lo ** 3, hi ** 3)
    wlo, whi, dl, dh = map(arb, (g['wlo'], g['whi'], g['dl'], g['dh']))
    t6 = t3hi * t3hi
    polys = [(_hi(arb(C) ** 2), _hi(arb(C) * t3hi / 3 + 2 * whi * t6 / 27), _hi(t6 / 18))]
    k = arb(_hi(K * t3hi))
    sigma = arb(max(0.0, _lo(arb(g['vlo']) * t2lo / 2)))
    h = -sigma + k * (arb(g['Lhi']) + arb(g['th']) - 2 * wlo)
    lines = [(_hi(h), float(k))]
    if mode in ('dual_cubic', 'dual_tyurin'):
        polys.append((max(0.0, _hi(1 - dl * t2lo + 4 * k * whi)), _hi(2 * k), 0.0))
    if mode == 'dual_tyurin':
        dmid = (g['dl'] + g['dh']) / 2
        vmid = (g['vlo'] + g['vhi']) / 2
        wmid = (g['wlo'] + g['whi']) / 2
        emid = (g['elo'] + g['ehi']) / 2
        es = (g['elo'], emid, g['ehi'])
        single_args = [thi * (2 * wmid + e) / dmid for e in es] if dmid > 0 else []
        rest_args = [thi * (g['Lhi'] + g['th'] - 2 * wmid - e) / vmid for e in es] if vmid > 0 else []
        for z, a, b in _selected_tyurin(single_args):
            if z == 0:
                continue
            A, B = (arb(a), arb(b))
            kk = arb(_hi(B * t3hi))
            polys.append((max(0.0, _hi(1 - 2 * A * dl * t2lo + 4 * kk * whi)), _hi(2 * kk), 0.0))
        for z, a, b in _selected_tyurin(rest_args):
            if z == 0:
                continue
            A, B = (arb(a), arb(b))
            kk = arb(_hi(B * t3hi))
            hh = -A * arb(g['vlo']) * t2lo + kk * (arb(g['Lhi']) + arb(g['th']) - 2 * wlo)
            lines.append((_hi(hh), float(kk)))
    if mode == 'simple' or g['vlo'] <= 0 or g['d2'] <= 0:
        return (polys, lines)
    argument = _hi(4 * K * hi * arb(g['d2']).sqrt())
    guide = 4 * float(K) * tlo * g['srlo'] / g['vhi'] if g['vhi'] > 0 else 0.0
    for price, u, aa, bb, domain, *_ in support_table():
        if price not in PRICES:
            continue
        wanted = np.searchsorted(u, guide, side='right') - 1
        permitted = np.searchsorted(-domain, -argument, side='right') - 1
        idx = min(np.clip(wanted, 0, len(u) - 1), max(0, permitted))
        if argument > domain[idx]:
            continue
        positive = arb(max(0.0, _lo(arb(float(aa[idx])) * t3lo)))
        negative = arb(_hi(arb(float(bb[idx])) * t2hi))
        a = arb(_lo(positive + (arb(price) - 1) * k))
        r = arb(_hi((1 + arb(price)) * k))
        tau = arb(g['tl'] if a >= 0 else g['th'])
        w = whi if a - r >= 0 else wlo
        hh = -sigma + negative * arb(g['vhi']) + r * arb(g['Lhi']) - a * tau + (a - r) * w
        lines.append((_hi(hh), float(r)))
    return (polys, lines)

def coupled_cf_cells(tlo, thi, Llo, Lhi, box, *, mode='dual_tyurin', return_details=False):
    if mode not in MODES:
        raise ValueError('Unknown coupled CF mode')
    tlo, thi = np.broadcast_arrays(np.asarray(tlo, dtype=float), np.asarray(thi, dtype=float))
    if not np.all(np.isfinite(tlo) & np.isfinite(thi) & (0 <= tlo) & (tlo <= thi)):
        raise ValueError('Finite ordered nonnegative frequency cells required')
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128)
        g = _projection(Llo, Lhi, box)
        if g is None:
            return None
        xlo = np.maximum(0, np.nextafter(tlo * _lo(arb(g['dl']).sqrt()), -np.inf))
        xhi = np.nextafter(thi * _hi(arb(g['dh']).sqrt()), np.inf)
        C = cos_abs_upper(xlo, xhi)
        result = np.empty(tlo.shape)
        candidates = active = polys_total = 0
        for index in np.ndindex(tlo.shape):
            polys, lines = _cell_portfolio(float(tlo[index]), float(thi[index]), float(C[index]), g, mode)
            value, metadata = quadratic_exponential_max(polys, lines, g['elo'], g['ehi'], details=True)
            result[index] = value
            candidates += metadata['candidate_count']
            active += metadata['active_lines']
            polys_total += metadata['polynomials']
        if return_details:
            return (result, dict(mode=mode, projection=g, frequency_cells=result.size, candidates=candidates, active_lines=active, polynomials=polys_total))
        return result
    finally:
        ctx.prec = previous

def tighten_envelopes(weight, Llo, Lhi, box, envelopes, *, mode='dual_tyurin'):
    if envelopes is None:
        return None
    B, F, low = envelopes
    high = coupled_cf_cells(weight.ht[0][:-1], weight.ht[1][1:], Llo, Lhi, box, mode=mode)
    small = coupled_cf_cells(weight.tlo, weight.thi, Llo, Lhi, box, mode=mode)
    if high is None or small is None:
        return envelopes
    return (B, np.minimum(F, high), np.minimum(low, small))

def full_b_cf_cells(tlo, thi, *, L, d, tau, b=None, mode='dual_tyurin', return_details=False):
    Llo, Lhi = map(float, L)
    dlo, dhi = map(float, d)
    taulo, tauhi = map(float, tau)
    blo, bhi = (0.0, Lhi) if b is None else map(float, b)
    return coupled_cf_cells(tlo, thi, Llo, Lhi, (dlo, dhi, blo, bhi, taulo, tauhi), mode=mode, return_details=return_details)
