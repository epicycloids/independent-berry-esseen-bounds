from pathlib import Path
import sys
import math
import time
import json
_THIS = Path(__file__).resolve().parent
_R26 = _THIS.parent.parent / '26'
_source_paths = [_THIS] + ([_R26] if (_R26 / 'paired_bounds.py').is_file() else [])
for _path in map(str, _source_paths):
    if _path not in sys.path:
        sys.path.insert(0, _path)
import numpy as np
from flint import arb, ctx
from interval_bounds import up, au, ia, sqpos
from paired_bounds import rational_poly, alpha_upper, cosine_cells, CosineWeights

def odd_cells(thi, Bhi, dhi, taulo, tauhi):
    tau = min(Bhi, max(0, taulo))
    x = up(thi * au(arb(dhi).sqrt()))
    z = sqpos(ia(np.minimum(x, 2)))
    p = rational_poly(z, [arb(1), arb(1) / 6, -arb(1) / 40, arb(1) / 1008])[1]
    q = up(1 + up(up(x * x) / 3))
    alpha = alpha_upper(x)
    beta2 = up(up(8 * up(p * p)) / 3)
    concavity = up(up(au(arb(8).sqrt()) * q) * alpha)
    D = np.maximum(beta2, concavity)
    excess = max(0, au(arb(Bhi) - arb(tau)))
    first = up(alpha * tau)
    linear = up(up(D * tau) * excess)
    quadratic = up(up(2 * up(q * q)) * up(excess * excess))
    squared = up(up(up(first * first) + linear) + quadratic)
    candidate = np.minimum(Bhi, up(np.sqrt(np.maximum(0, squared))))
    candidate = np.where(x <= 2, candidate, Bhi)
    return np.minimum(cosine_cells(thi, Bhi, dhi, taulo, tauhi), candidate)

class OddWeights(CosineWeights):

    @staticmethod
    def prefactor(thi, Bhi, dhi, taulo, tauhi):
        return odd_cells(thi, Bhi, dhi, taulo, tauhi)

def exact_prefactor(t, B, d, tau):
    t, B, d, tau = map(arb, [t, B, d, tau])
    x = t * d.sqrt()
    if x > 2:
        return B
    q = 1 + x * x / 3
    alpha = 2 * x / 3 - x ** 3 / 15 + x ** 5 / 420
    p = 1 + x * x / 6 - x ** 4 / 40 + x ** 6 / 1008
    D = max(8 * p * p / 3, arb(8).sqrt() * q * alpha)
    v = (alpha ** 2 * tau ** 2 + D * tau * (B - tau) + 2 * q * q * (B - tau) ** 2).sqrt()
    return min(B, v)

def checks():
    oldprec = ctx.prec
    ctx.prec = 160
    count = 0
    strict = 0
    try:
        ts = np.array([0.0, 0.01, 0.1, 0.5, 1.0, 1.5, 1.99, 2.0, 2.01, 3.0, 10.0])
        for B, d, taulo, tauhi in [(0.5, 0.25, 0.35, 0.45), (0.5, 1.0, 0.45, 0.5), (0.5, 0.04, 0.0, 0.5), (0.5, 0.25, 0.5, 0.5), (0.5, 0.0, 0.35, 0.5)]:
            out = odd_cells(ts, B, d, taulo, tauhi)
            old = cosine_cells(ts, B, d, taulo, tauhi)
            assert np.all(out <= old)
            strict += int(np.count_nonzero(out < old))
            for t, y in zip(ts, out):
                for frac in [0.0, 0.5, 1.0]:
                    for tau in [taulo, (taulo + tauhi) / 2, tauhi]:
                        ref = exact_prefactor(t * frac, B, d, tau)
                        assert arb(float(y)) >= ref, (t, B, d, tau, y, ref)
                        count += 1
        return {'outward_point_comparisons': count, 'strict_baseline_improvements': strict, 'all_passed': True}
    finally:
        ctx.prec = oldprec
