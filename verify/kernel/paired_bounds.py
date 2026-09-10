import math
import time
import numpy as np
from flint import arb
from interval_bounds import up, down, au, al, ia, add, mul, sqpos, cubepos
from even_bounds import EvenWeights, even_checks
from stable_bounds import finite_law_checks

def rational_poly(z, coefficients):
    p = (al(coefficients[-1]), au(coefficients[-1]))
    for c in coefficients[-2::-1]:
        p = add(mul(p, z), (al(c), au(c)))
    return p

def cos_abs_upper(xlo, xhi):
    safe_lo = np.clip(xlo, 0, 3)
    safe_hi = np.clip(xhi, 0, 3)
    co = [arb((-1) ** j) / math.factorial(2 * j) for j in range(10)]
    upper = rational_poly(sqpos(ia(safe_lo)), co[:-1])[1]
    lower = rational_poly(sqpos(ia(safe_hi)), co)[0]
    return np.minimum(1, np.where(xhi <= 3, np.maximum(upper, -lower), 1))

def paired_cells(t, dlo, dhi, bhi):
    x = mul(t, (al(arb(dlo).sqrt()), au(arb(dhi).sqrt())))
    C = cos_abs_upper(np.maximum(0, x[0][:-1]), x[1][1:])
    tau = min(bhi, max(0, al(arb(dlo) ** (arb(3) / 2))))
    excess = max(0, au(arb(bhi) - arb(tau)))
    other = au(2 * arb(bhi) + 2 * arb(tau) / 3)
    t3 = cubepos(ia(t[1][1:]))[1]
    term1 = up(C * C)
    term2 = up(up(up(C * excess) * t3) / 3)
    term3 = up(up(up(up(t3 * t3) * excess) * other) / 36)
    return np.minimum(1, up(np.sqrt(up(up(term1 + term2) + term3))))

def alpha_upper(x):
    z = sqpos(ia(np.minimum(x, 2)))
    co = [arb(2) / 3, -arb(1) / 15, arb(1) / 420]
    p = rational_poly(z, co)[1]
    old = up(up(2 * x) / 3)
    return np.where(x <= 2, np.minimum(old, up(p * x)), old)

def cosine_cells(thi, Bhi, dhi, taulo, tauhi):
    tau = min(Bhi, max(0, taulo))
    x = up(thi * au(arb(dhi).sqrt()))
    q = up(1 + up(up(x * x) / 3))
    alpha = alpha_upper(x)
    excess = max(0, au(arb(Bhi) - arb(tau)))
    other = au(2 * arb(Bhi) + 2 * arb(tau) / 3)
    first = up(alpha * tau)
    squared = up(up(first * first) + up(up(up(q * q) * excess) * other))
    return np.minimum(Bhi, up(np.sqrt(np.maximum(0, squared))))

class PairedWeights(EvenWeights):

    def single_cells(self, t, dlo, dhi, bhi):
        return np.minimum(super().single_cells(t, dlo, dhi, bhi), paired_cells(t, dlo, dhi, bhi))

class CosineWeights(PairedWeights):
    direct_sum = True

    @staticmethod
    def prefactor(thi, Bhi, dhi, taulo, tauhi):
        return cosine_cells(thi, Bhi, dhi, taulo, tauhi)

def checks():
    rng = np.random.default_rng(260929)
    cs = 0
    ps = 0
    zs = 0
    for _ in range(100):
        v = float(rng.uniform(0.0001, 1))
        b = float(v ** 1.5 * rng.uniform(1, 3))
        ts = np.array([0.0, 0.001, 0.1, 0.5, 1.0, 2.0, 4.0, 10.0])
        x = ts * math.sqrt(v)
        C = cos_abs_upper(x, x)
        for xx, y in zip(x, C):
            assert arb(float(y)) >= arb(float(xx)).cos().abs_upper()
            cs += 1
        out = paired_cells((ts, ts), v, v, b)
        for a, z, y in zip(ts[:-1], ts[1:], out):
            for tt in [a, (a + z) / 2, z]:
                t = arb(float(tt))
                tau = arb(v) ** (arb(3) / 2)
                bb = arb(b)
                c = (t * arb(v).sqrt()).cos().abs_upper()
                ref = (c * c + c * (bb - tau) * t ** 3 / 3 + t ** 6 * (bb - tau) * (2 * bb + 2 * tau / 3) / 36).sqrt()
                assert y == 1 or arb(float(y)) >= ref, (v, b, tt, y, ref)
                ps += 1
        tl = float(rng.uniform(0, b))
        th = float(rng.uniform(tl, b))
        out = cosine_cells(ts, b, v, tl, th)
        for tau in [tl, (tl + th) / 2, th]:
            for tt, y in zip(ts, out):
                xx = arb(float(tt)) * arb(v).sqrt()
                q = 1 + xx * xx / 3
                alpha = 2 * xx / 3 - xx ** 3 / 15 + xx ** 5 / 420 if xx <= 2 else 2 * xx / 3
                tau0 = arb(tau)
                bb = arb(b)
                ref = (alpha ** 2 * tau0 ** 2 + q * q * (bb - tau0) * (2 * bb + 2 * tau0 / 3)).sqrt()
                assert y == b or arb(float(y)) >= ref, (v, b, tt, tau, y, ref)
                zs += 1
    return {'cosine_endpoint_enclosures': cs, 'paired_cf_cell_checks': ps, 'cosine_zero_bias_enclosures': zs, 'all_passed': True}

def probe():
    import log_probe
    preceding_single = log_probe.cf_single

    def paired_single(t, v, b):
        tau = v ** 1.5
        excess = max(0, b - tau)
        c = np.abs(np.cos(t * math.sqrt(v)))
        bound = np.sqrt(c * c + c * excess * t ** 3 / 3 + t ** 6 * excess * (2 * b + 2 * tau / 3) / 36)
        return np.minimum(preceding_single(t, v, b), bound)

    def cosine(t, B, d, tau, directional=False):
        excess = max(0, B - tau)
        x = t * math.sqrt(d)
        q = 1 + x * x / 3
        alpha = np.where(x <= 2, 2 * x / 3 - x ** 3 / 15 + x ** 5 / 420, 2 * x / 3)
        return np.minimum(B, np.sqrt((alpha * tau) ** 2 + q * q * excess * (2 * B + 2 * tau / 3)))
    log_probe.cf_single = paired_single
    log_probe.stability = cosine
    result = log_probe.run_stable(True)
    result['refinements'] = ['paired CF', 'cosine zero bias']
    return result
