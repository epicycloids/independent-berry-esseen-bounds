import time
import numpy as np
from flint import arb
from interval_bounds import up, au
from stable_bounds import StableWeights, finite_law_checks

def directional_cells(thi, Bhi, dhi, taulo, tauhi):
    tau = min(Bhi, max(0, taulo))
    x = up(thi * au(arb(dhi).sqrt()))
    q = up(1 + up(up(x * x) / 3))
    alpha = up(up(2 * x) / 3)
    excess = max(0, float(up(Bhi - tau)))
    other = au(arb(Bhi) + 5 * arb(tau) / 3)
    first = up(up(alpha * tau) + up(q * excess))
    squared = up(up(first * first) + up(up(up(q * q) * excess) * other))
    return np.minimum(Bhi, up(np.sqrt(np.maximum(0, squared))))

class DirectionalWeights(StableWeights):

    @staticmethod
    def prefactor(thi, Bhi, dhi, taulo, tauhi):
        return directional_cells(thi, Bhi, dhi, taulo, tauhi)

def directional_checks():
    rng = np.random.default_rng(260927)
    count = 0
    for _ in range(100):
        B = float(rng.uniform(0.001, 2))
        tl = float(rng.uniform(0, B))
        th = float(rng.uniform(tl, B))
        d = float(rng.uniform(0, 1))
        ts = np.array([0.0, 0.001, 0.1, 0.5, 1.0, 2.0, 4.0, 10.0])
        out = directional_cells(ts, B, d, tl, th)
        for tau in [tl, (tl + th) / 2, th]:
            bb, tt, dd = (arb(B), arb(tau), arb(d))
            for t, y in zip(ts, out):
                x = arb(float(t)) * dd.sqrt()
                q = 1 + x * x / 3
                ref = ((2 * x * tt / 3 + q * (bb - tt)) ** 2 + q * q * (bb - tt) * (bb + 5 * tt / 3)).sqrt()
                assert y == B or arb(float(y)) >= ref, (B, d, tl, th, tau, t, y, ref)
                count += 1
    return {'directional_enclosures': count, 'all_passed': True}
