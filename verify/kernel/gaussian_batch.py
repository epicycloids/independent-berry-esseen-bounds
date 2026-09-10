import math
import time
from functools import lru_cache
import numpy as np
from flint import arb, acb
from scipy.integrate import quad
from scipy.special import exp1
from interval_bounds import au
from batch_weights import BatchWeights

@lru_cache(maxsize=2048)
def gaussian_G(T_value, s_value, M=8):
    T, s = (arb(T_value), arb(s_value))
    p = arb.pi()
    assert T > 0 and 0 < s < 1
    c = au((p * p / 6 - sum((arb(1) / (n * n) for n in range(1, M + 1)))) / (1 - s * s / (M + 1) ** 2))

    def integrand(u, analytic):
        h = 2 * u / p * (sum((1 / (n * n - u * u) for n in range(1, M + 1))) + c)
        return (1 - u) * (1 + h * h).sqrt(analytic=analytic) * (-T * T * u * u / 2).exp()
    value = acb.integral(integrand, acb(0), acb(s), rel_tol=1e-12, abs_tol=1e-14, eval_limit=15000, depth_limit=30, use_heap=True)
    assert value.is_finite() and value.imag.contains(0), value
    total = value.real + (T * T * s * s / 2).expint(1) / (2 * p)
    assert total > 0
    return au(total)

class GaussianBatchWeights(BatchWeights):

    def __init__(self, Lhi, N=512, tau_cap=None):
        super().__init__(Lhi, N, tau_cap)
        self.G = np.minimum(self.G, [gaussian_G(w.T, w.s) for w in self.columns])

def install_factory():
    import stable_cover
    if getattr(stable_cover, '_gaussian_installed', False):
        return
    original = stable_cover.weight_factory

    def factory(hi, N, kind, target=None):
        if kind != 'gaussian':
            return original(hi, N, kind, target)
        cache = {}

        def choose(tauhi):
            bucket = min(20, max(1, round(20 * tauhi / hi)))
            if bucket not in cache:
                cache[bucket] = GaussianBatchWeights(hi, N, hi * bucket / 20)
            return cache[bucket]
        return choose
    stable_cover.weight_factory = factory
    stable_cover._gaussian_installed = True
