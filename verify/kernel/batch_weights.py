import time
import numpy as np
from interval_bounds import up, U, positive_dot
from paired_bounds import CosineWeights
from stable_bounds import finite_law_checks
CHOICES = [(1.0, 0.0), (1.02, -0.02), (1.04, -0.02)]

class BatchWeights(CosineWeights):

    def __init__(self, Lhi, N=512, tau_cap=None):
        self.columns = [CosineWeights(Lhi, N, tau_cap, scale, shift) for scale, shift in CHOICES]
        self.N = N
        self.T = self.columns[0].T
        self.s = self.columns[0].s
        self.log_correction = True
        for name in ['t', 'ht']:
            setattr(self, name, tuple((np.column_stack([getattr(w, name)[j] for w in self.columns]) for j in [0, 1])))
        for name in ['tlo', 'thi', 'gauss', 'inner', 'low', 'outer', 'local', 'high']:
            setattr(self, name, np.column_stack([getattr(w, name) for w in self.columns]))
        self.G = np.array([w.G for w in self.columns])

    def first_terms(self, B, Llo, lowcf=None):
        cumul = up(np.cumsum(up(B * self.inner), axis=0) * (1 + 4 * self.N * U))
        previous = np.vstack([np.zeros((1, len(CHOICES))), cumul[:-1, :]])
        contributions = up(up(self.outer * previous) + up(self.local * B))
        if lowcf is not None:
            cap = up(up(up(self.low * lowcf) + self.outer) / Llo)
            contributions[1:, :] = np.minimum(contributions[1:, :], cap[1:, :])
        return np.array([positive_dot(np.ones(self.N), contributions[:, j]) for j in range(len(CHOICES))])

    def integrals(self, B, F, Llo, lowcf=None):
        i1 = self.first_terms(B, Llo, lowcf)
        i2 = np.array([positive_dot(self.high[:, j], F[:, j]) for j in range(len(CHOICES))])
        return up(i1 + up(up(i2 + self.G) / Llo))

    def integrate(self, B, F, Llo, lowcf=None):
        return float(np.min(self.integrals(B, F, Llo, lowcf)))

def install_factory():
    import stable_cover
    if getattr(stable_cover, '_batch_installed', False):
        return
    original = stable_cover.weight_factory

    def factory(hi, N, kind, target=None):
        if kind != 'batch':
            return original(hi, N, kind, target)
        cache = {}

        def choose(tauhi):
            bucket = min(20, max(1, round(20 * tauhi / hi)))
            if bucket not in cache:
                cache[bucket] = BatchWeights(hi, N, hi * bucket / 20)
            return cache[bucket]
        return choose
    stable_cover.weight_factory = factory
    stable_cover._batch_installed = True
