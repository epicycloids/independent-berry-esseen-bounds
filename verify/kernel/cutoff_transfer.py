from functools import lru_cache
from pathlib import Path
import math
import sys
import numpy as np
from flint import arb
HERE = Path(__file__).resolve().parent
for path in (HERE.parent.parent.parent / '26', HERE.parent.parent, HERE.parent):
    if path.exists():
        sys.path.insert(0, str(path))
from interval_bounds import U, up, au, positive_dot, two_kernel
from paired_bounds import CosineWeights
from gaussian_batch import gaussian_G
from batch_weights import CHOICES
from pilot_signed_gaussian import signed_certificate
DEFAULT_CUTOFFS = (0.75, 0.8125, 0.875, 0.9375, 1.0)

@lru_cache(maxsize=16384)
def grid_cutoff_certificate(T_value, s_value, k, N):
    if not 1 <= k <= N:
        raise ValueError('The first frequency cell must remain below the split.')
    exact = arb(float(s_value)) * int(k) / int(N)
    rounded = float(float(s_value) * int(k) / int(N))
    approximate = arb(rounded)
    interval = exact.union(approximate)
    correction = au(abs(exact - approximate) * two_kernel(interval))
    scalar = signed_certificate(float(T_value), rounded)
    result = au(arb(scalar['upper']) + arb(correction))
    return dict(upper=result, exact_split='s*k/N', rounded_split=rounded, endpoint_correction=correction, k=int(k), N=int(N), scalar_certificate=scalar)

def cell_contributions(weights, B, Llo, lowcf):
    if not (math.isfinite(Llo) and Llo > 0):
        raise ValueError('A positive Llo is required.')
    B, lowcf = (np.asarray(B), np.asarray(lowcf))
    if B.shape != weights.inner.shape or lowcf.shape != weights.low.shape:
        raise ValueError('Envelope and weight shapes must agree.')
    cumul = up(np.cumsum(up(B * weights.inner), axis=0) * (1 + 4 * weights.N * U))
    previous = np.vstack([np.zeros((1, B.shape[1])), cumul[:-1, :]])
    contributions = up(up(weights.outer * previous) + up(weights.local * B))
    cap = up(up(up(weights.low * lowcf) + weights.outer) / Llo)
    contributions[1:, :] = np.minimum(contributions[1:, :], cap[1:, :])
    return contributions

class CutoffBatchWeights(CosineWeights):

    def __init__(self, Lhi, N=512, tau_cap=None, choices=None, cutoff_fractions=DEFAULT_CUTOFFS):
        self.choices = tuple(CHOICES if choices is None else choices)
        if not self.choices:
            raise ValueError('At least one smoothing choice is required.')
        self.columns = [CosineWeights(Lhi, N, tau_cap, scale, shift) for scale, shift in self.choices]
        self.N = int(N)
        self.T, self.s = (self.columns[0].T, self.columns[0].s)
        self.log_correction = True
        for name in ('t', 'ht'):
            setattr(self, name, tuple((np.column_stack([getattr(w, name)[j] for w in self.columns]) for j in (0, 1))))
        for name in ('tlo', 'thi', 'gauss', 'inner', 'low', 'outer', 'local', 'high'):
            setattr(self, name, np.column_stack([getattr(w, name) for w in self.columns]))
        self.G = np.array([min(w.G, gaussian_G(w.T, w.s), signed_certificate(w.T, w.s)['upper']) for w in self.columns])
        fractions = tuple((float(z) for z in cutoff_fractions))
        if not fractions or not all((math.isfinite(z) and 0 < z <= 1 for z in fractions)):
            raise ValueError('Cutoff fractions must lie in (0,1].')
        self.cutoff_indices = tuple(sorted(set([self.N] + [min(self.N, max(1, int(round(self.N * z)))) for z in fractions])))
        self.cutoff_G = np.empty((len(self.cutoff_indices), len(self.columns)))
        self.cutoff_records = []
        for r, k in enumerate(self.cutoff_indices):
            records = []
            for j, w in enumerate(self.columns):
                if k == self.N:
                    self.cutoff_G[r, j] = self.G[j]
                    records.append(dict(k=k, N=self.N, upper=float(self.G[j]), rounded_split=w.s, endpoint_correction=0.0, scalar_certificate=signed_certificate(w.T, w.s)))
                else:
                    certificate = grid_cutoff_certificate(w.T, w.s, k, self.N)
                    self.cutoff_G[r, j] = certificate['upper']
                    records.append(certificate)
            self.cutoff_records.append(records)

    def first_terms(self, B, Llo, lowcf):
        contributions = cell_contributions(self, B, Llo, lowcf)
        return np.array([positive_dot(np.ones(self.N), contributions[:, j]) for j in range(len(self.columns))])

    def cutoff_integrals(self, B, F, Llo, lowcf):
        contributions = cell_contributions(self, B, Llo, lowcf)
        prefix = np.vstack([np.zeros((1, len(self.columns))), up(np.cumsum(contributions, axis=0) * (1 + 4 * self.N * U))])
        cf_cells = up(self.low * lowcf)
        suffix = np.vstack([up(np.cumsum(cf_cells[::-1, :], axis=0)[::-1, :] * (1 + 4 * self.N * U)), np.zeros((1, len(self.columns)))])
        high = np.array([positive_dot(self.high[:, j], F[:, j]) for j in range(len(self.columns))])
        scores = []
        for r, k in enumerate(self.cutoff_indices):
            numerator = up(up(suffix[k, :] + high) + self.cutoff_G[r, :])
            scores.append(up(prefix[k, :] + up(numerator / Llo)))
        return np.vstack(scores)

    def integrals(self, B, F, Llo, lowcf):
        return np.min(self.cutoff_integrals(B, F, Llo, lowcf), axis=0)

    def integrate(self, B, F, Llo, lowcf):
        return float(np.min(self.integrals(B, F, Llo, lowcf)))
