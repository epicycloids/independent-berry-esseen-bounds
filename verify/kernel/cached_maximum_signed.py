from collections import OrderedDict
import time
if __package__:
    from .full_b_paths import HERE
else:
    from full_b_paths import HERE
import numpy as np
from flint import ctx
from interval_signed import Q, _upper, up, down, normalize, directed_prefix, positive_dot, cell_contributions
from fast_maximum_variance_signed import FastMaximumVarianceRectangles
import stable_bounds

class CachedMaximumVarianceRectangles(FastMaximumVarianceRectangles):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._initialize_score_cache()

    def _initialize_score_cache(self):
        self.width_low = np.column_stack([[float(down(float(2 * (r - l)))) for l, r in zip(nodes, nodes[1:])] for nodes in self.nodes])
        self.width_low.flags.writeable = False
        self.score_cache = OrderedDict()
        self.score_cache_hits = self.score_cache_misses = 0

    def _proposal(self, Llo, Lhi, dlo, dhi, taulo, tauhi):
        key = (Llo, Lhi, dlo, dhi, taulo, tauhi, ctx.prec, max(ctx.prec, 128))
        if key in self.score_cache:
            self.score_cache.move_to_end(key)
            self.score_cache_hits += 1
            return self.score_cache[key]
        self.score_cache_misses += 1
        w = self.weights
        E, offset_lo, offset_hi, eligible = self._reference(Llo, Lhi, dlo, dhi, taulo, tauhi)
        eligible = eligible.copy()
        eligible[0, :] = False
        E = np.where(eligible, E, 0.0)
        offset_lo = np.where(eligible, offset_lo, 0.0)
        offset_hi = np.where(eligible, offset_hi, 0.0)
        product = np.maximum.reduce([up(self.qlo * offset_lo), up(self.qlo * offset_hi), up(self.qhi * offset_lo), up(self.qhi * offset_hi)])
        center = up(product * np.where(product >= 0, self.width_twice, self.width_low))
        numerator = up(up(E * w.low) + center)
        proposed = normalize(numerator, Llo, Lhi)
        proposed.flags.writeable = False
        eligible.flags.writeable = False
        self.score_cache[key] = (proposed, eligible)
        if len(self.score_cache) > 16:
            self.score_cache.popitem(last=False)
        return (proposed, eligible)

    def score(self, Llo, Lhi, box, envelopes, *, detail=False):
        w = self.weights
        clipped = stable_bounds.feasible_clip(box, Lhi)
        if clipped is None:
            return None
        dlo, dhi, _, _, taulo, tauhi = clipped
        B, F, lowcf = envelopes
        old_cells = cell_contributions(w, B, Llo, lowcf)
        old = float(np.min(w.cutoff_integrals(B, F, Llo, lowcf)))
        proposed, eligible = self._proposal(Llo, Lhi, dlo, dhi, taulo, tauhi)
        cells = np.where(eligible, np.minimum(old_cells, proposed), old_cells)
        prefix = directed_prefix(cells)
        cf = up(w.low * lowcf)
        suffix = np.zeros((self.N + 1, self.columns))
        for i in range(self.N - 1, -1, -1):
            suffix[i, :] = up(suffix[i + 1, :] + cf[i, :])
        high = np.array([positive_dot(w.high[:, j], F[:, j]) for j in range(self.columns)])
        scores = np.empty_like(self.gamma)
        for r, k in enumerate(w.cutoff_indices):
            common = up(up(suffix[k, :] + high) + self.gamma[:, r, :])
            scores[:, r, :] = up(prefix[:, k, :] + normalize(common, Llo, Lhi))
        per_interval = np.min(scores, axis=(1, 2))
        inside = float(np.max(per_interval))
        X = self.theta[-1]
        outside = _upper(1 / (Q(Llo) * (1 + X * X)))
        new = max(inside, outside)
        answer = min(old, new)
        if not detail:
            return answer
        peak = int(np.argmax(per_interval))
        return dict(upper=answer, old_upper=old, signed_upper=new, inside_upper=inside, outside_upper=outside, peak_threshold_interval=list(map(str, self.theta[peak:peak + 2])), threshold_intervals=len(per_interval), per_interval=per_interval.tolist(), eligible_cells=int(np.count_nonzero(eligible)), improved_cells=int(np.count_nonzero(cells < old_cells)), negative_cells=int(np.count_nonzero(cells < 0)), build_seconds=self.build_seconds, status='outward signed reference bound')
