from functools import lru_cache
import json
from pathlib import Path
import time
import numpy as np
from interval_signed import SignedRectangles, Q, _upper, up, down, normalize, directed_prefix, positive_dot, cell_contributions
from packed_energy_reference import packed_energy_reference_bounds
try:
    from maximum_variance_bounds import refine_reference
except ModuleNotFoundError as error:
    if error.name != 'maximum_variance_bounds':
        raise
    import importlib.util
    path = Path(__file__).resolve().parent.parent / 'farfield_phase/new_target/maximum_variance_reference/bounds.py'
    spec = importlib.util.spec_from_file_location('local_maximum_variance_bounds', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    refine_reference = module.refine_reference
from complete_energy_signed import load_complete_energy
from packed_energy_signed import PackedEnergyRectangles
from disk_reference_cover_engine import DiskReferenceWeights, selected, CUTOFFS
import stable_bounds

class MaximumVarianceRectangles(SignedRectangles):

    def _reference(self, Llo, Lhi, dlo, dhi, taulo, tauhi):
        key = (Llo, Lhi, dlo, dhi, taulo, tauhi)
        if key in self.reference_cache:
            self.reference_cache.move_to_end(key)
            return self.reference_cache[key]
        w = self.weights
        deleted = selected(w.tlo, w.thi, 1.0, 1.0, Lhi, dhi, taulo, tauhi, deleted=True)
        E, lower, upper, eligible = ([], [], [], [])
        for j, nodes in enumerate(self.nodes):
            result = packed_energy_reference_bounds(nodes, deleted[:, j], energy_table=load_complete_energy(), L=(Llo, Lhi), d=(dlo, dhi), tau=(taulo, tauhi))
            result = refine_reference(result, nodes, d=(dlo, dhi), tau=(taulo, tauhi))
            E.append(result['E_upper'])
            lower.append(result['lower_offset'])
            upper.append(result['upper_offset'])
            eligible.append(result['eligible'])
        result = tuple((np.asarray(values).T for values in (E, lower, upper, eligible)))
        self.reference_cache[key] = result
        if len(self.reference_cache) > 16:
            self.reference_cache.popitem(last=False)
        return result

    def score(self, Llo, Lhi, box, envelopes, *, detail=False):
        w = self.weights
        clipped = stable_bounds.feasible_clip(box, Lhi)
        if clipped is None:
            return None
        dlo, dhi, _, _, taulo, tauhi = clipped
        B, F, lowcf = envelopes
        old_cells = cell_contributions(w, B, Llo, lowcf)
        old = float(np.min(w.cutoff_integrals(B, F, Llo, lowcf)))
        E, offset_lo, offset_hi, eligible = self._reference(Llo, Lhi, dlo, dhi, taulo, tauhi)
        eligible = eligible.copy()
        eligible[0, :] = False
        E = np.where(eligible, E, 0.0)
        offset_lo = np.where(eligible, offset_lo, 0.0)
        offset_hi = np.where(eligible, offset_hi, 0.0)
        product = np.maximum.reduce([up(self.qlo * offset_lo), up(self.qlo * offset_hi), up(self.qhi * offset_lo), up(self.qhi * offset_hi)])
        width_low = np.column_stack([[float(down(float(2 * (r - l)))) for l, r in zip(nodes, nodes[1:])] for nodes in self.nodes])
        center = up(product * np.where(product >= 0, self.width_twice, width_low))
        numerator = up(up(E * w.low) + center)
        proposed = normalize(numerator, Llo, Lhi)
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
        return dict(upper=answer, old_upper=old, signed_upper=new, inside_upper=inside, outside_upper=outside, peak_threshold_interval=list(map(str, self.theta[peak:peak + 2])), threshold_intervals=len(per_interval), per_interval=per_interval.tolist(), eligible_cells=int(np.count_nonzero(eligible)), improved_cells=int(np.count_nonzero(cells < old_cells)), negative_cells=int(np.count_nonzero(cells < 0)), build_seconds=self.build_seconds, status='author outward signed reference comparison; independent review required')
