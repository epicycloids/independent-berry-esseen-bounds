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
        return dict(upper=answer, old_upper=old, signed_upper=new, inside_upper=inside, outside_upper=outside, peak_threshold_interval=list(map(str, self.theta[peak:peak + 2])), threshold_intervals=len(per_interval), per_interval=per_interval.tolist(), eligible_cells=int(np.count_nonzero(eligible)), improved_cells=int(np.count_nonzero(cells < old_cells)), negative_cells=int(np.count_nonzero(cells < 0)), build_seconds=self.build_seconds, status='author outward signed reference comparison; independent review required')

def profile_boxes(Llo, Lhi, boxes, mesh=512, partitions=(1, 2, 4, 8), seconds_limit=50):
    if not 0 < Llo <= Lhi or not 1 <= len(boxes) <= 3:
        raise ValueError('Require a positive L interval and one to three boxes')
    if mesh not in (128, 256, 512, 1024) or not 0 < seconds_limit <= 120:
        raise ValueError('Require a bounded mesh and a limit of at most 120 seconds')
    if not partitions or any((n not in (1, 2, 4, 8) for n in partitions)):
        raise ValueError('Use one-, two-, four-, or eight-slice partitions')
    from full_b_consumer import full_b_box, closed_b_slices
    from fast_maximum_cover_engine import FastMaximumVarianceWeights, selected, CUTOFFS
    began = time.monotonic()
    previous = stable_bounds.log_range_cells
    totals = dict(envelope_cpu=0.0, original_score_cpu=0.0, cached_score_cpu=0.0)
    reference = {name: dict(cold_calls=0, hit_calls=0, cold_cpu=0.0, hit_cpu=0.0) for name in ('original', 'cached')}
    rows = []
    stopped = False

    def timed(call):
        started = time.process_time()
        result = call()
        return (result, time.process_time() - started)

    def observe(consumer, name):
        method = consumer._reference

        def wrapped(*args):
            key = (*args, ctx.prec, max(ctx.prec, 128))
            label = 'hit' if key in consumer.reference_cache else 'cold'
            result, seconds = timed(lambda: method(*args))
            reference[name][label + '_calls'] += 1
            reference[name][label + '_cpu'] += seconds
            return result
        consumer._reference = wrapped
    try:
        stable_bounds.log_range_cells = selected
        weights, weight_cpu = timed(lambda: FastMaximumVarianceWeights(Lhi, mesh, 0.85 * Lhi, cutoff_fractions=CUTOFFS))
        original, original_build = timed(lambda: FastMaximumVarianceRectangles(weights))
        cached, cached_build = timed(lambda: CachedMaximumVarianceRectangles(weights))
        assert original.nodes == cached.nodes and original.theta == cached.theta
        for name in ('qlo', 'qhi', 'width_twice', 'gamma'):
            assert np.array_equal(getattr(original, name), getattr(cached, name)), name
        observe(original, 'original')
        observe(cached, 'cached')
        for index, box in enumerate(boxes):
            full = full_b_box(Lhi, box)
            if full is None:
                rows.append(dict(box_index=index, excluded=True))
                continue
            for count in partitions:
                values = []
                complete = True
                for piece in closed_b_slices(full, count):
                    if time.monotonic() - began >= seconds_limit:
                        stopped = True
                        complete = False
                        break
                    env, cpu = timed(lambda: weights.envelopes(Llo, Lhi, piece))
                    totals['envelope_cpu'] += cpu
                    if env is None:
                        continue
                    x, cpu = timed(lambda: original.score(Llo, Lhi, full, env, detail=True))
                    totals['original_score_cpu'] += cpu
                    y, cpu = timed(lambda: cached.score(Llo, Lhi, full, env, detail=True))
                    totals['cached_score_cpu'] += cpu
                    x.pop('build_seconds')
                    y.pop('build_seconds')
                    assert x == y, (index, count, piece, x, y)
                    values.append(x['upper'])
                rows.append(dict(box_index=index, full_box=full, parts=count, complete=complete, evaluated_slices=len(values), upper=max(values) if complete and values else None))
                if stopped:
                    break
            if stopped:
                break
        return dict(status='whole-box exact score comparison; no complete continuum cover', L=[Llo, Lhi], mesh=mesh, threshold_grid='default 1/32 on [-2,2]', rows=rows, exact_outputs_equal=True, stopped_at_limit=stopped, construction_cpu=dict(weights=weight_cpu, original=original_build, cached=cached_build), component_cpu=totals, reference=reference, score_cache=dict(hits=cached.score_cache_hits, misses=cached.score_cache_misses), elapsed_seconds=time.monotonic() - began)
    finally:
        stable_bounds.log_range_cells = previous

def profile_point(point, mesh=512, halfwidth=0.0005, **kwargs):
    L, d, tau = map(float, point)
    if not 0 <= halfwidth < L:
        raise ValueError('Require a nonnegative halfwidth smaller than L')
    return profile_boxes(L - halfwidth, L + halfwidth, [(d - halfwidth, d + halfwidth, 0.0, L + halfwidth, tau - halfwidth, tau + halfwidth)], mesh=mesh, **kwargs)
