from collections import OrderedDict
from fractions import Fraction as Q
import math
from pathlib import Path
import sys
import time
import numpy as np
from flint import arb, ctx
HERE = Path(__file__).resolve().parent
if (HERE.parent / 'scalar_table_cover_engine.py').is_file():
    R27 = HERE.parent
    for path in (R27.parent / '26', R27, R27 / 'smoothing_geometry', R27 / 'smoothing_geometry/adjoint', R27 / 'scalar_table', R27 / 'dual_logarithm', R27 / 'farfield_phase/new_target'):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))
from interval_bounds import up, down, au, al, positive_dot
from cutoff_transfer import cell_contributions
from scalar_table_cover_engine import ScalarCutoffWeights, selected, CUTOFFS
from reference_bounds import reference_bounds
from signed_gaussian import verify_partition
import stable_bounds

def _arb(x):
    x = Q(x)
    return arb(x.numerator) / x.denominator

def _upper(x):
    x = Q(x)
    value = float(x)
    return math.nextafter(value, math.inf) if Q(value) < x else value

def normalize(upper_numerator, Llo, Lhi):
    if not (math.isfinite(Llo) and math.isfinite(Lhi) and (0 < Llo <= Lhi)):
        raise ValueError('Require a positive ordered L interval')
    value = np.asarray(upper_numerator)
    return up(value / np.where(value >= 0, Llo, Lhi))

def directed_prefix(cells):
    cells = np.asarray(cells, dtype=float)
    result = np.empty((cells.shape[0], cells.shape[1] + 1, cells.shape[2]))
    result[:, 0, :] = 0.0
    for j in range(cells.shape[1]):
        result[:, j + 1, :] = up(result[:, j, :] + cells[:, j, :])
    return result

def gaussian_subinterval(scalar, left, right):
    left, right = (Q(left), Q(right))
    if not -Q(scalar['X']) <= left < right <= Q(scalar['X']):
        raise ValueError('Requested threshold interval lies outside the saved partition')
    nodes = {Q(row['x']): Q(row['upper']) for row in scalar['nodes']}
    curvature = Q(scalar['curvature_upper'])
    if curvature < 0:
        raise ValueError('A curvature upper bound must be nonnegative')
    values, cursor = ([], left)
    for cell in scalar['intervals']:
        a, b = (Q(cell['left']), Q(cell['right']))
        lo, hi = (max(left, a), min(right, b))
        if hi <= lo:
            continue
        if lo != cursor:
            raise ValueError('The saved threshold partition has a gap')
        ya, yb = (nodes[a], nodes[b])

        def line(x):
            return ((b - x) * ya + (x - a) * yb) / (b - a)
        vertex = min(hi, max(lo, (a + b) / 2))
        remainder = curvature * (vertex - a) * (b - vertex) / 2
        values.append(max(line(lo), line(hi)) + remainder)
        cursor = hi
    if cursor != right or not values:
        raise ValueError('The saved threshold partition did not cover the request')
    return max(values)

def real_kernel_range(T, left, right, theta_left, theta_right):
    if not (0 < left < right < T and theta_left < theta_right):
        raise ValueError('The kernel rectangle must avoid zero and the cutoff pole')
    tt = _arb(left).union(_arb(right))
    theta = _arb(theta_left).union(_arb(theta_right))
    T = _arb(T)
    u, phase = (tt / T, theta * tt)
    z = arb.pi() * u
    h = 1 / z - z.cos() / z.sin()
    twice = (1 - u) / T * (phase.cos() - h * phase.sin()) + theta / arb.pi() * phase.sinc()
    if not twice.is_finite():
        raise ArithmeticError('Nonfinite signed-kernel enclosure')
    return (al(twice / 2), au(twice / 2))

class SignedRectangles:

    def __init__(self, weights, theta_edges=None):
        began = time.monotonic()
        self.weights = weights
        self.theta = tuple((Q(j, 32) for j in range(-64, 65))) if theta_edges is None else tuple(map(Q, theta_edges))
        if len(self.theta) < 2 or self.theta[0] != -self.theta[-1] or self.theta[-1] <= 0 or any((a >= b for a, b in zip(self.theta, self.theta[1:]))):
            raise ValueError('Require increasing thresholds covering a symmetric interval')
        self.N, self.columns = (weights.N, len(weights.columns))
        H = len(self.theta) - 1
        self.qlo = np.zeros((H, self.N, self.columns))
        self.qhi = np.zeros_like(self.qlo)
        self.width_twice = np.empty((self.N, self.columns))
        self.nodes = []
        self.gamma = np.empty((H, len(weights.cutoff_indices), self.columns))
        self.reference_cache = OrderedDict()
        previous_precision = ctx.prec
        try:
            ctx.prec = max(previous_precision, 128)
            for j, column in enumerate(weights.columns):
                T, s = (Q(column.T), Q(column.s))
                nodes = tuple((T * s * i / self.N for i in range(self.N + 1)))
                self.nodes.append(nodes)
                for i, (left, right) in enumerate(zip(nodes, nodes[1:])):
                    self.width_twice[i, j] = _upper(2 * (right - left))
                    if i == 0:
                        continue
                    for h, (a, b) in enumerate(zip(self.theta, self.theta[1:])):
                        self.qlo[h, i, j], self.qhi[h, i, j] = real_kernel_range(T, left, right, a, b)
                for r, _ in enumerate(weights.cutoff_indices):
                    record = weights.cutoff_records[r][j]
                    scalar = record['scalar_certificate']
                    verify_partition(scalar)
                    correction = Q(record['endpoint_correction'])
                    if correction < 0:
                        raise ValueError('A cutoff adjustment must be nonnegative')
                    for h, (a, b) in enumerate(zip(self.theta, self.theta[1:])):
                        gamma = gaussian_subinterval(scalar, a, b) + correction
                        self.gamma[h, r, j] = min(float(weights.cutoff_G[r, j]), _upper(gamma))
        finally:
            ctx.prec = previous_precision
        self.build_seconds = time.monotonic() - began

    def _reference(self, Llo, Lhi, dhi, taulo, tauhi):
        key = (Llo, Lhi, dhi, taulo, tauhi)
        if key in self.reference_cache:
            self.reference_cache.move_to_end(key)
            return self.reference_cache[key]
        w = self.weights
        deleted = selected(w.tlo, w.thi, 1.0, 1.0, Lhi, dhi, taulo, tauhi, deleted=True)
        E, lower, upper, eligible = ([], [], [], [])
        for j, nodes in enumerate(self.nodes):
            result = reference_bounds(nodes, deleted[:, j], L=(Llo, Lhi), d=(0.0, dhi), tau=(taulo, tauhi))
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
        _, dhi, _, _, taulo, tauhi = clipped
        B, F, lowcf = envelopes
        old_cells = cell_contributions(w, B, Llo, lowcf)
        old = float(np.min(w.cutoff_integrals(B, F, Llo, lowcf)))
        E, offset_lo, offset_hi, eligible = self._reference(Llo, Lhi, dhi, taulo, tauhi)
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
