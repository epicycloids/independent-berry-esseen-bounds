from bisect import bisect_right
from collections import OrderedDict
from fractions import Fraction as Q
import time
import numpy as np
from flint import arb, ctx
import interval_signed as signed
from packed_energy_signed import PackedEnergyRectangles, selected, load_complete_energy
from fast_packed_reference import PreparedPrefix, PreparedGeometry, PreparedBox, fast_packed_energy_reference_bounds, rational
import packed_energy_reference as frozen_reference

def exact(value):
    value = rational(value)
    return arb(value.numerator) / value.denominator

class PreparedGaussian:

    def __init__(self, scalar):
        self.X = rational(scalar['X'])
        self.curvature = rational(scalar['curvature_upper'])
        if self.curvature < 0:
            raise ValueError('A curvature upper bound must be nonnegative')
        nodes = {rational(row['x']): rational(row['upper']) for row in scalar['nodes']}
        rows, cursor = ([], -self.X)
        for cell in scalar['intervals']:
            a, b = (rational(cell['left']), rational(cell['right']))
            if not (a == cursor and a < b):
                raise ValueError('The saved threshold partition has a gap or overlap')
            rows.append((a, b, nodes[a], nodes[b], b - a, (a + b) / 2))
            cursor = b
        if cursor != self.X or not rows:
            raise ValueError('The saved threshold partition is incomplete')
        self.rows = tuple(rows)
        self.rights = tuple((row[1] for row in rows))

    def upper(self, left, right):
        left, right = (rational(left), rational(right))
        if not -self.X <= left < right <= self.X:
            raise ValueError('Requested threshold interval lies outside the saved partition')
        values, cursor = ([], left)
        start = bisect_right(self.rights, left)
        for a, b, ya, yb, width, middle in self.rows[start:]:
            if a >= right:
                break
            lo, hi = (max(left, a), min(right, b))
            if hi <= lo:
                continue
            if lo != cursor:
                raise ValueError('The saved threshold partition has a gap')
            low_line = ((b - lo) * ya + (lo - a) * yb) / width
            high_line = ((b - hi) * ya + (hi - a) * yb) / width
            vertex = min(hi, max(lo, middle))
            remainder = self.curvature * (vertex - a) * (b - vertex) / 2
            values.append(max(low_line, high_line) + remainder)
            cursor = hi
        if cursor != right or not values:
            raise ValueError('The saved threshold partition did not cover the request')
        return max(values)

class PreparedThresholds:

    def __init__(self, edges):
        self.precision = ctx.prec
        self.edges = tuple(map(rational, edges))
        self.rows = tuple((exact(a).union(exact(b)) for a, b in zip(self.edges, self.edges[1:])))
        self.over_pi = tuple((theta / arb.pi() for theta in self.rows))

def prepared_kernel_column(T, nodes, thresholds):
    if thresholds.precision != ctx.prec:
        raise ValueError('Prepared threshold precision differs from active precision')
    T = rational(T)
    nodes = tuple(map(rational, nodes))
    count = len(nodes) - 1
    qlo = np.zeros((len(thresholds.rows), count))
    qhi = np.zeros_like(qlo)
    T_arb = exact(T)
    for i, (left, right) in enumerate(zip(nodes, nodes[1:])):
        if i == 0:
            continue
        if not 0 < left < right < T:
            raise ValueError('The kernel rectangle must avoid zero and the cutoff pole')
        tt = exact(left).union(exact(right))
        u = tt / T_arb
        z = arb.pi() * u
        h = 1 / z - z.cos() / z.sin()
        amplitude = (1 - u) / T_arb
        for j, (theta, over_pi) in enumerate(zip(thresholds.rows, thresholds.over_pi)):
            phase = theta * tt
            twice = amplitude * (phase.cos() - h * phase.sin()) + over_pi * phase.sinc()
            if not twice.is_finite():
                raise ArithmeticError('Nonfinite signed-kernel enclosure')
            qlo[j, i], qhi[j, i] = (signed.al(twice / 2), signed.au(twice / 2))
    return (qlo, qhi)

class FastSignedPreparation:

    def __init__(self, weights, theta_edges=None):
        began = time.monotonic()
        self.weights = weights
        self.theta = tuple((Q(j, 32) for j in range(-64, 65))) if theta_edges is None else tuple(map(rational, theta_edges))
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
            thresholds = PreparedThresholds(self.theta)
            for j, column in enumerate(weights.columns):
                T, s = (rational(column.T), rational(column.s))
                nodes = tuple((T * s * i / self.N for i in range(self.N + 1)))
                self.nodes.append(nodes)
                self.width_twice[:, j] = [signed._upper(2 * (right - left)) for left, right in zip(nodes, nodes[1:])]
                self.qlo[:, :, j], self.qhi[:, :, j] = prepared_kernel_column(T, nodes, thresholds)
                for r, _ in enumerate(weights.cutoff_indices):
                    record = weights.cutoff_records[r][j]
                    scalar = record['scalar_certificate']
                    signed.verify_partition(scalar)
                    parsed = PreparedGaussian(scalar)
                    correction = rational(record['endpoint_correction'])
                    if correction < 0:
                        raise ValueError('A cutoff adjustment must be nonnegative')
                    for h, (a, b) in enumerate(zip(self.theta, self.theta[1:])):
                        gamma = parsed.upper(a, b) + correction
                        self.gamma[h, r, j] = min(float(weights.cutoff_G[r, j]), signed._upper(gamma))
        finally:
            ctx.prec = previous_precision
        self.build_seconds = time.monotonic() - began

class FastPackedEnergyRectangles(FastSignedPreparation, PackedEnergyRectangles):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.prepared_reference = OrderedDict()

    def _reference(self, Llo, Lhi, dhi, taulo, tauhi):
        effective_precision = max(ctx.prec, 128)
        key = (Llo, Lhi, dhi, taulo, tauhi, ctx.prec, effective_precision)
        if key in self.reference_cache:
            self.reference_cache.move_to_end(key)
            return self.reference_cache[key]
        w = self.weights
        deleted = selected(w.tlo, w.thi, 1.0, 1.0, Lhi, dhi, taulo, tauhi, deleted=True)
        table = load_complete_energy()
        base = frozen_reference.reviewed_reference()
        previous_precision = ctx.prec
        E, lower, upper, eligible = ([], [], [], [])
        try:
            ctx.prec = effective_precision
            if effective_precision not in self.prepared_reference:
                self.prepared_reference[effective_precision] = (PreparedPrefix(table, base), tuple((PreparedGeometry(nodes, base) for nodes in self.nodes)))
                if len(self.prepared_reference) > 2:
                    self.prepared_reference.popitem(last=False)
            self.prepared_reference.move_to_end(effective_precision)
            prefix, geometries = self.prepared_reference[effective_precision]
            box = PreparedBox((Llo, Lhi), (0.0, dhi), (taulo, tauhi), base)
            for j, (nodes, geometry) in enumerate(zip(self.nodes, geometries)):
                result = fast_packed_energy_reference_bounds(nodes, deleted[:, j], energy_table=table, L=(Llo, Lhi), d=(0.0, dhi), tau=(taulo, tauhi), prepared_prefix=prefix, prepared_geometry=geometry, prepared_box=box)
                E.append(result['E_upper'])
                lower.append(result['lower_offset'])
                upper.append(result['upper_offset'])
                eligible.append(result['eligible'])
        finally:
            ctx.prec = previous_precision
        result = tuple((np.asarray(values).T for values in (E, lower, upper, eligible)))
        self.reference_cache[key] = result
        if len(self.reference_cache) > 16:
            self.reference_cache.popitem(last=False)
        return result
