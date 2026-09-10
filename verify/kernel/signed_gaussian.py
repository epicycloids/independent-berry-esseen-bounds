import heapq
import math
import time
from flint import arb, acb, ctx

def upper(value):
    return math.nextafter(float(arb(value).upper()), math.inf)

def lower(value):
    return math.nextafter(float(arb(value).lower()), -math.inf)

class GaussianSupCertificate:

    def __init__(self, T, s, tolerance=2e-06, X=8.0, max_nodes=161, polynomial_tolerance=2e-12):
        self.T_value, self.s_value = (float(T), float(s))
        self.T, self.s = (arb(self.T_value), arb(self.s_value))
        if not (self.T > 0 and self.s > 0 and (self.s < 1)):
            raise ValueError('Require T>0 and 0<s<1.')
        self.X = float(X)
        if not (math.isfinite(self.X) and self.X > 0):
            raise ValueError('Require finite X>0.')
        if not (math.isfinite(float(tolerance)) and tolerance > 0):
            raise ValueError('Require a positive finite tolerance.')
        if int(max_nodes) < 3:
            raise ValueError('Require at least three threshold nodes.')
        self.tolerance = float(tolerance)
        self.max_nodes = int(max_nodes)
        self.pi, self.sqrt2 = (arb.pi(), arb(2).sqrt())
        self.a = self.T * self.s
        m = 8
        while True:
            error = self.pi * self.s ** (2 * m + 2) / (3 * (1 - self.s * self.s) * (2 * m + 2))
            if error < polynomial_tolerance:
                break
            m += 1
            if m > 512:
                raise ValueError('s is too close to one for this polynomial certifier.')
        self.degree = m
        self.polynomial_error = upper(error)
        self.coefficients = [2 * arb(2 * j + 2).zeta() / self.pi for j in range(m)]
        ea = (-self.a * self.a / 2).exp()
        erf_a = (self.a / self.sqrt2).erf()
        H = self.pi * self.s / (3 * (1 - self.s * self.s))
        curvature = (1 + H) * ((self.pi / 2).sqrt() / self.T * erf_a - self.s * ea) + ea / self.pi
        if not curvature > 0:
            raise ArithmeticError('Nonpositive curvature enclosure.')
        self.curvature_upper = upper(curvature)
        xx = arb(self.X)
        sinc_tail = min(arb(1), 16 / (self.T * self.T * xx * xx))
        self.outside_upper = upper(sinc_tail + (xx / (2 * self.sqrt2)).erfc() + (self.pi / 2).sqrt() / self.T * (self.a / self.sqrt2).erfc() + (self.a * self.a / 2).expint(1) / (2 * self.pi))
        self.central_exact = (self.pi / 2).sqrt() / self.T * erf_a - (1 - ea) / (self.T * self.T)
        self.nodes = {}
        self.evaluations = 0

    def polynomial(self, u):
        u2 = u * u
        result = acb(self.coefficients[-1])
        for coefficient in reversed(self.coefficients[:-1]):
            result = result * u2 + coefficient
        return u * result

    def value(self, x):
        x = float(x)
        if x in self.nodes:
            return self.nodes[x]
        ax = arb(x)

        def integrand(t, analytic):
            u, phase = (t / self.T, ax * t)
            return (-t * t / 2).exp() * ((1 - u) / self.T * (phase.cos() - self.polynomial(u) * phase.sin()) + ax / self.pi * phase.sinc())
        integral = acb.integral(integrand, acb(0), acb(self.a), rel_tol=2e-13, abs_tol=2e-14, eval_limit=20000, depth_limit=30, use_heap=True)
        if not (integral.is_finite() and integral.imag.contains(0)):
            raise ArithmeticError(f'Nonfinite or nonreal integration at x={x}.')
        approximate = integral.real - (ax / self.sqrt2).erf() / 2
        lo = lower(approximate - self.polynomial_error)
        hi = upper(approximate + self.polynomial_error)
        if x == 0.0:
            if not (arb(lo) <= self.central_exact and self.central_exact <= arb(hi)):
                raise ArithmeticError('Central exact formula disagrees with quadrature.')
        self.nodes[x] = (lo, hi)
        self.evaluations += 1
        return (lo, hi)

    def interval_upper(self, left, right):
        endpoint = max(self.value(left)[1], self.value(right)[1])
        return upper(arb(endpoint) + arb(self.curvature_upper) * (arb(right) - arb(left)) ** 2 / 8)

    def run(self):
        began = time.monotonic()
        for x in (-self.X, 0.0, self.X):
            self.value(x)
        heap = []
        for left, right in [(-self.X, 0.0), (0.0, self.X)]:
            heapq.heappush(heap, (-self.interval_upper(left, right), left, right))
        while len(self.nodes) < self.max_nodes:
            best_lower = max((value[0] for value in self.nodes.values()))
            current_upper = -heap[0][0]
            if current_upper <= best_lower + self.tolerance:
                break
            if self.outside_upper >= current_upper:
                break
            _, left, right = heapq.heappop(heap)
            midpoint = (left + right) / 2
            if midpoint in (left, right):
                heapq.heappush(heap, (-self.interval_upper(left, right), left, right))
                break
            self.value(midpoint)
            for l, r in [(left, midpoint), (midpoint, right)]:
                heapq.heappush(heap, (-self.interval_upper(l, r), l, r))
        intervals = sorted(((left, right, -negative) for negative, left, right in heap))
        inside_upper = max((row[2] for row in intervals))
        best_x = max(self.nodes, key=lambda x: self.nodes[x][0])
        best_lower = self.nodes[best_x][0]
        return dict(status='certified Gaussian coefficient for one fixed (T,s); no array supremum', T=self.T_value, s=self.s_value, X=self.X, upper=max(inside_upper, self.outside_upper), lower=best_lower, inside_upper=inside_upper, outside_upper=self.outside_upper, best_threshold_node=best_x, polynomial_terms=self.degree, polynomial_error=self.polynomial_error, curvature_upper=self.curvature_upper, requested_gap=self.tolerance, achieved_gap=upper(arb(max(inside_upper, self.outside_upper)) - arb(best_lower)), max_nodes=self.max_nodes, evaluations=self.evaluations, nodes=[dict(x=x, lower=lo, upper=hi) for x, (lo, hi) in sorted(self.nodes.items())], intervals=[dict(left=l, right=r, upper=v) for l, r, v in intervals], precision_bits=ctx.prec, elapsed_seconds=time.monotonic() - began)

def certify(T, s, **kwargs):
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128)
        return GaussianSupCertificate(T, s, **kwargs).run()
    finally:
        ctx.prec = previous

def verify_partition(record):
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128)
        rows = record['intervals']
        if not rows or rows[0]['left'] != -record['X'] or rows[-1]['right'] != record['X']:
            raise AssertionError('Incomplete threshold interval.')
        nodes = {row['x']: row for row in record['nodes']}
        for i, row in enumerate(rows):
            l, r = (row['left'], row['right'])
            if not l < r or (i and rows[i - 1]['right'] != l):
                raise AssertionError('Threshold partition gap or overlap.')
            value = upper(arb(max(nodes[l]['upper'], nodes[r]['upper'])) + arb(record['curvature_upper']) * (arb(r) - arb(l)) ** 2 / 8)
            if value > row['upper']:
                raise AssertionError('Understated interval bound.')
        if max(record['outside_upper'], *(row['upper'] for row in rows)) > record['upper']:
            raise AssertionError('Understated global coefficient.')
        return dict(all_passed=True, threshold_nodes=len(nodes), intervals=len(rows), scope='partition and stored interpolation bounds only')
    finally:
        ctx.prec = previous
