from functools import lru_cache
from pathlib import Path
import math
import sys
import time
import numpy as np
HERE = Path(__file__).resolve().parent
for path in (HERE.parent, HERE.parent.parent / '26'):
    if path.exists():
        sys.path.insert(0, str(path))
from signed_gaussian import certify, verify_partition
from gaussian_batch import GaussianBatchWeights, gaussian_G

@lru_cache(maxsize=4096)
def signed_certificate(T, s):
    return certify(float(T), float(s))

def signed_G(T, s):
    return signed_certificate(T, s)['upper']

class SignedGaussianBatchWeights(GaussianBatchWeights):

    def __init__(self, Lhi, N=512, tau_cap=None):
        super().__init__(Lhi, N, tau_cap)
        self.previous_G = self.G.copy()
        self.signed_certificates = [signed_certificate(w.T, w.s) for w in self.columns]
        self.G = np.minimum(self.G, [row['upper'] for row in self.signed_certificates])

def floating_gamma(T, s, x):
    from scipy.integrate import quad
    a = T * s

    def integrand(t):
        u, z = (t / T, x * t)
        p = math.pi * u
        if abs(p) < 0.001:
            h = p / 3 + p ** 3 / 45 + 2 * p ** 5 / 945 + p ** 7 / 4725
        else:
            h = 1 / p - 1 / math.tan(p)
        sinc = 1.0 if z == 0 else math.sin(z) / z
        return math.exp(-t * t / 2) * ((1 - u) / T * (math.cos(z) - h * math.sin(z)) + x / math.pi * sinc)
    value, error = quad(integrand, 0.0, a, epsabs=2e-13, epsrel=2e-13, limit=150)
    return (value - math.erf(x / math.sqrt(2)) / 2, error)

def scalar_pilot():
    began = time.monotonic()
    rows = []
    for T, s in [(3.0, 0.6), (5.8, 0.376), (10.0, 0.3), (100.0, 0.05), (400.0, 0.0125)]:
        certificate = signed_certificate(T, s)
        structural = verify_partition(certificate)
        references = []
        for x in sorted(set([0.0, certificate['best_threshold_node'], -2.0, 2.0, -8.0, 8.0])):
            value, error = floating_gamma(T, s, x)
            if value > certificate['upper'] + abs(error) + 2e-11:
                raise AssertionError((T, s, x, value, certificate['upper']))
            references.append(dict(x=x, value=value, quadrature_error_estimate=error))
        old = gaussian_G(T, s)
        rows.append(dict(certificate=certificate, partition_check=structural, previous_G=old, coefficient_gain=old - certificate['upper'], floating_diagnostics=references))
    return dict(status='fixed-parameter scalar certificates and floating diagnostics only', rows=rows, elapsed_seconds=time.monotonic() - began)

def point_pilot(meshes=(512, 1024)):
    import stable_bounds
    from moment_fast import moment_range
    began = time.monotonic()
    original = stable_bounds.log_range_cells
    points = [(0.4, 0.195600146724, 0.092507500581, 0.376), (0.55, 0.299414822466, 0.191336229373, 0.495), (0.5945776691627666, 0.3543870444748948, 0.2470204829050391, 0.5203934540269098), (0.6, 0.397750191359, 0.280850860934, 0.54), (0.7, 0.454995299361, 0.359409523445, 0.595)]
    rows, certificates = ([], {})
    try:
        stable_bounds.log_range_cells = lambda *a, **kw: moment_range(*a, **kw, charge='tangent', baseline='original')
        for L, d, b, tau in points:
            cap = L * min(20, max(1, round(20 * tau / L))) / 20
            for N in meshes:
                weights = SignedGaussianBatchWeights(L, N, cap)
                for certificate in weights.signed_certificates:
                    key = repr((certificate['T'], certificate['s']))
                    certificates[key] = certificate
                    verify_partition(certificate)
                for eps in [2e-12, 0.0005]:
                    box = (d - eps, d + eps, b - eps, b + eps, tau - eps, tau + eps)
                    bounds = weights.envelopes(L - eps, L + eps, box)
                    if bounds is None:
                        raise AssertionError('A selected feasible point box was rejected.')
                    new_G = weights.G.copy()
                    new_values = weights.integrals(bounds[0], bounds[1], L - eps, bounds[2])
                    weights.G = weights.previous_G.copy()
                    old_values = weights.integrals(bounds[0], bounds[1], L - eps, bounds[2])
                    weights.G = new_G
                    expected_change = (weights.previous_G - new_G) / (L - eps)
                    if not np.all(np.abs(old_values - new_values - expected_change) < 2e-12):
                        raise AssertionError('A quantity other than the Gaussian coefficient changed.')
                    rows.append(dict(L=L, d=d, b=b, tau=tau, N=N, halfwidth=eps, previous_G=weights.previous_G.tolist(), signed_G=new_G.tolist(), previous_columns=old_values.tolist(), signed_columns=new_values.tolist(), previous_upper=float(np.min(old_values)), signed_upper=float(np.min(new_values))))
    finally:
        stable_bounds.log_range_cells = original
    return dict(status='interval enclosures around the sampled parameter points', rows=rows, coefficient_certificates=certificates, elapsed_seconds=time.monotonic() - began)
