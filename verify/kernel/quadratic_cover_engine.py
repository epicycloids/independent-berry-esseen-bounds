from functools import lru_cache
from further_bounds_pilot import CutoffBatchWeights, CUTOFFS, maxwell_cells, LatestWeights, certificate, quick_checks
from dual_log_bounds import dual_range, warmup, ORIGINAL
from signed_gaussian import verify_partition
import stable_bounds
import stable_cover
KIND = 'quadratic_maxwell_dual5_cutoff'
PRODUCT = 'direct whole-log dual prices .75,1,1.5,2,3 with original baseline; quadratic x_squared_over_4 odd majorant and least concave envelope; signed Gaussian fixed-cutoff portfolio'
PRICES = (0.75, 1.0, 1.5, 2.0, 3.0)
INHERITED_CHECKS = stable_bounds.stable_checks
INSTALLED = False
CERTIFICATES = {}
CUTOFF_RECORDS = {}

def selected(*args, **kwargs):
    return dual_range(*args, **kwargs, baseline='original', prices=PRICES)

class QuadraticCutoffWeights(CutoffBatchWeights):
    prefactor = staticmethod(maxwell_cells)

@lru_cache(maxsize=1)
def validation():
    old = stable_bounds.log_range_cells
    stable_bounds.log_range_cells = ORIGINAL
    try:
        inherited = INHERITED_CHECKS()
        scalar = certificate()
        envelope = quick_checks()
        stable_bounds.log_range_cells = selected
        finite_laws = stable_bounds.finite_law_checks(weights_class=LatestWeights, extra=True)
    finally:
        stable_bounds.log_range_cells = old
    return dict(inherited=inherited, scalar=scalar, envelope=envelope, finite_laws=finite_laws, dual_support=warmup())

def install():
    global INSTALLED
    stable_bounds.log_range_cells = selected
    stable_cover.stable_checks = validation
    if INSTALLED:
        return
    previous = stable_cover.weight_factory

    def factory(hi, N, kind, target=None):
        if kind != KIND:
            return previous(hi, N, kind, target)
        cache = {}

        def choose(tauhi):
            bucket = min(20, max(1, round(20 * tauhi / hi)))
            if bucket not in cache:
                weight = QuadraticCutoffWeights(hi, N, hi * bucket / 20, cutoff_fractions=CUTOFFS)
                for row in weight.cutoff_records:
                    for record in row:
                        scalar = record['scalar_certificate']
                        key = repr((scalar['T'], scalar['s']))
                        if key not in CERTIFICATES:
                            verify_partition(scalar)
                            CERTIFICATES[key] = scalar
                        cutoff_key = repr((scalar['T'], scalar['s'], record['k'], record['N']))
                        CUTOFF_RECORDS[cutoff_key] = {**{k: v for k, v in record.items() if k != 'scalar_certificate'}, 'scalar_certificate_key': key}
                cache[bucket] = weight
            return cache[bucket]
        return choose
    stable_cover.weight_factory = factory
    INSTALLED = True
