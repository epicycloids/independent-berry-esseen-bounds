from functools import lru_cache
if __package__:
    from .full_b_paths import HERE
else:
    from full_b_paths import HERE
from fast_maximum_cover_engine import FastMaximumVarianceWeights, FastMaximumVarianceRectangles, validation as inherited_validation, selected, CUTOFFS, verify_partition
from full_b_consumer import FullBMaximumMixin, check_helper_source
from full_b_cover import KIND, KINDS, COUPLING, full_b_cover_band, SPLIT_AXES, B_POLICY
import stable_bounds
import stable_cover
PRODUCT = 'full feasible attained-coordinate third-moment interval at every D/tau node; closed b-partition score portfolio with maximum over every complete partition; complete centered disk prefactor and unchanged zero-bias fallback; prepared packed-energy source and attained-maximum Jensen signed reference; fixed smoothing cap .85*Lhi and signed Gaussian cutoff portfolio; closed D/tau split tree with no B splits, priorities (1,2)'
INSTALLED = False
CERTIFICATES = {}
CUTOFF_RECORDS = {}

class FullBMaximumWeights(FullBMaximumMixin, FastMaximumVarianceWeights):
    rectangle_class = FastMaximumVarianceRectangles

    def __init__(self, *args, b_partitions=(1,), coupling=False, **kwargs):
        b_partitions = tuple(b_partitions)
        if b_partitions not in KINDS.values():
            raise ValueError('Require a cumulative 1/2/4/8 b portfolio')
        self.b_partitions = b_partitions
        if type(coupling) is not bool:
            raise ValueError('coupling must be boolean')
        self.coupling = coupling
        super().__init__(*args, **kwargs)

def product(kind):
    return PRODUCT + ('; optional finite-root joint CF after the entire cheap portfolio' if COUPLING[kind] else '; inherited CF and signed reference only; coupling disabled')

@lru_cache(maxsize=1)
def validation():
    from full_b_cover_checks import run as structural_checks
    return dict(inherited=inherited_validation(), coupled_cf_helper_sha256=check_helper_source(), full_b_structure=structural_checks(), numerical_contract='frozen reviewed whole-frequency/whole-moment CF helper; prepared packed-energy/maximum reference score inherited unchanged')

def install():
    global INSTALLED
    check_helper_source()
    stable_bounds.log_range_cells = selected
    stable_cover.stable_checks = validation
    stable_cover.cover_band = full_b_cover_band
    if INSTALLED:
        return
    previous = stable_cover.weight_factory

    def factory(hi, N, kind, target=None):
        if kind not in KINDS:
            return previous(hi, N, kind, target)
        cache = {}

        def choose(tauhi):
            if 17 not in cache:
                weight = FullBMaximumWeights(hi, N, hi * 17 / 20, cutoff_fractions=CUTOFFS, target=target, b_partitions=KINDS[kind], coupling=COUPLING[kind])
                for row in weight.cutoff_records:
                    for record in row:
                        scalar = record['scalar_certificate']
                        key = repr((scalar['T'], scalar['s']))
                        if key not in CERTIFICATES:
                            verify_partition(scalar)
                            CERTIFICATES[key] = scalar
                        cutoff_key = repr((scalar['T'], scalar['s'], record['k'], record['N']))
                        CUTOFF_RECORDS[cutoff_key] = {**{k: v for k, v in record.items() if k != 'scalar_certificate'}, 'scalar_certificate_key': key}
                cache[17] = weight
            return cache[17]
        return choose
    stable_cover.weight_factory = factory
    INSTALLED = True

def run_section(start, end, step, target, N, pending=None, kind=KIND):
    if kind not in KINDS:
        raise ValueError('Unknown full-b source kind')
    install()
    CERTIFICATES.clear()
    CUTOFF_RECORDS.clear()
    result = stable_cover.run_section(start, end, step, target, N, pending, kind)
    result.update(product_variant=product(kind), gaussian_certificates=dict(CERTIFICATES), cutoff_records=dict(CUTOFF_RECORDS), split_axes=SPLIT_AXES, b_policy=B_POLICY, b_partitions=list(KINDS[kind]), coupling=COUPLING[kind], threshold_grid='inherited default 1/32 on [-2,2]')
    return result
