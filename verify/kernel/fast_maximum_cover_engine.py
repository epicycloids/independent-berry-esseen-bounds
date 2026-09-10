from functools import lru_cache
from flint import arb
from scalar_table_cover_engine import ScalarCutoffWeights, ScalarLawWeights, selected, CUTOFFS, INHERITED_CHECKS, ORIGINAL, certificate, quick_checks, warmup
from disk_table_bounds import table_cells, load_table
from fast_maximum_variance_signed import FastMaximumVarianceRectangles
from complete_energy_signed import load_complete_energy
from balanced_cover import balanced_cover_band
from interval_bounds import al, au
from signed_gaussian import verify_partition
import stable_bounds
import stable_cover
KIND = 'fast_maximum_variance_fixed085_disk_reference_balanced_dual5_cutoff'
PRODUCT = 'direct whole-log dual prices .75,1,1.5,2,3 with original baseline; complete direct centered disk table and previous scalar/quadratic fallbacks; energy majorants through zero at price9/8 with finite coefficient packing; attained-maximum Jensen reference; prepared exact arithmetic; fixed smoothing cap .85*Lhi; signed reference rectangles on [-2,2] and Cantelli outside; signed Gaussian fixed-cutoff portfolio; coordinate split priorities (1,1/4,2)'
INSTALLED = False
CERTIFICATES = {}
CUTOFF_RECORDS = {}

class FastMaximumVarianceWeights(ScalarCutoffWeights):
    prefactor = staticmethod(table_cells)

    def __init__(self, *args, target=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.target = None if target is None else al(arb(target))
        self.rectangles = None

    def box(self, Llo, Lhi, *box):
        envelopes = self.envelopes(Llo, Lhi, box)
        if envelopes is None:
            return None
        old = self.integrate(*envelopes[:2], Llo, envelopes[2])
        outside = au(1 / (5 * arb(Llo)))
        if self.target is not None and (old < self.target or outside >= self.target):
            return old
        if old > 0.48 or outside >= old:
            return old
        if self.rectangles is None:
            self.rectangles = FastMaximumVarianceRectangles(self)
        return min(old, self.rectangles.score(Llo, Lhi, box, envelopes))

class FastMaximumVarianceLawWeights(ScalarLawWeights):
    prefactor = staticmethod(table_cells)

@lru_cache(maxsize=1)
def validation():
    old = stable_bounds.log_range_cells
    stable_bounds.log_range_cells = ORIGINAL
    try:
        inherited = INHERITED_CHECKS()
        scalar = certificate()
        envelope = quick_checks()
        endpoints, constants, prices = load_table()
        energy = load_complete_energy()
        stable_bounds.log_range_cells = selected
        finite_laws = stable_bounds.finite_law_checks(weights_class=FastMaximumVarianceLawWeights, extra=True)
    finally:
        stable_bounds.log_range_cells = old
    return dict(inherited=inherited, scalar=scalar, envelope=envelope, table=dict(rows=len(endpoints), last_endpoint=float(endpoints[-1]), prices=prices.tolist()), finite_laws=finite_laws, dual_support=warmup(), signed_reference='closed threshold rectangles and reviewed conditional reference recurrence', small_energy_rectangles=energy['small_frequency_prefix']['verified_rectangles'])

def install():
    global INSTALLED
    stable_bounds.log_range_cells = selected
    stable_cover.stable_checks = validation
    stable_cover.cover_band = balanced_cover_band
    if INSTALLED:
        return
    previous = stable_cover.weight_factory

    def factory(hi, N, kind, target=None):
        if kind != KIND:
            return previous(hi, N, kind, target)
        cache = {}

        def choose(tauhi):
            bucket = 17
            if bucket not in cache:
                weight = FastMaximumVarianceWeights(hi, N, hi * bucket / 20, cutoff_fractions=CUTOFFS, target=target)
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

def run_section(start, end, step, target, N, pending=None):
    install()
    CERTIFICATES.clear()
    CUTOFF_RECORDS.clear()
    result = stable_cover.run_section(start, end, step, target, N, pending, KIND)
    result['product_variant'] = PRODUCT
    result['gaussian_certificates'] = dict(CERTIFICATES)
    result['cutoff_records'] = dict(CUTOFF_RECORDS)
    return result
