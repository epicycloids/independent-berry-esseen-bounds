from functools import lru_cache
from flint import arb
from scalar_table_cover_engine import ScalarCutoffWeights, ScalarLawWeights, selected, CUTOFFS, INHERITED_CHECKS, ORIGINAL, certificate, quick_checks, warmup
from disk_table_bounds import table_cells, load_table
from interval_signed import SignedRectangles
from interval_bounds import al, au
from signed_gaussian import verify_partition
import stable_bounds
import stable_cover
KIND = 'direct_disk_signed_reference_dual5_cutoff'
PRODUCT = 'direct whole-log dual prices .75,1,1.5,2,3 with original baseline; complete direct centered disk table and scalar and quadratic fallbacks; signed reference rectangles on [-2,2] and Cantelli outside; signed Gaussian fixed-cutoff portfolio'
INSTALLED = False
CERTIFICATES = {}
CUTOFF_RECORDS = {}

class DiskReferenceWeights(ScalarCutoffWeights):
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
            self.rectangles = SignedRectangles(self)
        return min(old, self.rectangles.score(Llo, Lhi, box, envelopes))

class DiskLawWeights(ScalarLawWeights):
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
        stable_bounds.log_range_cells = selected
        finite_laws = stable_bounds.finite_law_checks(weights_class=DiskLawWeights, extra=True)
    finally:
        stable_bounds.log_range_cells = old
    return dict(inherited=inherited, scalar=scalar, envelope=envelope, table=dict(rows=len(endpoints), last_endpoint=float(endpoints[-1]), prices=prices.tolist()), finite_laws=finite_laws, dual_support=warmup(), signed_reference='closed threshold rectangles and reviewed conditional reference recurrence')

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
                weight = DiskReferenceWeights(hi, N, hi * bucket / 20, cutoff_fractions=CUTOFFS, target=target)
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
