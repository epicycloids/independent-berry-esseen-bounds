from functools import lru_cache
from full_b_paths import HERE
import full_b_cover as tree
import full_b_maximum_cover_engine as parent
from adaptive_b_consumer import AdaptiveBMaximumMixin
from cached_maximum_signed import CachedMaximumVarianceRectangles
from fast_maximum_cover_engine import FastMaximumVarianceWeights, CUTOFFS, verify_partition
import stable_cover
KIND = 'adaptive_b_fast_maximum_fixed085_dual5_cutoff_v1'
PRODUCT = parent.PRODUCT + '; selective complete b partitions, depth at most six, at most 127 upper evaluations, worst leaf first; exact prepared signed-cell cache'
CERTIFICATES = {}
CUTOFF_RECORDS = {}
INSTALLED = False

class AdaptiveBMaximumWeights(AdaptiveBMaximumMixin, FastMaximumVarianceWeights):
    rectangle_class = CachedMaximumVarianceRectangles

@lru_cache(maxsize=1)
def validation():
    from adaptive_b_checks import run
    return dict(inherited=parent.validation(), adaptive_partition=run())

def install():
    global INSTALLED
    parent.install()
    tree.KINDS[KIND] = (1,)
    tree.COUPLING[KIND] = False
    stable_cover.stable_checks = validation
    if INSTALLED:
        return
    previous = stable_cover.weight_factory

    def factory(hi, N, kind, target=None):
        if kind != KIND:
            return previous(hi, N, kind, target)
        cache = {}

        def choose(tauhi):
            if not cache:
                weight = AdaptiveBMaximumWeights(hi, N, hi * 17 / 20, cutoff_fractions=CUTOFFS, target=target)
                for row in weight.cutoff_records:
                    for record in row:
                        scalar = record['scalar_certificate']
                        key = repr((scalar['T'], scalar['s']))
                        if key not in CERTIFICATES:
                            verify_partition(scalar)
                            CERTIFICATES[key] = scalar
                        cutoff_key = repr((scalar['T'], scalar['s'], record['k'], record['N']))
                        CUTOFF_RECORDS[cutoff_key] = {**{k: v for k, v in record.items() if k != 'scalar_certificate'}, 'scalar_certificate_key': key}
                cache[0] = weight
            return cache[0]
        return choose
    stable_cover.weight_factory = factory
    INSTALLED = True

def run_section(start, end, step, target, N, pending=None):
    install()
    CERTIFICATES.clear()
    CUTOFF_RECORDS.clear()
    result = stable_cover.run_section(start, end, step, target, N, pending, KIND)
    result.update(product_variant=PRODUCT, gaussian_certificates=dict(CERTIFICATES), cutoff_records=dict(CUTOFF_RECORDS), split_axes=tree.SPLIT_AXES, b_policy=tree.B_POLICY, b_partitions='selective depth<=6, evaluations<=127', coupling=False, threshold_grid='inherited default 1/32 on [-2,2]')
    return result
