from fractions import Fraction
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
if (HERE.parent / 'quadratic_cover_engine.py').is_file():
    sys.path.insert(0, str(HERE.parent))
from flint import arb
import quadratic_cover_engine as original
from interval_bounds import al, au
from signed_gaussian import verify_partition
import stable_cover
KIND = 'quadratic_maxwell_dual5_cutoff_lazy_ceil20'
PRODUCT = original.PRODUCT + '; quadratic smoothing portfolio first, then a cap rounded upward to a multiple of Lhi/20'
INSTALLED = False

def ceiling_bucket(hi, tauhi):
    if not hi > 0:
        raise ValueError('Require positive Lhi.')
    ratio = 20 * Fraction.from_float(float(tauhi)) / Fraction.from_float(float(hi))
    return min(20, max(1, -(-ratio.numerator // ratio.denominator)))

def ceiling_cap(hi, bucket):
    if not hi > 0 or not 1 <= bucket <= 20:
        raise ValueError('Require positive Lhi and a bucket from 1 through 20.')
    return au(arb(float(hi)) * int(bucket) / 20)

def collect(weights):
    for row in weights.cutoff_records:
        for record in row:
            scalar = record['scalar_certificate']
            key = repr((scalar['T'], scalar['s']))
            if key not in original.CERTIFICATES:
                verify_partition(scalar)
                original.CERTIFICATES[key] = scalar
            cutoff_key = repr((scalar['T'], scalar['s'], record['k'], record['N']))
            original.CUTOFF_RECORDS[cutoff_key] = {**{k: v for k, v in record.items() if k != 'scalar_certificate'}, 'scalar_certificate_key': key}

def install():
    global INSTALLED
    original.install()
    if INSTALLED:
        return
    previous = stable_cover.weight_factory

    def factory(hi, N, kind, target=None):
        if kind != KIND:
            return previous(hi, N, kind, target)
        old_choose = previous(hi, N, original.KIND, target)
        limit = None if target is None else al(arb(target))
        ceilings, portfolios = ({}, {})

        def extra(bucket):
            if bucket not in ceilings:
                weight = original.QuadraticCutoffWeights(hi, N, ceiling_cap(hi, bucket), cutoff_fractions=original.CUTOFFS)
                collect(weight)
                ceilings[bucket] = weight
            return ceilings[bucket]

        class Portfolio:

            def __init__(self, old_weight, bucket):
                self.old_weight, self.bucket = (old_weight, bucket)
                self.N, self.T, self.s = (N, old_weight.T, old_weight.s)

            def uniform(self, lo, high):
                return self.old_weight.uniform(lo, high)

            def box(self, *args):
                first = self.old_weight.box(*args)
                if first is None or (limit is not None and first < limit):
                    return first
                second = extra(self.bucket).box(*args)
                assert second is not None
                return min(first, second)

        def choose(tauhi):
            old_weight = old_choose(tauhi)
            bucket = ceiling_bucket(hi, tauhi)
            key = (id(old_weight), bucket)
            if key not in portfolios:
                portfolios[key] = Portfolio(old_weight, bucket)
            return portfolios[key]
        return choose
    stable_cover.weight_factory = factory
    INSTALLED = True
