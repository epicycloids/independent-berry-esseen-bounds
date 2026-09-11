from functools import lru_cache
from fractions import Fraction as Q
CONSTANT_UPPER = 0.17724986805182083

@lru_cache(maxsize=1)
def constant_upper():
    from flint import arb, ctx
    from interval_bounds import au
    previous = ctx.prec
    try:
        ctx.prec = 128
        computed = au(arb(1) / 5 - arb(2).sqrt().erfc() / 2)
        if computed > CONSTANT_UPPER or not arb(CONSTANT_UPPER) > arb(1) / 5 - arb(2).sqrt().erfc() / 2:
            raise ArithmeticError('The pinned outside constant is not an outward upper')
        return CONSTANT_UPPER
    finally:
        ctx.prec = previous

def outside_upper(Llo):
    if isinstance(Llo, bool) or not Q(Llo) > 0:
        raise ValueError('A positive exact L lower endpoint is required')
    from flint import arb, ctx
    from interval_bounds import au
    previous = ctx.prec
    try:
        ctx.prec = max(previous, 128)
        exact = Q(Llo)
        return au(arb(constant_upper()) / (arb(exact.numerator) / arb(exact.denominator)))
    finally:
        ctx.prec = previous

def rectangle_class():
    from cached_maximum_signed import CachedMaximumVarianceRectangles

    class CantelliNormalTailRectangles(CachedMaximumVarianceRectangles):

        def score(self, Llo, Lhi, box, envelopes, *, detail=False):
            if Q(self.theta[0]) != -2 or Q(self.theta[-1]) != 2:
                raise ValueError('The proven outside bound requires threshold domain [-2,2]')
            old = super().score(Llo, Lhi, box, envelopes, detail=True)
            if old is None:
                return None
            bound = outside_upper(Llo)
            signed = max(old['inside_upper'], bound)
            answer = min(old['old_upper'], signed)
            if not detail:
                return answer
            return {**old, 'upper': answer, 'signed_upper': signed, 'outside_upper': bound, 'old_cantelli_outside_upper': old['outside_upper'], 'outside_constant_upper': constant_upper(), 'status': 'Signed interior bound with the Cantelli-minus-normal exterior bound'}
    return CantelliNormalTailRectangles

def score_box(weight, Llo, Lhi, box, reference_box):
    from full_b_consumer import FullBMaximumMixin
    old_envelopes = super(FullBMaximumMixin, weight).envelopes(Llo, Lhi, box)
    if old_envelopes is None:
        return None
    old = weight.integrate(*old_envelopes[:2], Llo, old_envelopes[2])
    if weight.target is not None and old < weight.target:
        return old
    outside = outside_upper(Llo)
    if outside < old and (weight.target is None or outside < weight.target):
        if weight.rectangles is None:
            if weight.rectangle_class.__name__ != 'CantelliNormalTailRectangles':
                raise ValueError('The sharper guard requires the matching outside-score consumer')
            weight.rectangles = weight.rectangle_class(weight)
        if weight.rectangles.__class__.__name__ != 'CantelliNormalTailRectangles':
            raise ValueError('The outside guard requires CantelliNormalTailRectangles')
        return min(old, weight.rectangles.score(Llo, Lhi, reference_box, old_envelopes))
    return old
