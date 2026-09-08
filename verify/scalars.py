"""Short Arb enclosures used at the analytic endpoints of the proof."""
from fractions import Fraction
from flint import arb, ctx


def rational(value):
    q = Fraction(value)
    return arb(q.numerator) / q.denominator


def scalar_bounds(lower_endpoint, upper_endpoint):
    previous = ctx.prec
    ctx.prec = 192
    try:
        pi = arb.pi()
        ce = (arb(10).sqrt() + 3) / (6 * (2*pi).sqrt())
        def h(x):
            return x*x + 2*x*x.sin() + 6*(x.cos()-1)
        lo, hi = Fraction('3.99'), Fraction('4')
        assert h(rational(lo)) < 0 and h(rational(hi)) > 0
        for _ in range(100):
            mid = (lo+hi)/2
            sign = h(rational(mid))
            if sign < 0:
                lo = mid
            elif sign > 0:
                hi = mid
            else:
                raise ArithmeticError('Unresolved sign in kappa root isolation')
        theta = rational(lo).union(rational(hi))
        kappa = (theta-theta.sin())/(3*theta*theta)
        assert kappa < rational('0.099162') and 3.99 < theta < 4
        # The paper proves uniqueness of the Cantelli stationary point.
        def derivative(x):
            return (-x*x/2).exp()/(2*pi).sqrt() - 2*x/(1+x*x)**2
        lo, hi = Fraction('0.1'), Fraction('0.5')
        assert derivative(rational(lo)) > 0 and derivative(rational(hi)) < 0
        for _ in range(100):
            mid = (lo+hi)/2
            sign = derivative(rational(mid))
            if sign > 0:
                lo = mid
            elif sign < 0:
                hi = mid
            else:
                raise ArithmeticError('Unresolved sign in Cantelli root isolation')
        x = rational(lo).union(rational(hi))
        variance = (1+(x/arb(2).sqrt()).erf())/2 - x*x/(1+x*x)
        start, end = rational(lower_endpoint), rational(upper_endpoint)
        assert 0 < start <= rational('0.01') and end > start
        small = ce + rational('0.3413')*start.root(3)
        large = variance/end
        assert variance < rational('0.540936541549')
        assert small < rational('0.471750509536')
        assert large < rational('0.470379601347')
        assert small < rational('0.474999998') and large < rational('0.474999998')
        return {'arb_precision_bits': ctx.prec, 'lower_bound': str(ce),
                'theta': str(theta), 'kappa': str(kappa),
                'cantelli_stationary_point': str(x), 'variance_upper': str(variance),
                'small_L_upper': str(small), 'large_L_upper': str(large)}
    finally:
        ctx.prec = previous
