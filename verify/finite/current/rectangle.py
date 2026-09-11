"""Curvature inequalities for a rectangular event in a sum of two triples."""
from fractions import Fraction

import hashlib

import json

from pathlib import Path

import flint

from flint import arb, arb_series, ctx

def rational(n, d=1):
    return arb(n) / d


def verify():
    ctx.prec = 320
    count = 1024
    cells = []
    largest_upper = None
    for j in range(count):
        midpoint = Fraction(2*j+1, 2*count)
        radius = Fraction(1, 2*count)
        interval = arb(str(midpoint), str(radius))
        # The coefficient of e in T(k+e) encloses T'(k) for every k in the cell.
        k = arb_series([interval, 1], prec=2)
        a = 3*(1+k)**2/(2*(k**3+3*k+2))
        b = (1+(1+4*a-4*a*a).sqrt())/2
        c = k*a
        expression = (b-c)*(a+c)*(c+2*b-rational(3,2))/(3*(a+b)*(b+c-1))
        derivative = expression[1]
        assert derivative < rational(-1,4), (j, derivative)
        upper = derivative.upper()
        largest_upper = upper if largest_upper is None else max(largest_upper, upper)
        cells.append({'lo': str(Fraction(j,count)), 'hi': str(Fraction(j+1,count)),
                      'derivative_lower': derivative.lower().str(80),
                      'derivative_upper': upper.str(80)})
    return {'complete': True, 'method': 'independent interval formal-series differentiation',
            'precision_bits': 320, 'python_flint': flint.__version__,
            'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'cells': count, 'exact_union': ['0','1'], 'strict_derivative_cap': '-1/4',
            'largest_upper': largest_upper.str(70),
            'intervals': cells}


