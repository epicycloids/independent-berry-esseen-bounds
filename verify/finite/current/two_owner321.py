"""Scalar exclusion inequalities for a triangular event with two threshold representations."""
from fractions import Fraction as F

from pathlib import Path

import hashlib

import json

from flint import arb, arb_series, ctx

C0 = F(49, 120)

V0 = F(2401, 24964)

CASES = {
    'small_g': (40, F(0), F(-5, 4), F(3, 5), 32, 16),
    'large_g_negative': (22, F(0), F(-5, 4), F(-3, 4), 16, 16),
    'large_g_slack': (22, F(7, 200), F(-3, 4), F(3, 5), 48, 64),
}

def a(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def mesh(lo, hi, index, count):
    return lo + (hi-lo)*index/count, lo + (hi-lo)*(index+1)/count


def room(coordinates, mode, witness, coefficient=0, slack=F(0)):
    if mode == 'adjacent':
        r, u, z = coordinates
        bound = (32-8*r)*u**4/r**2 + 3*a(C0)*(1+r)/u
        if witness == 'bounded':
            event = 1-(u*u*(1+z*z)-z*u)/(r*(1+r))
        else:
            assert witness == 'plain'
            event = arb(1)
    else:
        u, z = coordinates
        bound = coefficient*u**4 + 6*a(C0)/u
        if witness == 'bounded':
            event = 1-(u*u*(1+z*z)-z*u+a(slack))/2
        elif witness == 'cantelli':
            event = 1/(1+z*z)
        else:
            assert witness == 'plain'
            event = arb(1)
    cdf = (-z/arb(2).sqrt()).erfc()/2
    density = (-z*z/2).exp()/(2*arb.pi()).sqrt()
    return 15*(event-cdf)+(3*z-z**3)*density-bound


