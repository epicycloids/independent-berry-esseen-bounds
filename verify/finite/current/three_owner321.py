"""Scalar exclusion inequalities for a triangular event with three threshold representations."""
from fractions import Fraction as F

from pathlib import Path

import hashlib

import json

from flint import arb, arb_series, ctx

def exact(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def original(u, z, coefficient, slack, witness):
    # Directly construct the upper event probability in the new coordinate.
    if witness == 'bounded':
        event = 1 - (u*u - z*u + z*z*u*u + exact(slack))/2
    elif witness == 'plain':
        event = arb(1)
    elif witness == 'cantelli':
        event = 1 / (1 + z*z)
    else:
        raise AssertionError('Unknown proof witness')
    normal_cdf = (-z / arb(2).sqrt()).erfc()/2
    normal_density = (-z*z/2).exp() / (2*arb.pi()).sqrt()
    return (15*(event-normal_cdf) + (3*z-z*z*z)*normal_density
            - coefficient*u**4 - exact(F(49,20))/u)


