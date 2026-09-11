"""Contact equations for a triple/pair triangular event with a unique lower representation."""
from fractions import Fraction as Q

import hashlib

import json

from pathlib import Path

import flint

from flint import arb, ctx

ROOT_DOMAIN = [["1/2", "3/4"], ["1/2", "3/2"], ["1/2", "1"]]

BASE = ((0, 0), (0, 0), (0, 0))

WIDTH = (Q(1, 4), Q(1), Q(1, 2))

def point(q):
    q = Q(q)
    return arb(q.numerator)/q.denominator


def support_box(box):
    ans = []
    for j, (n, k) in enumerate(box):
        lo = Q(1, 2)+WIDTH[j]*Q(n, 1 << k)
        hi = Q(1, 2)+WIDTH[j]*Q(n+1, 1 << k)
        center = (lo+hi)/2
        radius = (hi-lo)/2
        ans.append(arb(f"{center.numerator}/{center.denominator}",
                       f"{radius.numerator}/{radius.denominator}"))
    return ans


def independently_evaluate(box):
    a, w, d = support_box(box)
    b = (1+(1+4*a*(1-a)).sqrt())/2
    h = a+b
    # Minimum-level difference uses a^2+b^2=a+b.
    height0 = (b-a)*(h+4*a*b)/2
    height1 = w*w*(w+3*a-point(Q(3, 2)))
    scale = height0+height1
    # Recover endpoints directly from the mass ratio, not via p,q,H_B.
    e = d*height0/height1
    diameter = d+e
    resource_slope = (d*d+e*e)/diameter
    low = d*(1-resource_slope)/(w*(w+2*a-1))
    high = (a+low*w)/h
    middle = 1-low-high
    # Raw support second moments, independent of the producer's reduced
    # variance formula ab+low*w*(a+b+w).
    triple_var = low*(a+w)*(a+w)+middle*a*a+high*b*b
    bern_var = d*e
    variance = triple_var+bern_var
    bern_resource = d*e*resource_slope
    positive_cdf = (b*b*(2*b-point(Q(3, 2)))
                    -point(Q(3, 2))*triple_var-bern_resource)
    # Expanded physical support-height difference, with no M(p) polynomial.
    height_difference = (point(Q(3, 2))*diameter*diameter
                         -e*e*e-3*d*d*e-2*d*d*d)
    contact = middle*scale-height_difference
    ratio_room = point(Q(227, 500)**2)*scale*scale-variance*variance*variance
    return {
        "bern_slope": 1-resource_slope,
        "owner_gap_lower": diameter-w,
        "owner_gap_upper": h+w-diameter,
        "triple_low_mass": low,
        "triple_middle_mass": middle,
        "triple_high_mass": high,
        "variance": variance,
        "normal_cdf_positive": positive_cdf,
        "global_ratio_room": ratio_room,
        "contact_equation": contact,
    }


def independent_rejection(box):
    for name, value in independently_evaluate(box).items():
        if name == "contact_equation":
            if value > 0 or value < 0:
                return name, value.str(60)
        elif value <= 0:
            return name, value.str(60)
    return None


def children(box):
    axis = max(range(3), key=lambda j: WIDTH[j]/(1 << box[j][1]))
    n, k = box[axis]
    left, right = list(box), list(box)
    left[axis], right[axis] = (2*n, k+1), (2*n+1, k+1)
    return tuple(left), tuple(right)


def decode_path(path):
    box = BASE
    for bit in path:
        left, right = children(box)
        box = right if bit == "1" else left
    return box


