"""Contact equations for a triple/pair lower event."""
from fractions import Fraction as Q

import hashlib

import json

from pathlib import Path

import flint

from flint import arb, ctx

CASES = ("negative", "positive")

DOMAIN = [["0", "1"], ["0", "1"], ["0", "2"]]

BASE = ((0, 0), (0, 0), (0, 0))

WIDTHS = (Q(1), Q(1), Q(2))

def point(q):
    q = Q(q)
    return arb(q.numerator)/q.denominator


def intervals(box):
    out, endpoints = [], []
    for j, (n, level) in enumerate(box):
        lo = WIDTHS[j]*Q(n, 1 << level)
        hi = WIDTHS[j]*Q(n+1, 1 << level)
        mid, rad = (lo+hi)/2, (hi-lo)/2
        out.append(arb(f"{mid.numerator}/{mid.denominator}",
                       f"{rad.numerator}/{rad.denominator}"))
        endpoints.append((lo, hi))
    return out, endpoints


def chart_a(k, case):
    if case == "negative":
        return (2+k)/4
    return point(Q(3, 2))*(k*k+2*k+1)/(k*k*k+3*k+2)


def chart_b(a):
    return (1+(1+4*a-4*a*a).sqrt())/2


def gap_height(a, b):
    return (b-a)*(a+b+4*a*b)/2


def independently_evaluate(case, box):
    (k, middle, odds), ends = intervals(box)
    a0, a1 = chart_a(point(ends[0][0]), case), chart_a(point(ends[0][1]), case)
    b0, b1 = chart_b(a0), chart_b(a1)
    a, b = a0.union(a1), b0.union(b1)
    D = gap_height(a0, b0).union(gap_height(a1, b1))
    c = (k-1)/2 if case == "negative" else a*k
    span = a+b
    high = (a*(1-middle)-middle*c)/span
    low = 1-middle-high
    if case == "negative":
        H = 3*(a+c)*(1-a+c)
    else:
        H = 3*(b-c)*(b+c-1)
    sigma = middle*H
    p, q = odds/(1+odds), 1/(1+odds)
    square_mass_sum = p*p+q*q
    # sigma=3*q*x*(square_mass_sum*x-1), solved for diameter x.
    diameter = (1+(1+4*square_mass_sum*sigma/(3*q)).sqrt())/(2*square_mass_sum)
    d, e = p*diameter, q*diameter
    triple_var = low*a*a+middle*c*c+high*b*b
    bern_var = q*d*d+p*e*e
    variance = triple_var+bern_var
    bern_resource = q*d*d*d+p*e*e*e
    expanded_height = (d*d*d+3*d*e*e+2*e*e*e
                       -point(Q(3, 2))*diameter*diameter)
    contact = (low+middle)*D-q*expanded_height
    positive_cdf = (b*b*(2*b-point(Q(3, 2)))
                    -point(Q(3, 2))*triple_var-bern_resource)
    values = [
        ("positive_minimum_at_least_one", "nonnegative", e-1),
        ("owner_gap", "positive", diameter-a-c),
        ("triple_low_mass", "positive", low),
        ("triple_middle_mass", "positive", middle),
        ("triple_high_mass", "positive", high),
        ("owner_slope", "positive", sigma),
        ("variance", "positive", variance),
        ("normal_cdf_positive", "positive", positive_cdf),
        ("height_equation", "zero", contact),
    ]
    # q<=1 on the exact original odds domain. With D>=1, admissible
    # positive variance <=1/4 would force C=q*V**(3/2)/D<=1/8<49/120.
    # This guard needs no square root of the raw moment enclosure.
    if D >= 1:
        values.append(("small_variance_ratio_collar", "positive", variance-point(Q(1, 4))))
    if variance > 0:
        root_var = variance.sqrt()
        scaled_ratio = q*variance*root_var
        values.extend([
            ("global_ratio_upper_room", "positive", point(Q(227, 500))*D-scaled_ratio),
            ("global_ratio_lower_room", "positive", scaled_ratio-point(Q(49, 120))*D),
        ])
        offset = c-d
        true_density = (-offset*offset/(2*variance)).exp()/(2*arb.pi()).sqrt()
        density_contact = sigma*q*root_var-D*true_density
        values.append(("unsquared_normal_density_equation", "zero", density_contact))
    return values


def rejection(case, box):
    for name, kind, value in independently_evaluate(case, box):
        if ((kind == "positive" and value <= 0)
            or (kind == "nonnegative" and value < 0)
            or (kind == "zero" and (value > 0 or value < 0))):
            return name, value.str(60)
    return None


def children(box):
    axis = max(range(3), key=lambda j: WIDTHS[j]/(1 << box[j][1]))
    n, level = box[axis]
    left, right = list(box), list(box)
    left[axis], right[axis] = (2*n, level+1), (2*n+1, level+1)
    return tuple(left), tuple(right)


def decode(path):
    box = BASE
    for bit in path:
        pair = children(box)
        box = pair[int(bit)]
    return box


