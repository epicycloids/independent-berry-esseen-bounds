"""Contact equations for a triple/pair upper-gap collision."""
from fractions import Fraction as Q

import hashlib

import json

from pathlib import Path

import flint

from flint import arb, ctx

CASES = ("negative", "positive")

DOMAIN = [["0", "1"] for _ in range(4)]

BASE = ((0, 0), (0, 0), (0, 0), (0, 0))

WIDTHS = (Q(1), Q(1), Q(1), Q(1))

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


def other_minimum(a):
    return (1+(1+4*a-4*a*a).sqrt())/2


def independently_evaluate(family, case, box):
    (first, second, third, middle), ends = intervals(box)
    if family == "221":
        a0 = chart_a(point(ends[0][0]), case)
        a1 = chart_a(point(ends[0][1]), case)
        a = a0.union(a1)
        B = other_minimum(a0).union(other_minimum(a1))
        c = (first-1)/2 if case == "negative" else a*first
        b = second*B
        p, q = third, 1-third
        diameter = b-c
        high = (a*(1-middle)-middle*c)/(a+b)
        low = 1-middle-high
        beta = B*(B-1)
        Hc = (3*(a+c)*(1-a+c) if case == "negative"
              else 3*(B-c)*(B+c-1))
        Hb = 3*(B-b)*(B+b-1)
        sigma = middle*Hc+high*Hb
        abs_c3 = -c*c*c if case == "negative" else c*c*c
        # Difference g(c)-g(b), using the actual middle contact.
        level = (abs_c3-b*b*b+point(Q(3, 2))*(b*b-c*c)
                 +3*beta*(b-c))
        factor = p
        common = [
            ("positive_high_support", "positive", b),
            ("positive_collision_gap", "positive", diameter),
            ("positive_level_difference", "positive", level),
            ("high_contact_slope", "nonnegative", Hb),
        ]
    else:
        b = 1+point(Q(3, 2))*first
        a = point(Q(7, 2))*second
        c = -a+a*third if case == "negative" else b*third
        high = (a*(1-middle)-middle*c)/(a+b)
        low = 1-middle-high
        beta = b*b-b
        constant = 2*b*b*b-point(Q(3, 2))*b*b
        abs_c3 = -c*c*c if case == "negative" else c*c*c
        level = a*a*a-point(Q(3, 2))*a*a+3*beta*a+constant
        value_middle = abs_c3-point(Q(3, 2))*c*c-3*beta*c+constant
        Hl = 3*(a*a-a+beta)
        Hc = 3*(c*c+c+beta) if case == "negative" else 3*(beta+c-c*c)
        sigma = low*Hl+middle*Hc
        diameter = a+c
        common = [
            ("positive_low_support_distance", "positive", a),
            ("positive_collision_gap", "positive", diameter),
            ("positive_upper_triple_gap", "positive", b-c),
            ("low_contact_slope", "nonnegative", Hl),
            ("middle_contact_slope", "nonnegative", Hc),
            ("contact_height_lower_collar", "positive", level-point(Q(1, 2))),
            ("contact_height_upper_collar", "positive", 20-level),
        ]
        if not level > 0:
            return common
        q = value_middle/level
        p = 1-q
        factor = point(1)
    d, e = p*diameter, q*diameter
    triple_var = low*a*a+middle*c*c+high*b*b
    bern_var = q*d*d+p*e*e
    variance = triple_var+bern_var
    triple_resource = low*a*a*a+middle*abs_c3+high*b*b*b
    bern_resource = q*d*d*d+p*e*e*e
    resource = triple_resource+bern_resource
    signed_second = -q*d*d+p*e*e
    fprime_high_bern = 3*(e*e-e-signed_second)-sigma
    selected_mass = middle if family == "221" else low
    selected_H = Hc if family == "221" else Hl
    support_equation = -selected_mass*selected_H-p*fprime_high_bern
    # Integral of the Bernoulli derivative between its two endpoints.
    integrated = (e*e*e-d*d*d-point(Q(3, 2))*(e*e-d*d)
                  -3*signed_second*diameter-sigma*diameter)
    lost_mass = high if family == "221" else middle
    height_equation = factor*integrated+lost_mass*level
    offset = c+e if family == "221" else c-d
    scale_equation = variance-resource-sigma*offset/3
    squared_ratio_num = factor*factor*variance*variance*variance
    # Obtain the normal-CDF positivity test from its actual value F-C*R.
    event_probability = 1-p*high if family == "221" else low+middle*q
    common.extend([
        ("triple_low_mass", "positive", low),
        ("triple_middle_mass", "positive", middle),
        ("triple_high_mass", "positive", high),
        ("bernoulli_high_mass", "positive", p),
        ("bernoulli_low_mass", "positive", q),
        ("positive_owner_slope", "positive", sigma),
        ("positive_variance", "positive", variance),
        ("triple_variance_fraction_collar", "positive", 200*triple_var-13*variance),
        ("bernoulli_variance_fraction_collar", "positive", 200*bern_var-13*variance),
        ("normal_cdf_positive", "positive", level*event_probability-factor*resource),
        ("ratio_upper_room", "positive", point(Q(227, 500)**2)*level*level-squared_ratio_num),
        ("ratio_lower_room", "positive", squared_ratio_num-point(Q(49, 120)**2)*level*level),
        ("direct_support_gradient_equation", "zero", support_equation),
        ("integrated_bernoulli_height_equation", "zero", height_equation),
        ("scale_normalization_equation", "zero", scale_equation),
    ])
    common.extend([
        ("variance_lower_collar", "positive", variance-point(Q(25, 64))),
        ("variance_upper_collar", "positive", 4-variance),
    ])
    if variance > 0:
        z = offset/variance.sqrt()
        density = (-z*z/2).exp()/(2*arb.pi()).sqrt()
        squared_density = (sigma*sigma*factor*factor*variance
                           -level*level*(-offset*offset/variance).exp()/(2*arb.pi()))
        Phi = (-z/point(2).sqrt()).erfc()/2
        cdf_equation = level*(Phi-event_probability)+factor*resource
        common.extend([
            ("squared_normal_density_equation", "zero", squared_density),
            ("erfc_normal_cdf_equation", "zero", cdf_equation),
        ])
        gamma_denominator = (a+b)*(b-c) if case == "negative" else (a+b)*(a+c)
        if gamma_denominator > 0:
            gamma = (a-c-b+2*b*b*b/gamma_denominator if case == "negative"
                     else b+c-a+2*a*a*a/gamma_denominator)
            phi_third = (3*z-z*z*z)*density
            original_variance_soc = (12*factor*variance*gamma
                                     -15*factor*resource-level*phi_third)
            common.append(("original_variance_hessian", "nonpositive", original_variance_soc))
    return common


def rejection(family, case, box):
    for name, kind, value in independently_evaluate(family, case, box):
        if ((kind == "positive" and value <= 0)
            or (kind == "nonnegative" and value < 0)
            or (kind == "nonpositive" and value > 0)
            or (kind == "zero" and (value > 0 or value < 0))):
            return name, value.str(60)
    return None


def children(box):
    axis = max(range(4), key=lambda j: WIDTHS[j]/(1 << box[j][1]))
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


