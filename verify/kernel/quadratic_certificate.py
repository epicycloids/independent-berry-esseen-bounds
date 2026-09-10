from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
import sys
import time
HERE = Path(__file__).resolve().parent
if (HERE.parent / 'odd_refinement').is_dir():
    sys.path.insert(0, str(HERE.parent / 'odd_refinement'))
from scalar_majorant_certificate import components, add, bernstein_lower

def certificate():
    began = time.monotonic()
    ps = components(24)
    for ij, value in [((0, 1), -F(1, 6)), ((1, 1), -F(1, 36)), ((1, 0), -F(1, 24)), ((2, 0), -F(1, 144))]:
        add(ps[0], ij, value)
    for ij, value in [((1, 0), -F(1, 6)), ((2, 0), -F(1, 36))]:
        add(ps[2], ij, value)
    B25 = F(49, 1764) * 25 ** 2 * 196 ** 25 / F(factorial(48))
    assert F(196 * 26 ** 2, 25 ** 2 * 50 * 49) < F(1, 10)
    tail = F(5, 18) * B25
    result = {}
    for name, p in zip(('F1', 'F0', 'G'), ps):
        lower = bernstein_lower(p, (F(0), F(4), F(0), F(49)))
        assert lower - tail > F(1, 250), (name, lower, tail)
        result[name] = dict(minimum_Bernstein_coefficient=str(lower), minimum_float=float(lower), absolute_tail_upper=str(tail), absolute_tail_float=float(tail), certified_component_lower='1/250')
    return dict(status='exact rational scalar certificate passed', order=24, rectangle=dict(x_squared=[0, 4], z_squared=[0, 49]), components=result, spatial_tail_ratio_upper='16/25', elapsed_seconds=time.monotonic() - began)
