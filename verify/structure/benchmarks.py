
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path


def atan_sum(x, last):
    return sum(((-1)**j*x**(2*j+1)/(2*j+1) for j in range(last+1)), Q(0))


def sqrt_bounds(v):
    lo, hi = Q(0), max(Q(1), v)
    for _ in range(48):
        mid = (lo+hi)/2
        if mid*mid <= v:
            lo = mid
        else:
            hi = mid
    assert lo*lo <= v <= hi*hi
    return lo, hi


def run():
                                                        
    pi_hi = 16*atan_sum(Q(1, 5), 10)-4*atan_sum(Q(1, 239), 1)
    normal_lo = Q(39894228, 10**8)
    assert normal_lo**2*2*pi_hi < 1
    rows = []
    for m, p, k, target in [(5, Q(9, 20), 2, Q(2, 5)),
                             (3, Q(43, 100), 1, Q(3, 8))]:
        q = 1-p
        V = m*p*q
        R = V*(p*p+q*q)
        F = sum((comb(m, j)*p**j*q**(m-j) for j in range(k+1)), Q(0))
        root_lo, root_hi = sqrt_bounds(V)
        x_lo = (m*p-k)/root_hi
        x_hi = (m*p-k)/root_lo
        assert 0 < x_lo <= x_hi < 1
                                                        
                                                        
        integral_lo = normal_lo*(x_lo-x_lo**3/6+x_lo**5/40-x_lo**7/336)
        D_lo = F-Q(1, 2)+integral_lo
        assert D_lo > 0
        C_lo = V*root_lo*D_lo/R
        assert C_lo > target
        rows.append(dict(m=m,p=str(p),owner=k,variance=str(V),third_budget=str(R),
                         cdf=str(F),ratio_lower=str(C_lo),target=str(target),
                         decimal_diagnostic=float(C_lo)))
    return dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                all_exact_checks_passed=True,pi_upper=str(pi_hi),
                gaussian_density_lower=str(normal_lo),benchmarks=rows,
                scope="Exact lower bounds at explicit admissible binomial laws; no optimization or all-owner upper bound.")


if __name__ == "__main__":
    print(json.dumps(run(),sort_keys=True,indent=2))
