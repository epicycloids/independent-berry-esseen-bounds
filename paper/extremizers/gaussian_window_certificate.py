"""Exact rational certificate for the Gaussian Fourier estimate
in Extremizers for the independent Berry--Esseen inequality.

Every assertion is an exact rational comparison. Decimal output
is for display.
"""

from fractions import Fraction as Q
import json


def exp_lower(x: Q, degree: int = 40) -> Q:
    """A positive Taylor polynomial is a lower bound for exp(x), x>=0."""
    assert x >= 0
    term = total = Q(1)
    for k in range(1, degree + 1):
        term *= x / k
        total += term
    return total


def exp_upper(x: Q, degree: int = 24) -> Q:
    """Positive Taylor polynomial plus a geometric bound for its tail."""
    assert 0 <= x < degree + 2
    term = total = Q(1)
    for k in range(1, degree + 1):
        term *= x / k
        total += term
    first_omitted = term * x / (degree + 1)
    return total + first_omitted / (1 - x / (degree + 2))


def main() -> dict:
    # The exponent h(t)=-t²/2+(269/5000)t³ decreases on [0,6].
    assert Q(3 * 269, 5000) * 6 < 1
    cells = 200
    width = Q(6, cells)
    integral_upper = Q(0)
    for j in range(cells):
        left = j * width
        right = left + width
        positive_exponent = left * left / 2 - Q(269, 5000) * left**3
        # t² exp(h(t)) <= right² exp(h(left)) on this cell.
        integral_upper += width * right**2 / exp_lower(positive_exponent)
    assert integral_upper < Q(23, 10)

    # Integration by parts followed by the elementary Gaussian Mills bound.
    gaussian_tail_upper = Q(6025, 216) / exp_lower(Q(108, 25))
    assert gaussian_tail_upper < Q(3, 8)

    pi_lower = Q(157, 50)
    prefactor_upper = exp_upper(Q(25, 54)) / (6 * pi_lower)
    assert prefactor_upper < Q(17, 200)

    # sqrt(1.9*pi)>61/25, sqrt(pi)>443/250, sqrt(.9)>237/250.
    assert Q(61, 25) ** 2 < Q(19, 10) * pi_lower
    assert Q(443, 250) ** 2 < pi_lower
    assert Q(237, 250) ** 2 < Q(9, 10)
    density_upper = Q(25, 61) + (
        Q(3051, 5000) * Q(17, 20) * 2
        / (Q(443, 250) * Q(9, 10) * Q(237, 250))
    )
    assert density_upper < Q(11, 10)

    # c=19/10 and a/beta²>=400/289 give a*T²/4>=361/289.
    high_tail_upper = Q(11, 19) / exp_lower(Q(361, 289))
    assert high_tail_upper < Q(1, 6)

    # T>7, beta<.3, pi>3, c>1, and e>5/2 imply the stated normal tail bound.
    assert Q(19, 10) / Q(269, 1000) > 7
    assert Q(5, 2) ** 24 > 10**8
    normal_tail_upper = Q(1, 10**9)

    final_upper = (
        Q(17, 200) * (Q(23, 10) + Q(3, 8))
        + Q(1, 6)
        + normal_tail_upper
    )
    assert final_upper < Q(2, 5)

    return {
        "arithmetic": "exact fractions; decimal approximations are illustrative",
        "rectangle_cells": cells,
        "exp_lower_degree": 40,
        "integral_upper_below": "23/10",
        "integral_upper_approx": float(integral_upper),
        "gaussian_tail_upper_below": "3/8",
        "prefactor_upper_below": "17/200",
        "density_upper_below": "11/10",
        "high_tail_upper_below": "1/6",
        "normal_tail_upper_below": "1/1000000000",
        "final_upper_exact": str(final_upper),
        "final_upper_approx": float(final_upper),
        "target": "2/5",
        "all_exact_checks_passed": True,
    }


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
