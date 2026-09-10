#!/usr/bin/env python3
"""Arb interval verification of the scalar atom inequality using even moments.

For 2m-2 <= nu <= 2m, log-convexity bounds the normalized integral J_nu
by the geometric interpolation of q_{m-1}=J_{2m-2} and q_m=J_{2m}.
Maximizing that interpolation over positive nu gives a one-variable
majorant, including the two endpoints of the order interval.
Every accepted box has a strictly positive interval lower bound.
"""

from __future__ import annotations

from math import comb
from sys import argv

from flint import arb, ctx


DEFAULT_CELLS = 2048


def constants() -> tuple[arb, arb, arb, arb]:
    root10 = arb(10).sqrt()
    ke = (root10 + 3) / 6
    star = root10 - 3
    a = 1 / root10
    t0 = star + (1 / (5 * ke)).sqrt()
    return ke, star, a, t0


def parameter_box(index: int, cells: int) -> arb:
    """The exact rational interval [index/cells,(index+1)/cells]."""

    return arb(f"{2 * index + 1}/{2 * cells}", f"1/{2 * cells}")


def q_m(m: int, t: arb) -> arb:
    """Exact positive-coefficient polynomial for J_{2m}(t)."""

    x = t * t
    denominator = arb(4**m)
    answer = arb(0)
    for k in range(m, -1, -1):
        coefficient = arb(comb(2 * k, k) * comb(2 * (m - k), m - k)) / denominator
        answer = answer * x + coefficient
    return answer


def global_order_majorant(m: int, t: arb) -> arb:
    """Majorize the whole slab 2m-2 <= nu <= 2m.

    If r=q_m/q_{m-1}, the interpolating envelope is proportional to
    sqrt(nu) r^(nu/2).  Its global maximum over nu>0 is at
    nu_*=-1/log(r); taking the global rather than constrained maximum also
    covers the cases where nu_* lies outside the slab.
    """

    qlo = q_m(m - 1, t)
    qhi = q_m(m, t)
    ratio = qhi / qlo
    ell = -ratio.log()
    assert ell > 0
    lo = arb(2 * m - 2)
    c = 1 - t * t
    prefactor = (arb.pi() * c / (2 * arb(1).exp() * ell)).sqrt()
    return prefactor * qlo * (lo * ell / 2).exp()


def rhs(t: arb, baseline: arb | None = None) -> arb:
    ke, star, _, _ = constants()
    if baseline is None:
        baseline = arb(1)
    delta = t - star
    return baseline + ke * delta * delta


def fixed_controls() -> None:
    """Independent special-function and negative-branch controls."""

    _, star, _, _ = constants()
    nu = arb(5)
    t = arb("0.51")
    c = 1 - t * t
    g = (arb.pi() * nu * c / 2).sqrt() * c.hypgeom_2f1(-nu / 2, arb(1) / 2, 1)
    assert g > 1

    # The exact polynomial agrees with the hypergeometric value at an even
    # order, independently of the positive-coefficient evaluator.
    m = 7
    tcheck = arb("0.47")
    ccheck = 1 - tcheck * tcheck
    hyper = ccheck.hypgeom_2f1(-m, arb(1) / 2, 1)
    poly = q_m(m, tcheck)
    assert (hyper - poly).contains(0)

    near_star = arb(1000)
    cstar = 1 - star * star
    gstar = (arb.pi() * near_star * cstar / 2).sqrt() * cstar.hypgeom_2f1(
        -near_star / 2, arb(1) / 2, 1
    )
    false_gap = rhs(star, arb(99) / 100) - gstar
    assert false_gap < 0

    print("compact-controls: PASS")
    print("  fixed-G-over-one", g)
    print("  even-moment-crosscheck", hyper - poly)
    print("  lowered-baseline-control", false_gap)


def certify(cells: int) -> None:
    ctx.prec = 256
    fixed_controls()
    _, _, a, t0 = constants()

    weakest_float = float("inf")
    weakest = None
    accepted = 0
    uniform_floor = arb(19) / 1000
    for index in range(cells):
        s = parameter_box(index, cells)
        t = a + (t0 - a) * s
        right = rhs(t)
        for m in range(2, 65):
            majorant = global_order_majorant(m, t)
            gap = right - majorant
            if not gap > 0:
                raise AssertionError(
                    f"unresolved box index={index} m={m} t={t} "
                    f"rhs={right} majorant={majorant} gap={gap}"
                )
            if not gap > uniform_floor:
                raise AssertionError(
                    f"gap does not clear registered report floor 0.019: "
                    f"index={index} m={m} gap={gap}"
                )
            accepted += 1
            diagnostic = float(gap.lower())
            if diagnostic < weakest_float:
                weakest_float = diagnostic
                weakest = (index, m, t, gap, majorant, right)

    assert weakest is not None
    index, m, t, gap, majorant, right = weakest
    print("compact-certificate: PASS")
    print("  cells", cells)
    print("  slabs", 63)
    print("  accepted-box-slabs", accepted)
    print("  weakest-index", index)
    print("  weakest-m", m)
    print("  weakest-t-box", t)
    print("  weakest-majorant", majorant)
    print("  weakest-rhs", right)
    print("  weakest-certified-gap", gap)
    print("  weakest-gap-lower", gap.lower().str(30))
    print("  weakest-gap-upper", gap.upper().str(30))
    print("  uniform-certified-floor", uniform_floor)
    print("  coverage", "2 <= nu <= 128, 1/sqrt(10) <= t <= t0")


def main() -> None:
    cells = int(argv[1]) if len(argv) > 1 else DEFAULT_CELLS
    if cells <= 0:
        raise SystemExit("cell count must be positive")
    certify(cells)


if __name__ == "__main__":
    main()
