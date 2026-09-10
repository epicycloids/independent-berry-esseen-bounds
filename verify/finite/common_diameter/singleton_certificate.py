#!/usr/bin/env python3
"""Arb interval verification of the curvature inequality for a binomial sum
and one Bernoulli variable.

The three parametrizations cover 0<p,x<1, p(1-p)<=W<=3,
0<=f<=1, and W>=(1-f)(1-p)^2. The proof bounds three neighboring
binomial masses and uses

    z phi(z) (2+z^2) >= sqrt(2/pi) min(z,0).

Every accepted box has a strictly positive interval lower bound.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction

from flint import arb, ctx


CASE_HIGH_P = "p>=1/2"
CASE_LOW_P_HIGH_W = "p<=1/2,W>=q^2"
CASE_LOW_P_LOW_W = "p<=1/2,W<=q^2"
CASES = (CASE_HIGH_P, CASE_LOW_P_HIGH_W, CASE_LOW_P_LOW_W)
DELETION_VARIANCE_CAP = 3


@dataclass(frozen=True)
class Box:
    case: str
    bounds: tuple[tuple[Fraction, Fraction], ...]
    depth: int = 0


def interval(lo: Fraction, hi: Fraction) -> arb:
    """Return the exact outward ball with endpoints ``lo`` and ``hi``."""

    midpoint = (lo + hi) / 2
    radius = (hi - lo) / 2
    return arb(
        f"{midpoint.numerator}/{midpoint.denominator}",
        f"{radius.numerator}/{radius.denominator}",
    )


def unit_box(case: str) -> Box:
    return Box(case, ((Fraction(0), Fraction(1)),) * 4)


def split_box(box: Box) -> tuple[Box, Box]:
    widths = [hi - lo for lo, hi in box.bounds]
    coordinate = max(range(4), key=lambda index: widths[index])
    lo, hi = box.bounds[coordinate]
    midpoint = (lo + hi) / 2
    left = list(box.bounds)
    right = list(box.bounds)
    left[coordinate] = (lo, midpoint)
    right[coordinate] = (midpoint, hi)
    return (
        Box(box.case, tuple(left), box.depth + 1),
        Box(box.case, tuple(right), box.depth + 1),
    )


def physical_parameters(box: Box) -> tuple[arb, arb, arb, arb, arb, arb]:
    """Map a unit box onto one exact piece of the physical domain."""

    sp, sx, sw, sf = (interval(lo, hi) for lo, hi in box.bounds)
    x = sx
    if box.case == CASE_HIGH_P:
        p = (1 + sp) / 2
        q = 1 - p
        v = p * q
        w = v + sw * (DELETION_VARIANCE_CAP - v)
        phase = sf
        one_minus_phase = 1 - sf
        previous_numerator = q * sp + sw * (DELETION_VARIANCE_CAP - v) + phase * q * q
    elif box.case == CASE_LOW_P_HIGH_W:
        p = sp / 2
        q = 1 - p
        v = p * q
        q_squared = q * q
        w = q_squared + sw * (DELETION_VARIANCE_CAP - q_squared)
        phase = sf
        one_minus_phase = 1 - sf
        previous_numerator = (
            sw * (DELETION_VARIANCE_CAP - q_squared) + phase * q_squared
        )
    elif box.case == CASE_LOW_P_LOW_W:
        p = sp / 2
        q = 1 - p
        v = p * q
        q_squared = q * q
        w = v + sw * (q_squared - v)
        phase = ((1 - sw) * (q_squared - v) + sf * w) / q_squared
        one_minus_phase = (1 - sf) * w / q_squared
        previous_numerator = w * sf
    else:
        raise ValueError(f"unknown domain case: {box.case}")
    return p, x, w, phase, one_minus_phase, previous_numerator


def gap_enclosure(box: Box, *, penalty_multiplier: arb | None = None) -> arb | None:
    """Return an outward lower bound for the scalar gap on ``box``.

    ``None`` means that interval dependency has hidden positivity of a
    denominator; subdivision is then mandatory.
    """

    p, x, w, phase, one_minus_phase, previous_numerator = physical_parameters(box)
    q = 1 - p
    v = p * q
    pair_scale = w + v
    variance = x * (1 - x) + pair_scale
    variance_lower = variance.lower()
    variance_upper = variance.upper()
    ce = (arb(10).sqrt() + 3) / (6 * (2 * arb.pi()).sqrt())
    if variance_upper < arb(1) / 5:
        # Universally max(V,S_2)>=1/3, hence the smooth term is at least
        # 8 C_E/sqrt(3).  Also h(z)>=-4/(e sqrt(pi)) and ell<=1.
        return (
            8 * ce.lower() / arb(3).sqrt().upper()
            - 4 / (arb(1).exp().lower() * arb.pi().sqrt().lower())
            - 5 * variance_upper
        )
    if not variance_lower > 0:
        return None
    root_lower = variance_lower.sqrt()
    root_upper = variance_upper.sqrt()

    tx = 1 - 2 * x
    tp = 1 - 2 * p
    smooth_numerator = 4 * w + 2 + tx * tp
    if smooth_numerator > 1:
        smooth_numerator_lower = smooth_numerator.lower()
    else:
        # Since tx,tp are in [-1,1] and W>=0, this is an exact fallback.
        smooth_numerator_lower = arb(1)
    smooth_term_lower = 2 * ce.lower() * smooth_numerator_lower / root_upper

    # The phase parametrization is y=f(W+v)/(1-qf).  Its denominator can
    # vanish only at the deleted p=0 boundary.  There we use the weaker but
    # continuous lower bound y>=0 instead of dividing an indeterminate ball.
    phase_denominator = p + q * one_minus_phase
    phase_numerator = pair_scale - p * x - one_minus_phase * (pair_scale + q * x)
    normal_slope = (arb(2) / arb.pi()).sqrt()
    if phase_numerator > 0:
        normal_term_lower = arb(0)
    else:
        numerator_lower = phase_numerator.lower()
        if not numerator_lower < 0:
            normal_term_lower = arb(0)
        elif phase_denominator > 0:
            denominator_lower = phase_denominator.lower()
            normal_term_lower = (
                normal_slope.upper()
                * numerator_lower
                / (denominator_lower * root_lower)
            )
        else:
            normal_term_lower = -normal_slope.upper() * x.upper() / root_lower

    # Let d_j be the upper neighboring atom and f=(d_j-d_{j+1})/d_j.
    # The recurrence gives d_{j-1}=r_- d_j, where the following rational
    # expression is the exact r_- after eliminating m and k.
    previous_denominator = p * p + one_minus_phase * (w + 2 * v)

    # ``previous_numerator`` is nonnegative by the case parametrizations;
    # a ball touching zero need not certify that fact through comparison.
    def nonnegative_lower(value: arb) -> arb:
        return value.lower() if value > 0 else arb(0)

    p_lower = nonnegative_lower(p)
    one_minus_lower = nonnegative_lower(one_minus_phase)
    w_plus_lower = nonnegative_lower(w + 2 * v)
    previous_lower = p_lower * p_lower + one_minus_lower * w_plus_lower
    previous_upper = previous_denominator.upper()
    atom_lower = nonnegative_lower(previous_numerator)
    neighbor_lower = (1 + one_minus_lower) * previous_lower + atom_lower
    if neighbor_lower > 0:
        phase_upper = phase.upper()
        if not phase_upper < 1:
            phase_upper = arb(1)
        adjacent_upper = phase_upper * previous_upper / neighbor_lower
        if not adjacent_upper < 1:
            adjacent_upper = arb(1)
    else:
        # Every signed adjacent difference of a probability mass function is
        # at most one.  This handles the deleted deterministic boundary.
        adjacent_upper = arb(1)
    if penalty_multiplier is None:
        penalty_multiplier = arb(1)
    penalty_upper = 5 * penalty_multiplier.upper() * variance_upper * adjacent_upper
    return smooth_term_lower + normal_term_lower - penalty_upper


def fixed_controls() -> None:
    """Check constants, parametrizations, and fixed perturbation controls."""

    ce = (arb(10).sqrt() + 3) / (6 * (2 * arb.pi()).sqrt())
    assert ce > arb(409) / 1000
    assert ce < arb(410) / 1000
    assert ce > arb(2) / 5

    normal_floor = 4 / (arb(1).exp() * arb.pi().sqrt())
    assert normal_floor < 1
    small_variance_margin = 8 * ce / arb(3).sqrt() - normal_floor - 1
    assert small_variance_margin > arb(1) / 247
    assert arb(3).sqrt() > arb(5) / 3
    assert (arb(7) / 2).sqrt() > arb(11) / 6
    assert arb.pi() > 3
    tail_left = 8 * ce * arb(3).sqrt() - 10 / arb.pi() - normal_floor
    tail_right = 8 * ce * (arb(7) / 2).sqrt() - arb(35) / (3 * arb.pi()) - normal_floor
    assert tail_left > 1
    assert tail_right > arb(3) / 2

    # The same coarse Fourier tail estimate is deliberately insufficient at
    # W=1, so the compact certificate cannot silently lose its lower slabs.
    premature_tail = 8 * ce * (arb(3) / 2).sqrt() - arb(15) / arb.pi() - normal_floor
    assert premature_tail < 0

    # Exact physical anchor near the floating minimizer of the relaxation.
    anchor = Box(
        CASE_HIGH_P,
        (
            (Fraction(4633, 10000), Fraction(4633, 10000)),
            (Fraction(4186, 10000), Fraction(4186, 10000)),
            (Fraction(0), Fraction(0)),
            (Fraction(1), Fraction(1)),
        ),
    )
    anchor_gap = gap_enclosure(anchor)
    if anchor_gap is None or not anchor_gap > arb(45) / 100:
        raise AssertionError(f"anchor gap not reproduced: {anchor_gap}")

    # Increasing the adjacent-mass penalty by 25 percent must fail at an
    # unrelated box; this catches a reversed inequality or missing factor.
    false_gap = gap_enclosure(anchor, penalty_multiplier=arb(5) / 4)
    if false_gap is None or not false_gap < 0:
        raise AssertionError(f"inflated-penalty control did not fire: {false_gap}")

    print("fixed-controls: PASS")
    print("  C_E", ce)
    print("  anchor-gap", anchor_gap)
    print("  inflated-penalty-control", false_gap)
    print("  W>=3-tail-left-endpoint", tail_left)
    print("  W>=3-tail-right-endpoint", tail_right)
    print("  premature-W=1-tail-control", premature_tail)
    print("  V<=1/5-fallback-margin", small_variance_margin)


def certify(max_depth: int, max_boxes: int) -> None:
    ctx.prec = 192
    fixed_controls()
    stack = [unit_box(case) for case in CASES]
    accepted = 0
    processed = 0
    deepest = 0
    weakest_lower = float("inf")
    weakest: tuple[Box, arb] | None = None

    while stack:
        box = stack.pop()
        processed += 1
        if processed > max_boxes:
            raise AssertionError(
                f"box budget exceeded: processed={processed}, pending={len(stack)}"
            )
        gap = gap_enclosure(box)
        if gap is not None and gap > 0:
            accepted += 1
            deepest = max(deepest, box.depth)
            diagnostic = float(gap.lower())
            if diagnostic < weakest_lower:
                weakest_lower = diagnostic
                weakest = (box, gap)
            continue
        if box.depth >= max_depth:
            raise AssertionError(
                f"unresolved box at depth {box.depth}: case={box.case}, "
                f"bounds={box.bounds}, gap={gap}"
            )
        stack.extend(split_box(box))

    if weakest is None:
        raise AssertionError("no boxes accepted")
    weakest_box, weakest_gap = weakest
    print("singleton-compact-certificate: PASS")
    print("  cases", len(CASES))
    print("  processed-boxes", processed)
    print("  accepted-boxes", accepted)
    print("  deepest-accepted", deepest)
    print("  weakest-case", weakest_box.case)
    print("  weakest-normalized-box", weakest_box.bounds)
    print("  weakest-gap", weakest_gap)
    print("  weakest-gap-lower", weakest_gap.lower().str(30))
    print(
        "  coverage",
        "0<p<1, 0<x<1, p(1-p)<=W<=3, valid 0<f<=1",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-depth", type=int, default=36)
    parser.add_argument("--max-boxes", type=int, default=2_000_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_depth <= 0 or args.max_boxes <= 0:
        raise SystemExit("depth and box budget must be positive")
    certify(args.max_depth, args.max_boxes)


if __name__ == "__main__":
    main()
