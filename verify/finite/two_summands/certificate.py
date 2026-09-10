#!/usr/bin/env python3
"""Arb bounds for the two interior closed-CDF expressions.

The default run verifies the upper interior atom, indexed by (0,1).
The lower interior expression can also be evaluated; its supplied complete
verification uses row10_certificate.py.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction

from flint import arb, ctx


@dataclass(frozen=True)
class Box:
    p1_bounds: tuple[Fraction, Fraction]
    p2_bounds: tuple[Fraction, Fraction]
    t_bounds: tuple[Fraction, Fraction]
    depth: int = 0


def interval(lo: Fraction, hi: Fraction) -> arb:
    midpoint = (lo + hi) / 2
    radius = (hi - lo) / 2
    return arb(
        f"{midpoint.numerator}/{midpoint.denominator}",
        f"{radius.numerator}/{radius.denominator}",
    )


def normal_cdf(x: arb) -> arb:
    return (1 + (x / arb(2).sqrt()).erf()) / 2


def bernoulli_variance(bounds: tuple[Fraction, Fraction]) -> arb:
    lo, hi = bounds
    endpoint_values = (lo * (1 - lo), hi * (1 - hi))
    lower = min(endpoint_values)
    upper = Fraction(1, 4) if lo <= Fraction(1, 2) <= hi else max(endpoint_values)
    return interval(lower, upper)


def split_box_axis(box: Box, axis: int) -> tuple[Box, Box]:
    bounds = (box.p1_bounds, box.p2_bounds, box.t_bounds)
    lo, hi = bounds[axis]
    midpoint = (lo + hi) / 2
    left = list(bounds)
    right = list(bounds)
    left[axis] = (lo, midpoint)
    right[axis] = (midpoint, hi)
    return (
        Box(*left, box.depth + 1),
        Box(*right, box.depth + 1),
    )


def choose_split(box: Box) -> tuple[Box, Box]:
    bounds = (box.p1_bounds, box.p2_bounds, box.t_bounds)
    widths = tuple(hi - lo for lo, hi in bounds)
    p2_touches_boundary = box.p2_bounds[0] == 0 or box.p2_bounds[1] == 1
    t_upper = box.t_bounds[1]
    if (
        p2_touches_boundary
        and t_upper < Fraction(1, 8)
        and widths[1] > t_upper**3 / 64
    ):
        return split_box_axis(box, 1)
    width_floor = max(widths) / 8
    candidates: list[tuple[float, Fraction, tuple[Box, Box]]] = []
    for axis in range(3):
        if widths[axis] < width_floor:
            continue
        children = split_box_axis(box, axis)
        score = min(
            float(gap.lower())
            for child in children
            for gap in gap_lowers(child)
        )
        candidates.append((score, widths[axis], children))
    return max(candidates, key=lambda item: (item[0], item[1]))[2]


def resource_lower(
    a1: arb,
    a2: arb,
    r1: arb,
    r2: arb,
    variance: arb,
    budget: arb,
) -> arb:
    reciprocal_denominator = a1 / (r1 * r1) + a2 / (r2 * r2)
    denominator_upper = reciprocal_denominator.upper()
    if denominator_upper == 0:
        raise AssertionError("resource denominator vanished on a nondegenerate box")
    else:
        reciprocal = 1 / denominator_upper.sqrt()

    lower = reciprocal.lower()
    if variance.lower() > 0:
        direct = budget / (variance * variance.sqrt())
        if direct.lower() > lower:
            lower = direct.lower()
    return lower


def normal_lower(numerator: arb, variance: arb) -> arb:
    variance_upper = variance.upper()
    if not variance_upper > 0:
        return arb(0)
    if variance.lower() > 0:
        return normal_cdf(numerator / variance.sqrt()).lower()
    if numerator.lower() >= 0:
        z_lower = numerator.lower() / variance_upper.sqrt()
        return normal_cdf(z_lower).lower()
    return arb(0)


def one_bit_normal_lower(bounds: tuple[Fraction, Fraction]) -> arb:
    p_upper = bounds[1]
    if p_upper >= 1:
        return arb(0)
    odds_upper = arb(p_upper.numerator) / arb(p_upper.denominator)
    odds_upper /= 1 - odds_upper
    return normal_cdf(-odds_upper.sqrt()).lower()


def gap_lowers(
    box: Box,
    *,
    constant_shift: arb | None = None,
    wrong_row01_mass: bool = False,
) -> tuple[arb, arb]:
    p1 = interval(*box.p1_bounds)
    p2 = interval(*box.p2_bounds)
    t = interval(*box.t_bounds)
    q1 = 1 - p1
    q2 = 1 - p2
    a1 = bernoulli_variance(box.p1_bounds)
    a2 = bernoulli_variance(box.p2_bounds)
    r1 = 1 - 2 * a1
    r2 = 1 - 2 * a2
    variance = a1 * t * t + a2
    budget = a1 * r1 * t * t * t + a2 * r2

    resource = resource_lower(a1, a2, r1, r2, variance, budget)
    ce = (arb(10).sqrt() + 3) / (6 * (2 * arb.pi()).sqrt())
    if constant_shift is not None:
        ce += constant_shift
    penalty = ce.lower() * resource

    x10_numerator = q1 * t - p2
    x01_numerator = -p1 * t + q2
    phi10 = normal_lower(x10_numerator, variance)
    phi01 = normal_lower(x01_numerator, variance)
    fallback10 = one_bit_normal_lower(box.p2_bounds)
    fallback01 = one_bit_normal_lower(box.p1_bounds)
    if fallback10 > phi10:
        phi10 = fallback10
    if fallback01 > phi01:
        phi01 = fallback01
    cdf10 = q2.upper()
    if wrong_row01_mass:
        cdf01 = (1 - p1 * q2).upper()
    else:
        cdf01 = (1 - p1 * p2).upper()
    return penalty + phi10 - cdf10, penalty + phi01 - cdf01


def point_box(p1: Fraction, p2: Fraction, t: Fraction) -> Box:
    return Box((p1, p1), (p2, p2), (t, t))


def direct_closed_masses(
    p1: Fraction,
    p2: Fraction,
    t: Fraction,
) -> tuple[Fraction, Fraction]:
    q1 = 1 - p1
    q2 = 1 - p2
    atoms = (
        (-p1 * t - p2, q1 * q2),
        (q1 * t - p2, p1 * q2),
        (-p1 * t + q2, q1 * p2),
        (q1 * t + q2, p1 * p2),
    )
    state10 = atoms[1][0]
    state01 = atoms[2][0]
    closed10 = sum(mass for value, mass in atoms if value <= state10)
    closed01 = sum(mass for value, mass in atoms if value <= state01)
    return closed10, closed01


def fixed_controls() -> None:
    p1 = Fraction(2, 5)
    p2 = Fraction(3, 5)
    t = Fraction(1, 2)
    closed10, closed01 = direct_closed_masses(p1, p2, t)
    if closed10 != 1 - p2 or closed01 != 1 - p1 * p2:
        raise AssertionError(
            f"four-atom control failed: closed10={closed10}, closed01={closed01}"
        )

    anchor = point_box(Fraction(131, 250), Fraction(131, 250), Fraction(1))
    anchor10, anchor01 = gap_lowers(anchor)
    if not anchor10 > arb(1) / 4 or not anchor01 > arb(3) / 100:
        raise AssertionError(f"diagonal anchor failed: {anchor10}, {anchor01}")
    _, lowered01 = gap_lowers(anchor, constant_shift=-arb(3) / 50)
    if not lowered01 < 0:
        raise AssertionError(f"lowered-constant control did not fire: {lowered01}")

    mass_anchor = point_box(Fraction(3, 5), Fraction(3, 5), Fraction(1))
    _, correct01 = gap_lowers(mass_anchor)
    _, wrong01 = gap_lowers(mass_anchor, wrong_row01_mass=True)
    if not correct01 > arb(1) / 20 or not wrong01 < 0:
        raise AssertionError(
            f"row-01 mass control failed: correct={correct01}, wrong={wrong01}"
        )
    print("fixed-controls: PASS")
    print("  direct-closed-masses", closed10, closed01)
    print("  diagonal-anchor", anchor10, anchor01)
    print("  lowered-constant-row01", lowered01)
    print("  row01-mass-control", correct01, wrong01)


def certify(max_depth: int, max_boxes: int, owners: tuple[int, ...]) -> None:
    ctx.prec = 192
    fixed_controls()
    stack = [
        Box(
            (Fraction(0), Fraction(1)),
            (Fraction(0), Fraction(1)),
            (Fraction(0), Fraction(1)),
        )
    ]
    processed = 0
    accepted = 0
    deepest = 0
    weakest = [float("inf"), float("inf")]
    weakest_data: list[tuple[Box, arb] | None] = [None, None]

    while stack:
        box = stack.pop()
        processed += 1
        if processed > max_boxes:
            print("  budget-current-box", box, gap_lowers(box))
            for pending_box in stack[-5:]:
                print("  budget-pending-box", pending_box, gap_lowers(pending_box))
            raise AssertionError(
                f"box budget exceeded: processed={processed}, pending={len(stack)}"
            )
        gaps = gap_lowers(box)
        if all(gaps[owner] > 0 for owner in owners):
            accepted += 1
            deepest = max(deepest, box.depth)
            for owner in owners:
                gap = gaps[owner]
                diagnostic = float(gap.lower())
                if diagnostic < weakest[owner]:
                    weakest[owner] = diagnostic
                    weakest_data[owner] = (box, gap)
            continue
        if box.depth >= max_depth:
            raise AssertionError(
                f"unresolved box: depth={box.depth}, p1={box.p1_bounds}, "
                f"p2={box.p2_bounds}, t={box.t_bounds}, gaps={gaps}"
            )
        stack.extend(choose_split(box))

    verdict = {
        (0,): "GLOBAL CLOSED-10 GAP-RATIO CERTIFICATE: PASS",
        (1,): "GLOBAL CLOSED-01 GAP-RATIO CERTIFICATE: PASS",
        (0, 1): "ARBITRARY-GAP TWO-BIT CROSSED OWNERS: PASS",
    }[owners]
    print(verdict)
    print("  processed-boxes", processed)
    print("  accepted-boxes", accepted)
    print("  deepest-accepted", deepest)
    for owner in owners:
        data = weakest_data[owner]
        if data is None:
            raise AssertionError(f"no accepted box for owner {owner}")
        box, gap = data
        print(f"  weakest-owner-{owner}", box, gap.lower().str(30))
    print("  coverage", f"0<=p1,p2,t<=1, V>0, owners={owners}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-depth", type=int, default=45)
    parser.add_argument("--max-boxes", type=int, default=2_000_000)
    parser.add_argument("--owner", choices=("10", "01", "both"), default="01")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_depth <= 0 or args.max_boxes <= 0:
        raise SystemExit("depth and box budget must be positive")
    owner_map = {"10": (0,), "01": (1,), "both": (0, 1)}
    certify(args.max_depth, args.max_boxes, owner_map[args.owner])


if __name__ == "__main__":
    main()
