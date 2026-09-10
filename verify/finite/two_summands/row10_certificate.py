#!/usr/bin/env python3
"""Arb interval verification of the lower interior inequality, parametrized
by the two variance contributions.
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
    angle_bounds: tuple[Fraction, Fraction]
    depth: int = 0


def rational(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def normal_cdf(x: arb) -> arb:
    return (1 + (x / arb(2).sqrt()).erf()) / 2


def shape_at(p_value: Fraction) -> arb:
    p = rational(p_value)
    q = 1 - p
    variance = p * q
    if not variance > 0:
        raise AssertionError("shape evaluated at a deterministic probability")
    return (p * p + q * q) / variance.sqrt()


def shape_lower(bounds: tuple[Fraction, Fraction]) -> arb:
    lo, hi = bounds
    half = Fraction(1, 2)
    if lo <= half <= hi:
        return arb(1)
    nearest = hi if hi < half else lo
    return shape_at(nearest).lower()


def row10_lower(box: Box, *, constant_shift: arb | None = None) -> arb:
    p1_lo, p1_hi = box.p1_bounds
    p2_lo, p2_hi = box.p2_bounds
    angle_lo, angle_hi = box.angle_bounds
    theta_lo = arb.pi() * rational(angle_lo) / 2
    theta_hi = arb.pi() * rational(angle_hi) / 2
    alpha_lower = theta_hi.cos().lower()
    beta_lower = theta_lo.sin().lower()
    beta_upper = theta_hi.sin().upper()
    if alpha_lower < 0:
        alpha_lower = arb(0)
    if beta_lower < 0:
        beta_lower = arb(0)

    resource = alpha_lower**3 * shape_lower(box.p1_bounds)
    resource += beta_lower**3 * shape_lower(box.p2_bounds)
    ce = (arb(10).sqrt() + 3) / (6 * (2 * arb.pi()).sqrt())
    if constant_shift is not None:
        ce += constant_shift
    penalty = ce.lower() * resource

    if p2_hi == 1:
        normal = arb(0)
    else:
        if p1_hi == 1:
            positive_phase = arb(0)
        else:
            positive_phase = alpha_lower * (
                (1 - rational(p1_hi)) / rational(p1_hi)
            ).sqrt().lower()
        negative_phase = beta_upper * (
            rational(p2_hi) / (1 - rational(p2_hi))
        ).sqrt().upper()
        normal = normal_cdf(positive_phase - negative_phase).lower()

    cdf_upper = 1 - rational(p2_lo)
    return penalty + normal - cdf_upper


def split_box(box: Box) -> tuple[Box, Box]:
    bounds = (box.p1_bounds, box.p2_bounds, box.angle_bounds)
    widths = tuple(hi - lo for lo, hi in bounds)
    axis = max(range(3), key=widths.__getitem__)
    lo, hi = bounds[axis]
    midpoint = (lo + hi) / 2
    left = list(bounds)
    right = list(bounds)
    left[axis] = (lo, midpoint)
    right[axis] = (midpoint, hi)
    return Box(*left, box.depth + 1), Box(*right, box.depth + 1)


def point_box(p1: Fraction, p2: Fraction, angle: Fraction) -> Box:
    return Box((p1, p1), (p2, p2), (angle, angle))


def fixed_controls() -> None:
    anchor = point_box(Fraction(1, 2), Fraction(1, 2), Fraction(1, 2))
    gap = row10_lower(anchor)
    if not gap > arb(1) / 4:
        raise AssertionError(f"row-10 anchor failed: {gap}")

    lower_anchor = point_box(Fraction(2, 5), Fraction(2, 5), Fraction(1))
    lower_gap = row10_lower(lower_anchor)
    if not lower_gap > arb(1) / 25:
        raise AssertionError(f"one-bit lower-support anchor failed: {lower_gap}")
    shifted = row10_lower(lower_anchor, constant_shift=-arb(1) / 20)
    if not shifted < 0:
        raise AssertionError(f"lowered-constant control did not fire: {shifted}")
    print("fixed-controls: PASS")
    print("  interior-row10-anchor", gap)
    print("  one-bit-lower-anchor", lower_gap)
    print("  lowered-constant-control", shifted)


def certify(max_depth: int, max_boxes: int) -> None:
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
    weakest = float("inf")
    weakest_data: tuple[Box, arb] | None = None

    while stack:
        box = stack.pop()
        processed += 1
        if processed > max_boxes:
            raise AssertionError(
                f"box budget exceeded: processed={processed}, pending={len(stack)}"
            )
        gap = row10_lower(box)
        if gap > 0:
            accepted += 1
            deepest = max(deepest, box.depth)
            diagnostic = float(gap.lower())
            if diagnostic < weakest:
                weakest = diagnostic
                weakest_data = (box, gap)
            continue
        if box.depth >= max_depth:
            raise AssertionError(
                f"unresolved box: depth={box.depth}, p1={box.p1_bounds}, "
                f"p2={box.p2_bounds}, angle={box.angle_bounds}, gap={gap}"
            )
        stack.extend(split_box(box))

    if weakest_data is None:
        raise AssertionError("no box accepted")
    weakest_box, weakest_gap = weakest_data
    print("GLOBAL CLOSED-10 VARIANCE-SHARE CERTIFICATE: PASS")
    print("  processed-boxes", processed)
    print("  accepted-boxes", accepted)
    print("  deepest-accepted", deepest)
    print("  weakest-box", weakest_box)
    print("  weakest-gap", weakest_gap.lower().str(30))
    print("  coverage", "0<=p1,p2,theta/(pi/2)<=1")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-depth", type=int, default=42)
    parser.add_argument("--max-boxes", type=int, default=1_000_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_depth <= 0 or args.max_boxes <= 0:
        raise SystemExit("depth and box budget must be positive")
    certify(args.max_depth, args.max_boxes)


if __name__ == "__main__":
    main()
