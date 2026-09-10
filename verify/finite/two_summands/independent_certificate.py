#!/usr/bin/env python3
"""Separate Arb verification of the upper interior inequality for two Bernoulli summands."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction

from flint import arb, ctx


@dataclass(frozen=True)
class Box:
    p1: tuple[Fraction, Fraction]
    p2: tuple[Fraction, Fraction]
    t: tuple[Fraction, Fraction]
    depth: int = 0


def ball(bounds: tuple[Fraction, Fraction]) -> arb:
    lo, hi = bounds
    middle = (lo + hi) / 2
    radius = (hi - lo) / 2
    return arb(
        f"{middle.numerator}/{middle.denominator}",
        f"{radius.numerator}/{radius.denominator}",
    )


def variance_hull(bounds: tuple[Fraction, Fraction]) -> arb:
    lo, hi = bounds
    values = (lo * (1 - lo), hi * (1 - hi))
    lower = min(values)
    upper = Fraction(1, 4) if lo <= Fraction(1, 2) <= hi else max(values)
    return ball((lower, upper))


def phi(x: arb) -> arb:
    return (1 + (x / arb(2).sqrt()).erf()) / 2


def fallback(p1_bounds: tuple[Fraction, Fraction]) -> arb:
    upper = p1_bounds[1]
    if upper == 1:
        return arb(0)
    probability = arb(upper.numerator) / upper.denominator
    odds = probability / (1 - probability)
    return phi(-odds.sqrt()).lower()


def gap(box: Box, *, wrong_mass: bool = False) -> arb:
    p1 = ball(box.p1)
    p2 = ball(box.p2)
    t = ball(box.t)
    q2 = 1 - p2
    a1 = variance_hull(box.p1)
    a2 = variance_hull(box.p2)
    r1 = 1 - 2 * a1
    r2 = 1 - 2 * a2
    variance = a1 * t * t + a2
    budget = a1 * r1 * t * t * t + a2 * r2

    reciprocal_sum = a1 / (r1 * r1) + a2 / (r2 * r2)
    reciprocal_upper = reciprocal_sum.upper()
    if not reciprocal_upper > 0:
        raise AssertionError("reciprocal resource vanished")
    resource = (1 / reciprocal_upper.sqrt()).lower()
    if variance.lower() > 0:
        direct = budget / (variance * variance.sqrt())
        if direct.lower() > resource:
            resource = direct.lower()

    numerator = -p1 * t + q2
    if variance.lower() > 0:
        gaussian = phi(numerator / variance.sqrt()).lower()
    elif numerator.lower() >= 0 and variance.upper() > 0:
        standardized = numerator.lower() / variance.upper().sqrt()
        gaussian = phi(standardized).lower()
    else:
        gaussian = arb(0)
    one_bit = fallback(box.p1)
    if one_bit > gaussian:
        gaussian = one_bit

    ce = (arb(10).sqrt() + 3) / (6 * (2 * arb.pi()).sqrt())
    if wrong_mass:
        owned = (1 - p1 * q2).upper()
    else:
        owned = (1 - p1 * p2).upper()
    return ce.lower() * resource + gaussian - owned


def split(box: Box) -> tuple[Box, Box]:
    bounds = (box.p1, box.p2, box.t)
    widths = [hi - lo for lo, hi in bounds]
    axis = max(range(3), key=lambda index: (widths[index], -index))
    lo, hi = bounds[axis]
    middle = (lo + hi) / 2
    left = list(bounds)
    right = list(bounds)
    left[axis] = (lo, middle)
    right[axis] = (middle, hi)
    return Box(*left, box.depth + 1), Box(*right, box.depth + 1)


def collision_control() -> None:
    for p1, p2 in (
        (Fraction(2, 5), Fraction(3, 5)),
        (Fraction(1, 3), Fraction(4, 7)),
        (Fraction(7, 11), Fraction(5, 13)),
    ):
        q1 = 1 - p1
        q2 = 1 - p2
        atoms = (
            (-p1 - p2, q1 * q2),
            (q1 - p2, p1 * q2),
            (-p1 + q2, q1 * p2),
            (q1 + q2, p1 * p2),
        )
        collision = atoms[2][0]
        closed = sum(mass for value, mass in atoms if value <= collision)
        assert atoms[1][0] == collision
        assert closed == 1 - p1 * p2
    print("COLLISION AND OWNED MASS: PASS")


def certify(
    *,
    precision: int,
    max_depth: int,
    max_boxes: int,
    wrong_mass: bool,
) -> None:
    ctx.prec = precision
    collision_control()
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
    weakest: arb | None = None
    while stack:
        box = stack.pop()
        processed += 1
        if processed > max_boxes:
            raise AssertionError(
                f"box budget exceeded with {len(stack) + 1} pending"
            )
        lower = gap(box, wrong_mass=wrong_mass)
        if lower > 0:
            accepted += 1
            deepest = max(deepest, box.depth)
            if weakest is None or lower < weakest:
                weakest = lower
            continue
        if box.depth >= max_depth:
            raise AssertionError(f"unresolved box: {box}, gap={lower}")
        stack.extend(split(box))
    if weakest is None:
        raise AssertionError("no accepted box")
    print("INDEPENDENT CLOSED-01 CERTIFICATE: PASS")
    print("  processed", processed)
    print("  accepted", accepted)
    print("  deepest", deepest)
    print("  weakest", weakest.lower().str(30))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--precision", type=int, default=192)
    parser.add_argument("--max-depth", type=int, default=60)
    parser.add_argument("--max-boxes", type=int, default=2_000_000)
    parser.add_argument("--wrong-mass", action="store_true")
    args = parser.parse_args()
    certify(
        precision=args.precision,
        max_depth=args.max_depth,
        max_boxes=args.max_boxes,
        wrong_mass=args.wrong_mass,
    )


if __name__ == "__main__":
    main()
