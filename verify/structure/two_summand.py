                      


from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction

from flint import arb, ctx


@dataclass(frozen=True)
class Box:
    threshold: int
    x_bounds: tuple[Fraction, Fraction]
    p_bounds: tuple[Fraction, Fraction]
    depth: int = 0


def interval(lo: Fraction, hi: Fraction) -> arb:
    midpoint = (lo + hi) / 2
    radius = (hi - lo) / 2
    return arb(
        f"{midpoint.numerator}/{midpoint.denominator}",
        f"{radius.numerator}/{radius.denominator}",
    )


def split_box(box: Box) -> tuple[Box, Box]:
    x_width = box.x_bounds[1] - box.x_bounds[0]
    p_width = box.p_bounds[1] - box.p_bounds[0]
    if x_width >= p_width:
        lo, hi = box.x_bounds
        midpoint = (lo + hi) / 2
        return (
            Box(box.threshold, (lo, midpoint), box.p_bounds, box.depth + 1),
            Box(box.threshold, (midpoint, hi), box.p_bounds, box.depth + 1),
        )
    lo, hi = box.p_bounds
    midpoint = (lo + hi) / 2
    return (
        Box(box.threshold, box.x_bounds, (lo, midpoint), box.depth + 1),
        Box(box.threshold, box.x_bounds, (midpoint, hi), box.depth + 1),
    )


def gap_lower(box: Box, *, constant_shift: arb | None = None) -> arb | None:
    x = interval(*box.x_bounds)
    p = interval(*box.p_bounds)
    ux = x * (1 - x)
    up = p * (1 - p)
    variance = ux + up
    variance_upper = variance.upper()
    ce = (arb(10).sqrt() + 3) / (6 * (2 * arb.pi()).sqrt())
    if constant_shift is not None:
        ce += constant_shift

    if variance_upper < arb(1) / 25:
                                                          
        return ce.lower() / (2 * variance_upper.sqrt()) - 1

    variance_lower = variance.lower()
    if not variance_lower > 0:
        return None
    budget = ux * (1 - 2 * ux) + up * (1 - 2 * up)
    if budget > 0:
        budget_lower = budget.lower()
    else:
        budget_lower = variance_lower / 2
    penalty_lower = ce.lower() * budget_lower / (variance_upper * variance_upper.sqrt())

    z = (box.threshold - x - p) / variance.sqrt()
    normal = (1 + (z / arb(2).sqrt()).erf()) / 2
    normal_lower = normal.lower()
    if box.threshold == 0:
        cdf = (1 - x) * (1 - p)
    elif box.threshold == 1:
        cdf = 1 - x * p
    elif box.threshold == 2:
        cdf = arb(1)
    else:
        raise ValueError(f"unexpected threshold: {box.threshold}")
    return penalty_lower + normal_lower - cdf.upper()


def fixed_controls() -> None:
    anchor = Box(
        1,
        (Fraction(131, 250), Fraction(131, 250)),
        (Fraction(131, 250), Fraction(131, 250)),
    )
    gap = gap_lower(anchor)
    if gap is None or not gap > arb(3) / 100:
        raise AssertionError(f"two-summand anchor not reproduced: {gap}")
    false_gap = gap_lower(anchor, constant_shift=-arb(3) / 50)
    if false_gap is None or not false_gap < 0:
        raise AssertionError(f"lowered-constant control did not fire: {false_gap}")
    print("fixed-controls: PASS")
    print("  diagonal-anchor-gap", gap)
    print("  lowered-constant-control", false_gap)


def certify(max_depth: int, max_boxes: int) -> None:
    ctx.prec = 192
    fixed_controls()
    stack = [
        Box(threshold, (Fraction(0), Fraction(1)), (Fraction(0), Fraction(1)))
        for threshold in range(3)
    ]
    processed = 0
    accepted = 0
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
        gap = gap_lower(box)
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
                f"unresolved box: threshold={box.threshold}, depth={box.depth}, "
                f"x={box.x_bounds}, p={box.p_bounds}, gap={gap}"
            )
        stack.extend(split_box(box))

    if weakest is None:
        raise AssertionError("no boxes accepted")
    weakest_box, weakest_gap = weakest
    print("two-summand-certificate: PASS")
    print("  processed-boxes", processed)
    print("  accepted-boxes", accepted)
    print("  deepest-accepted", deepest)
    print("  weakest-threshold", weakest_box.threshold)
    print("  weakest-x-box", weakest_box.x_bounds)
    print("  weakest-p-box", weakest_box.p_bounds)
    print("  weakest-gap", weakest_gap)
    print("  weakest-gap-lower", weakest_gap.lower().str(30))
    print("  coverage", "0<=x,p<=1, V>0, and closed thresholds k=0,1,2")
    print("  deterministic-corners", "boundary deletion or coercive extended deficit")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-depth", type=int, default=36)
    parser.add_argument("--max-boxes", type=int, default=1_000_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_depth <= 0 or args.max_boxes <= 0:
        raise SystemExit("depth and box budget must be positive")
    certify(args.max_depth, args.max_boxes)


if __name__ == "__main__":
    main()
