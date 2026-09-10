#!/usr/bin/env python3
"""Exact checks of the upper interior atom and its normal-CDF lower bound."""

from __future__ import annotations

import sympy as sp


def main() -> None:
    p1, p2, t = sp.symbols(
        "p1 p2 t",
        nonnegative=True,
    )
    q1 = 1 - p1
    q2 = 1 - p2
    atoms = (
        -p1 * t - p2,
        q1 * t - p2,
        -p1 * t + q2,
        q1 * t + q2,
    )
    assert sp.expand(atoms[1] - atoms[0]) == t
    assert sp.expand(atoms[2] - atoms[1]) == 1 - t
    assert sp.expand(atoms[3] - atoms[2]) == t
    assert sp.expand((q1 * q2 + p1 * q2 + q1 * p2) - (1 - p1 * p2)) == 0
    assert sp.expand((atoms[2] - atoms[1]).subs(t, 1)) == 0
    print("ORDERED FOUR-ATOM OWNER: PASS")
    print("COLLISION t=1: PASS")

    a1 = p1 * q1
    a2 = p2 * q2
    variance = a1 * t**2 + a2
    fallback_margin = sp.factor(
        p1 * variance - q1 * (p1 * t - q2) ** 2
    )
    assert sp.expand(
        fallback_margin - q2 * (2 * p1 * q1 * t + p1 - q2)
    ) == 0
    secondary_margin = sp.expand(
        1 + t * (1 - 2 * p1)
    )
    assert sp.expand(
        secondary_margin
        - ((1 - t) + 2 * t * q1)
    ) == 0
    print("ONE-BIT NORMAL FALLBACK ALGEBRA: PASS")

    r1 = p1**2 + q1**2
    r2 = p2**2 + q2**2
    assert sp.expand(r1 - (1 - 2 * a1)) == 0
    assert sp.expand(r2 - (1 - 2 * a2)) == 0
    print("RESOURCE SHAPE IDENTITIES: PASS")


if __name__ == "__main__":
    main()
