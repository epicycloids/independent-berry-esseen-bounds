#!/usr/bin/env python3
"""Exact checks of atom order, collisions, deletion, and reflection for two
Bernoulli summands.

The script does not import or execute either interval verification program.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import product


Q = Fraction


def atom_law(
    p1: Fraction, p2: Fraction, h1: Fraction, h2: Fraction
) -> dict[Fraction, Fraction]:
    q1, q2 = 1 - p1, 1 - p2
    mean = h1 * p1 + h2 * p2
    grouped: dict[Fraction, Fraction] = defaultdict(Fraction)
    for a1, a2 in product((0, 1), repeat=2):
        value = h1 * a1 + h2 * a2 - mean
        mass = (p1 if a1 else q1) * (p2 if a2 else q2)
        if mass:
            grouped[value] += mass
    return dict(grouped)


def raw_atom(
    a1: int,
    a2: int,
    p1: Fraction,
    p2: Fraction,
    h1: Fraction,
    h2: Fraction,
) -> Fraction:
    return h1 * (a1 - p1) + h2 * (a2 - p2)


def cdf_pair(
    law: dict[Fraction, Fraction], threshold: Fraction
) -> tuple[Fraction, Fraction]:
    strict = sum((mass for value, mass in law.items() if value < threshold), Q(0))
    closed = sum((mass for value, mass in law.items() if value <= threshold), Q(0))
    return strict, closed


def variance_resource(
    p1: Fraction, p2: Fraction, h1: Fraction, h2: Fraction
) -> tuple[Fraction, Fraction]:
    variance = Q(0)
    resource = Q(0)
    for p, h in ((p1, h1), (p2, h2)):
        q = 1 - p
        variance += p * q * h * h
        resource += p * q * (p * p + q * q) * h**3
    return variance, resource


def check_strict_chamber() -> None:
    probabilities = (Q(1, 7), Q(1, 3), Q(1, 2), Q(4, 5))
    gaps = ((Q(1), Q(2)), (Q(2), Q(5)), (Q(3, 7), Q(11, 9)))
    for p1, p2, (h1, h2) in product(probabilities, probabilities, gaps):
        assert 0 < h1 < h2
        q1, q2 = 1 - p1, 1 - p2
        law = atom_law(p1, p2, h1, h2)
        expected = {
            (0, 0): (Q(0), q1 * q2),
            (1, 0): (q1 * q2, q2),
            (0, 1): (q2, 1 - p1 * p2),
            (1, 1): (1 - p1 * p2, Q(1)),
        }
        ordered_values = [
            raw_atom(a1, a2, p1, p2, h1, h2)
            for a1, a2 in ((0, 0), (1, 0), (0, 1), (1, 1))
        ]
        assert ordered_values == sorted(ordered_values)
        for atom, wanted in expected.items():
            threshold = raw_atom(*atom, p1, p2, h1, h2)
            assert cdf_pair(law, threshold) == wanted


def check_equal_gap_collision() -> None:
    probabilities = (Q(1, 7), Q(2, 5), Q(4, 5))
    for p1, p2, h in product(probabilities, probabilities, (Q(1), Q(7, 3))):
        q1, q2 = 1 - p1, 1 - p2
        law = atom_law(p1, p2, h, h)
        x10 = raw_atom(1, 0, p1, p2, h, h)
        x01 = raw_atom(0, 1, p1, p2, h, h)
        assert x10 == x01
        assert cdf_pair(law, x10) == (q1 * q2, 1 - p1 * p2)
        assert q2 != 1 - p1 * p2


def check_zero_gap_face() -> None:
    probabilities = (Q(1, 7), Q(2, 5), Q(4, 5))
    for p1, p2, h2 in product(probabilities, probabilities, (Q(1), Q(7, 3))):
        q2 = 1 - p2
        law = atom_law(p1, p2, Q(0), h2)
        lower = raw_atom(0, 0, p1, p2, Q(0), h2)
        upper = raw_atom(0, 1, p1, p2, Q(0), h2)
        assert len(law) == 2
        assert cdf_pair(law, lower) == (Q(0), q2)
        assert cdf_pair(law, upper) == (q2, Q(1))


def check_deterministic_faces() -> None:
    interior = (Q(1, 7), Q(2, 5), Q(4, 5))
    for deterministic, p, h1, h2 in product(
        (Q(0), Q(1)), interior, (Q(1), Q(2)), (Q(3), Q(5))
    ):
        for p1, p2 in ((deterministic, p), (p, deterministic)):
            law = atom_law(p1, p2, h1, h2)
            variance, _ = variance_resource(p1, p2, h1, h2)
            assert len(law) == 2
            assert variance > 0
    for p1, p2 in product((Q(0), Q(1)), repeat=2):
        law = atom_law(p1, p2, Q(1), Q(2))
        variance, resource = variance_resource(p1, p2, Q(1), Q(2))
        assert len(law) == 1
        assert variance == resource == 0


def check_complement_reflection() -> None:
    probabilities = (Q(1, 7), Q(1, 3), Q(4, 5))
    gaps = ((Q(1), Q(2)), (Q(2), Q(5)))
    for p1, p2, (h1, h2) in product(probabilities, probabilities, gaps):
        original = atom_law(p1, p2, h1, h2)
        reflected = atom_law(1 - p1, 1 - p2, h1, h2)
        assert reflected == {-value: mass for value, mass in original.items()}
        assert variance_resource(p1, p2, h1, h2) == variance_resource(
            1 - p1, 1 - p2, h1, h2
        )
        for threshold in original:
            strict, _ = cdf_pair(original, threshold)
            _, reflected_closed = cdf_pair(reflected, -threshold)
            assert strict == 1 - reflected_closed


def main() -> None:
    check_strict_chamber()
    check_equal_gap_collision()
    check_zero_gap_face()
    check_deterministic_faces()
    check_complement_reflection()
    print("independent frame checks: PASS")


if __name__ == "__main__":
    main()
