#!/usr/bin/env python3
"""Separate exact Bernstein verification of the Bernoulli-odds inequalities."""

from __future__ import annotations

from math import comb

import sympy as sp


def bernstein_coefficients_2d(
    expression: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    degree_x: int | None = None,
    degree_y: int | None = None,
) -> list[list[sp.Expr]]:
    """Return the tensor Bernstein coefficients, with optional degree elevation."""
    polynomial = sp.Poly(expression, x, y, extension=sp.sqrt(5))
    monomial_degree_x = polynomial.degree(x)
    monomial_degree_y = polynomial.degree(y)
    target_degree_x = monomial_degree_x if degree_x is None else degree_x
    target_degree_y = monomial_degree_y if degree_y is None else degree_y
    assert target_degree_x >= monomial_degree_x
    assert target_degree_y >= monomial_degree_y

    monomial = dict(polynomial.terms())
    first_transform: dict[tuple[int, int], sp.Expr] = {}
    for p in range(target_degree_x + 1):
        for j in range(monomial_degree_y + 1):
            first_transform[p, j] = sum(
                monomial.get((i, j), sp.Rational(0))
                * sp.Rational(comb(p, i), comb(target_degree_x, i))
                for i in range(min(p, monomial_degree_x) + 1)
            )

    coefficients: list[list[sp.Expr]] = []
    for p in range(target_degree_x + 1):
        row: list[sp.Expr] = []
        for q in range(target_degree_y + 1):
            value = sum(
                first_transform[p, j]
                * sp.Rational(comb(q, j), comb(target_degree_y, j))
                for j in range(min(q, monomial_degree_y) + 1)
            )
            row.append(sp.expand(value))
        coefficients.append(row)
    return coefficients


def quadratic_parts(value: sp.Expr) -> tuple[sp.Rational, sp.Rational]:
    """Write an element of Q(sqrt(5)) uniquely as A + B sqrt(5)."""
    root_five = sp.sqrt(5)
    expanded = sp.expand(value)
    coefficient = sp.expand(expanded).coeff(root_five)
    constant = sp.expand(expanded - coefficient * root_five)
    assert constant.is_Rational
    assert coefficient.is_Rational
    return sp.Rational(constant), sp.Rational(coefficient)


def quadratic_sign(value: sp.Expr) -> int:
    """Decide the sign in Q(sqrt(5)) without numerical approximation."""
    constant, coefficient = quadratic_parts(value)
    if constant == 0:
        return int(sp.sign(coefficient))
    if coefficient == 0:
        return int(sp.sign(constant))
    if (constant > 0 and coefficient > 0) or (
        constant < 0 and coefficient < 0
    ):
        return int(sp.sign(constant))
    square_comparison = constant**2 - 5 * coefficient**2
    if square_comparison == 0:
        return 0
    if square_comparison > 0:
        return int(sp.sign(constant))
    return int(sp.sign(coefficient))


def print_table(name: str, table: list[list[sp.Expr]]) -> None:
    """Emit every coefficient of a tensor Bernstein table."""
    print(name)
    for row_index, row in enumerate(table):
        print(f"row {row_index}: " + " ".join(str(value) for value in row))


def main() -> None:
    a, k, r, t = sp.symbols("a k r t", nonnegative=True)
    x, y = sp.symbols("x y", nonnegative=True)
    root_five = sp.sqrt(5)
    crossover = 3 - root_five

    assert quadratic_sign(3 - root_five) == 1
    assert quadratic_sign(2 - root_five) == -1
    assert quadratic_sign(root_five - 2) == 1
    assert quadratic_sign(root_five - 3) == -1
    assert quadratic_sign(5 - root_five**2) == 0
    assert sp.expand(crossover**2 - 6 * crossover + 4) == 0
    assert quadratic_sign(5 * crossover / 2 - sp.Rational(16, 9)) == 1

    def g(o: sp.Expr) -> sp.Expr:
        return o * (1 + o) ** 2 / (1 + o**2) ** 2

    d = 1 + a**2 - 2 * k + k**2
    n = (
        a
        + 2 * a**2
        + a**3
        + (a**3 - a - 4) * k
        + (2 * a**2 - a + 8) * k**2
        + (a - 4) * k**3
    )
    atom_denominator = 1 + a + k
    pair_resource = n / (d**2 * atom_denominator**2)

    direct_pair = (
        g(x) + g(y)
    ) / ((1 + x) ** 2 * (1 + y) ** 2)
    reduced_pair = pair_resource.subs({a: x + y, k: x * y})
    assert sp.cancel(direct_pair - reduced_pair) == 0

    merge_resource = sp.factor(pair_resource.subs(k, 0))
    equal_product = a**2 / 4
    equal_resource = sp.factor(pair_resource.subs(k, equal_product))

    small_gap = sp.factor(merge_resource - pair_resource)
    small_numerator, small_denominator = sp.together(small_gap).as_numer_denom()
    small_quotient, small_remainder = sp.div(small_numerator, k, k)
    assert small_remainder == 0
    assert sp.factor(small_denominator) == (
        (1 + a**2) ** 2
        * (1 + a + k) ** 2
        * (1 + a**2 - 2 * k + k**2) ** 2
    )
    failed_split_witness = sp.factor(
        small_gap.subs({a: sp.Rational(4, 5), k: sp.Rational(4, 25)})
    )
    assert failed_split_witness == -sp.Rational(408000, 69272329)

    small_box = sp.expand(
        small_quotient.subs(
            {a: crossover * r, k: crossover**2 * r**2 * t / 4},
            simultaneous=True,
        )
    )
    small_polynomial = sp.Poly(small_box, r, t, extension=root_five)
    small_monomial_degree = (
        small_polynomial.degree(r),
        small_polynomial.degree(t),
    )
    assert small_monomial_degree == (11, 5)
    small_table = bernstein_coefficients_2d(
        small_box,
        r,
        t,
        degree_x=60,
        degree_y=20,
    )
    assert len(small_table) == 61
    assert all(len(row) == 21 for row in small_table)
    small_signs = [
        (row_index, column_index, quadratic_sign(value))
        for row_index, row in enumerate(small_table)
        for column_index, value in enumerate(row)
    ]
    small_negative = [entry for entry in small_signs if entry[2] < 0]
    assert not small_negative, small_negative[:5]
    small_zero = [entry[:2] for entry in small_signs if entry[2] == 0]
    assert small_zero == [(60, 20)], small_zero

    z_profile = 1 / atom_denominator**2
    z_merge = 1 / (1 + a) ** 2
    z_equal = 1 / (1 + a + equal_product) ** 2
    theta = sp.factor((z_profile - z_equal) / (z_merge - z_equal))
    chord_gap = sp.factor(
        theta * merge_resource + (1 - theta) * equal_resource - pair_resource
    )
    chord_numerator, chord_denominator = sp.together(chord_gap).as_numer_denom()
    chord_factor = k * (a**2 - 4 * k)
    chord_quotient, chord_remainder = sp.div(chord_numerator, chord_factor, k)
    assert chord_remainder == 0
    assert sp.factor(chord_denominator) == (
        (1 + a**2) ** 2
        * (4 + a**2) ** 2
        * (1 + a + k) ** 2
        * (8 + 8 * a + a**2)
        * (1 + a**2 - 2 * k + k**2) ** 2
    )
    large_a = crossover + (sp.Rational(16, 9) - crossover) * r
    large_box = (
        chord_quotient.subs(
            {a: large_a, k: large_a**2 * t / 4},
            simultaneous=True,
        )
    )
    large_polynomial = sp.Poly(large_box, r, t, extension=root_five)
    large_degree = (
        large_polynomial.degree(r),
        large_polynomial.degree(t),
    )
    assert large_degree == (13, 4), large_degree
    large_table = bernstein_coefficients_2d(large_box, r, t)
    assert len(large_table) == 14
    assert all(len(row) == 5 for row in large_table)
    large_signs = [quadratic_sign(value) for row in large_table for value in row]
    assert all(sign > 0 for sign in large_signs)

    s = sp.symbols("s", positive=True)

    def equal_block(block_size: int) -> sp.Expr:
        return s / (
            (1 + (s / block_size) ** 2) ** 2
            * (1 + s / block_size) ** (2 * (block_size - 1))
        )

    equal_two = equal_block(2)
    comparison_one = sp.factor(equal_two - equal_block(1))
    comparison_three = sp.factor(equal_two - equal_block(3))
    comparison_four = sp.factor(equal_two - equal_block(4))
    expected_one = -(
        s**2
        * (s**2 - 6 * s + 4)
        * (s**3 + 10 * s**2 + 4 * s + 16)
    ) / ((s + 2) ** 2 * (s**2 + 1) ** 2 * (s**2 + 4) ** 2)
    expected_three = (
        s**2
        * (8 * s**3 - 33 * s**2 - 18 * s + 108)
        * (8 * s**4 + 129 * s**3 + 306 * s**2 + 756 * s + 1296)
    ) / (
        (s + 2) ** 2
        * (s + 3) ** 4
        * (s**2 + 4) ** 2
        * (s**2 + 9) ** 2
    )
    expected_four = (
        64
        * s**2
        * (s**4 + 12 * s**3 - 64 * s**2 + 256)
        * (s**5 + 12 * s**4 + 192 * s**3 + 512 * s**2 + 1280 * s + 2048)
    ) / (
        (s + 2) ** 2
        * (s + 4) ** 6
        * (s**2 + 4) ** 2
        * (s**2 + 16) ** 2
    )
    assert sp.cancel(comparison_one - expected_one) == 0
    assert sp.cancel(comparison_three - expected_three) == 0
    assert sp.cancel(comparison_four - expected_four) == 0
    assert sp.factor((6 * s - s**2 - 4).subs(s, 1)) == 1
    assert sp.factor((6 * s - s**2 - 4).subs(s, sp.Rational(16, 9))) > 0
    assert sp.factor((8 * s**3 - 33 * s**2 - 18 * s + 108).subs(s, sp.Rational(16, 9))) == sp.Rational(12140, 729)
    assert sp.factor((s**4 + 12 * s**3 - 64 * s**2 + 256).subs(s, sp.Rational(16, 9))) == sp.Rational(860416, 6561)

    o = sp.symbols("o", positive=True)
    elasticity = sp.factor(o * sp.diff(g(o), o) / g(o))
    elasticity_gap = (
        9 * o**3 + 17 * o**2 - 7 * o + 1
    ) / (4 * (1 + o) * (1 + o**2))
    assert sp.cancel(sp.Rational(5, 4) - elasticity - elasticity_gap) == 0
    assert sp.expand(
        17 * o**2
        - 7 * o
        + 1
        - 17 * (o - sp.Rational(7, 34)) ** 2
        - sp.Rational(19, 68)
    ) == 0
    u, v = sp.symbols("u v", nonnegative=True)
    odds_fraction_gap = sp.factor(
        u / (1 + u) + v / (1 + v) - (u + v) / (1 + u + v)
    )
    assert odds_fraction_gap == u * v * (2 + u + v) / (
        (1 + u) * (1 + v) * (1 + u + v)
    )
    boundary_envelope = sp.factor(equal_two.subs(s, sp.Rational(16, 9)))
    assert boundary_envelope == sp.Rational(944784, 6076225)
    boundary_resource_margin = sp.factor(
        sp.Rational(2401, 14400) - boundary_envelope
    )
    assert boundary_resource_margin == sp.Rational(7873013, 699981120)

    print("PAIR PROFILE IDENTITIES: PASS")
    print("algebraic sign controls: PASS")
    print("failed 4/5 split witness:", failed_split_witness)
    print("small monomial degree:", small_monomial_degree)
    print("small elevated Bernstein degree: (60, 20)")
    print("small coefficient count: 1281")
    print("small negative coefficient count: 0")
    print("small zero coefficients:", small_zero)
    print_table("SMALL PAIR BERNSTEIN TABLE", small_table)
    print("large minimal Bernstein degree:", large_degree)
    print("large coefficient count: 70")
    print("large nonpositive coefficient count: 0")
    print_table("LARGE PAIR BERNSTEIN TABLE", large_table)
    print("EQUAL-BLOCK COMPARISONS d=1,3,4: PASS")
    print("RADIAL ELASTICITY BOUND: PASS")
    print("boundary resource margin:", boundary_resource_margin)
    print("ALL-DIMENSIONAL PAIR COMPRESSION: PASS")


if __name__ == "__main__":
    main()
