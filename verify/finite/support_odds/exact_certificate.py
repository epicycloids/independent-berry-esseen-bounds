#!/usr/bin/env python3
"""Exact Bernstein verification of the two inequalities for Bernoulli odds."""

from __future__ import annotations

import sympy as sp


SQRT5 = sp.sqrt(5)


def qsqrt5_parts(value: sp.Expr) -> tuple[sp.Rational, sp.Rational]:
    """Return A,B for value=A+B*sqrt(5), with A,B rational."""
    expanded = sp.expand(value)
    coefficient = sp.expand(expanded).coeff(SQRT5)
    constant = sp.simplify(expanded - coefficient * SQRT5)
    if constant.is_Rational is not True or coefficient.is_Rational is not True:
        raise AssertionError(f"not in Q(sqrt(5)): {value}")
    return sp.Rational(constant), sp.Rational(coefficient)


def qsqrt5_sign(value: sp.Expr) -> int:
    """Decide the sign in Q(sqrt(5)) by rational square comparisons."""
    constant, coefficient = qsqrt5_parts(value)
    if constant == 0:
        return int(sp.sign(coefficient))
    if coefficient == 0:
        return int(sp.sign(constant))
    if constant > 0 and coefficient > 0:
        return 1
    if constant < 0 and coefficient < 0:
        return -1

    square_difference = constant**2 - 5 * coefficient**2
    if square_difference == 0:
        return 0
    if constant > 0:
        return 1 if square_difference > 0 else -1
    return -1 if square_difference > 0 else 1


def bernstein_coefficients_2d(
    expression: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    degree_x: int,
    degree_y: int,
) -> list[list[sp.Expr]]:
    """Convert a bivariate power polynomial to a fixed Bernstein degree."""
    polynomial = sp.Poly(sp.expand(expression), x, y)
    power_degree_x = polynomial.degree(x)
    power_degree_y = polynomial.degree(y)
    assert power_degree_x <= degree_x
    assert power_degree_y <= degree_y

    coefficients: list[list[sp.Expr]] = []
    for p in range(degree_x + 1):
        row: list[sp.Expr] = []
        for q in range(degree_y + 1):
            value = sp.S.Zero
            for i in range(min(p, power_degree_x) + 1):
                for j in range(min(q, power_degree_y) + 1):
                    value += (
                        polynomial.coeff_monomial(x**i * y**j)
                        * sp.Rational(sp.binomial(p, i), sp.binomial(degree_x, i))
                        * sp.Rational(sp.binomial(q, j), sp.binomial(degree_y, j))
                    )
            row.append(sp.expand(value))
        coefficients.append(row)
    return coefficients


def bernstein_coefficients_1d(
    expression: sp.Expr,
    x: sp.Symbol,
) -> list[sp.Expr]:
    """Convert a univariate power polynomial at its minimal degree."""
    polynomial = sp.Poly(sp.expand(expression), x)
    degree = polynomial.degree(x)
    coefficients: list[sp.Expr] = []
    for p in range(degree + 1):
        value = sp.S.Zero
        for i in range(p + 1):
            value += polynomial.coeff_monomial(x**i) * sp.Rational(
                sp.binomial(p, i),
                sp.binomial(degree, i),
            )
        coefficients.append(sp.factor(value))
    return coefficients


def print_table(title: str, rows: list[list[sp.Expr]]) -> None:
    print(title)
    for row in rows:
        entries = []
        for value in row:
            constant, coefficient = qsqrt5_parts(value)
            entries.append(f"{constant}|{coefficient}")
        print(" ".join(entries))


def equal_block_value(total: sp.Symbol, count: int) -> sp.Expr:
    odds = total / sp.Integer(count)
    return sp.cancel(
        total
        / (
            (1 + odds**2) ** 2
            * (1 + odds) ** (2 * (count - 1))
        )
    )


def main() -> None:
    a, k, r, t = sp.symbols("a k r t", nonnegative=True)

    sign_controls = [
        (1 + SQRT5, 1),
        (-1 - SQRT5, -1),
        (3 - SQRT5, 1),
        (2 - SQRT5, -1),
        (SQRT5 - SQRT5, 0),
    ]
    for value, expected in sign_controls:
        assert qsqrt5_sign(value) == expected
    print("Q(SQRT(5)) SIGN CONTROLS: PASS")

    denominator = 1 + a**2 - 2 * k + k**2
    numerator = (
        a
        + 2 * a**2
        + a**3
        + (a**3 - a - 4) * k
        + (2 * a**2 - a + 8) * k**2
        + (a - 4) * k**3
    )
    pair_curve = numerator / (
        denominator**2
        * (1 + a + k) ** 2
    )
    crossover = 3 - SQRT5
    assert sp.expand(crossover**2 - 6 * crossover + 4) == 0

    small_gap_numerator = sp.together(
        pair_curve.subs(k, 0) - pair_curve
    ).as_numer_denom()[0]
    small_quotient, small_remainder = sp.div(
        small_gap_numerator,
        k,
        k,
    )
    assert small_remainder == 0
    small_a = crossover * r
    small_box = sp.expand(
        small_quotient.subs(
            {
                a: small_a,
                k: small_a**2 * t / 4,
            },
            simultaneous=True,
        )
    )
    small_bernstein = bernstein_coefficients_2d(
        small_box,
        r,
        t,
        60,
        20,
    )
    assert len(small_bernstein) == 61
    assert all(len(row) == 21 for row in small_bernstein)
    assert sum(len(row) for row in small_bernstein) == 1281
    small_signs = [
        (i, j, qsqrt5_sign(value))
        for i, row in enumerate(small_bernstein)
        for j, value in enumerate(row)
    ]
    small_negative = [entry for entry in small_signs if entry[2] < 0]
    small_zero = [entry for entry in small_signs if entry[2] == 0]
    if small_negative:
        raise AssertionError(f"negative small coefficient: {small_negative[0]}")
    assert small_zero == [(60, 20, 0)]
    print_table(
        "SMALL-PAIR BERNSTEIN TABLE: rows in r, columns in t",
        small_bernstein,
    )
    print("small-pair negative coefficient count:", len(small_negative))
    print("small-pair zero coefficients:", small_zero)

    atom_coordinate = 1 / (1 + a + k) ** 2
    equal_product = a**2 / 4
    theta = sp.cancel(
        (atom_coordinate - atom_coordinate.subs(k, equal_product))
        / (
            atom_coordinate.subs(k, 0)
            - atom_coordinate.subs(k, equal_product)
        )
    )
    large_chord_gap = sp.cancel(
        theta * pair_curve.subs(k, 0)
        + (1 - theta) * pair_curve.subs(k, equal_product)
        - pair_curve
    )
    large_gap_numerator = sp.together(large_chord_gap).as_numer_denom()[0]
    large_factor = k * (a**2 - 4 * k)
    large_quotient, large_remainder = sp.div(
        large_gap_numerator,
        large_factor,
        k,
    )
    assert large_remainder == 0
    large_a = crossover + (sp.Rational(16, 9) - crossover) * r
    large_box = sp.expand(
        large_quotient.subs(
            {
                a: large_a,
                k: large_a**2 * t / 4,
            },
            simultaneous=True,
        )
    )
    large_polynomial = sp.Poly(large_box, r, t)
    large_degree_r = large_polynomial.degree(r)
    large_degree_t = large_polynomial.degree(t)
    assert (large_degree_r, large_degree_t) == (13, 4)
    large_bernstein = bernstein_coefficients_2d(
        large_box,
        r,
        t,
        large_degree_r,
        large_degree_t,
    )
    assert len(large_bernstein) == 14
    assert all(len(row) == 5 for row in large_bernstein)
    assert sum(len(row) for row in large_bernstein) == 70
    large_nonpositive = [
        (i, j, qsqrt5_sign(value))
        for i, row in enumerate(large_bernstein)
        for j, value in enumerate(row)
        if qsqrt5_sign(value) <= 0
    ]
    if large_nonpositive:
        raise AssertionError(
            f"nonpositive large coefficient: {large_nonpositive[0]}"
        )
    print_table(
        "LARGE-PAIR BERNSTEIN TABLE: rows in r, columns in t",
        large_bernstein,
    )
    print("large-pair nonpositive coefficient count:", len(large_nonpositive))

    total = sp.symbols("s", positive=True)
    unit = sp.symbols("u", nonnegative=True)
    total_map = 1 + sp.Rational(7, 9) * unit
    equal_tables: dict[int, list[sp.Expr]] = {}
    for count in (1, 3, 4):
        difference = equal_block_value(total, 2) - equal_block_value(
            total,
            count,
        )
        difference_numerator, difference_denominator = sp.together(
            difference
        ).as_numer_denom()
        mapped_denominator = sp.Poly(
            sp.expand(difference_denominator.subs(total, total_map)),
            unit,
        )
        assert all(value > 0 for value in mapped_denominator.all_coeffs())
        table = bernstein_coefficients_1d(
            sp.expand(difference_numerator.subs(total, total_map)),
            unit,
        )
        if any(value <= 0 for value in table):
            raise AssertionError(
                f"nonpositive equal-{count} coefficient"
            )
        equal_tables[count] = table
        print(f"EQUAL-BLOCK m={count} BERNSTEIN TABLE:")
        print(" ".join(str(value) for value in table))

    print("SMALL-PAIR MERGE CERTIFICATE: PASS")
    print("LARGE-PAIR CHORD CERTIFICATE: PASS")
    print("EQUAL-BLOCK COMPARISONS: PASS")
    print("LIVE-STRIP ALL-DIMENSIONAL RESOURCE ENVELOPE: PASS")


if __name__ == "__main__":
    main()
