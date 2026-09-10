#!/usr/bin/env python3
"""Formula, boundary, and rejection checks for the lower interior inequality.

Higher-precision runs reuse the primary program.
"""

from __future__ import annotations

import hashlib
import io
from fractions import Fraction
from pathlib import Path
from contextlib import redirect_stdout
import sys
import types

import sympy as sp


SOURCE = Path(__file__).resolve().parent / "row10_certificate.py"
TERMINAL_PASS = "GLOBAL CLOSED-10 VARIANCE-SHARE CERTIFICATE: PASS"
SOURCE_SHA256 = "4875833df0eb05dea1fed76cc0a69513a06d466502dc7138883f8cd2192c3fcd"


def exact_reduction() -> None:
    p1, q1, p2, q2, alpha, beta = sp.symbols(
        "p1 q1 p2 q2 alpha beta",
        positive=True,
    )
    h1 = alpha / sp.sqrt(p1 * q1)
    h2 = beta / sp.sqrt(p2 * q2)
    resource = (
        p1 * q1 * (p1**2 + q1**2) * h1**3
        + p2 * q2 * (p2**2 + q2**2) * h2**3
    )
    shape1 = (p1**2 + q1**2) / sp.sqrt(p1 * q1)
    shape2 = (p2**2 + q2**2) / sp.sqrt(p2 * q2)
    assert sp.simplify(resource - alpha**3 * shape1 - beta**3 * shape2) == 0

    phase10 = q1 * h1 - p2 * h2
    expected_phase = (
        alpha * sp.sqrt(q1 / p1) - beta * sp.sqrt(p2 / q2)
    )
    assert sp.simplify(phase10 - expected_phase) == 0

    pp1, pp2 = sp.symbols("pp1 pp2")
    qq1 = 1 - pp1
    qq2 = 1 - pp2
    strict_row10_cdf = qq1 * qq2 + pp1 * qq2
    assert sp.expand(strict_row10_cdf - qq2) == 0
    collision_cdf = 1 - pp1 * pp2
    assert sp.expand(collision_cdf - qq2 - qq1 * pp2) == 0

    variance = sp.symbols("v", positive=True)
    shape_square = (1 - 2 * variance) ** 2 / variance
    assert sp.simplify(
        sp.diff(shape_square, variance) - (4 - 1 / variance**2)
    ) == 0
    print("EXACT VARIANCE-SHARE REDUCTION: PASS")
    print("EXACT STRICT-ORDER CDF MASS: PASS")
    print("EQUAL-GAP CDF CORRECTION: q1*p2")
    print("SHAPE MONOTONICITY REDUCTION: PASS")


def execute_source(source: str, argv: list[str]) -> tuple[bool, str]:
    stream = io.StringIO()
    failed = False
    old_argv = sys.argv
    try:
        sys.argv = [str(SOURCE), *argv]
        try:
            with redirect_stdout(stream):
                exec(
                    compile(source, str(SOURCE), "exec"),
                    {"__name__": "__main__", "__file__": str(SOURCE)},
                )
        except (AssertionError, SystemExit):
            failed = True
    finally:
        sys.argv = old_argv
    return failed, stream.getvalue()


def load_source() -> dict[str, object]:
    source = SOURCE.read_text()
    assert hashlib.sha256(source.encode()).hexdigest() == SOURCE_SHA256
    name = "row10_audit_target"
    module = types.ModuleType(name)
    module.__file__ = str(SOURCE)
    sys.modules[name] = module
    try:
        exec(compile(source, str(SOURCE), "exec"), module.__dict__)
    finally:
        del sys.modules[name]
    return module.__dict__


def boundary_and_split_checks() -> None:
    namespace = load_source()
    box_type = namespace["Box"]
    row10_lower = namespace["row10_lower"]
    shape_at = namespace["shape_at"]
    split_box = namespace["split_box"]

    try:
        shape_at(Fraction(0))
    except AssertionError:
        pass
    else:
        raise AssertionError("deterministic shape evaluation was not rejected")

    probes = (
        box_type((Fraction(0), Fraction(1)), (Fraction(0), Fraction(1)), (Fraction(0), Fraction(1))),
        box_type((Fraction(1, 2), Fraction(1)), (Fraction(1, 4), Fraction(1, 2)), (Fraction(1, 4), Fraction(1, 2))),
        box_type((Fraction(1, 4), Fraction(1, 2)), (Fraction(1, 4), Fraction(1, 2)), (Fraction(1, 4), Fraction(1, 2))),
    )
    for probe in probes:
        row10_lower(probe)

    for axis in range(3):
        bounds = [
            (Fraction(1, 7), Fraction(2, 7)),
            (Fraction(1, 9), Fraction(2, 9)),
            (Fraction(1, 11), Fraction(2, 11)),
        ]
        bounds[axis] = (Fraction(1, 20), Fraction(19, 20))
        parent = box_type(*bounds, depth=3)
        left, right = split_box(parent)
        midpoint = sum(bounds[axis], Fraction(0)) / 2
        left_bounds = (left.p1_bounds, left.p2_bounds, left.angle_bounds)
        right_bounds = (right.p1_bounds, right.p2_bounds, right.angle_bounds)
        assert left_bounds[axis] == (bounds[axis][0], midpoint)
        assert right_bounds[axis] == (midpoint, bounds[axis][1])
        assert left.depth == right.depth == 4
        for other in range(3):
            if other != axis:
                assert left_bounds[other] == right_bounds[other] == bounds[other]

    print("DETERMINISTIC/BOUNDARY BRANCH PROBES: PASS")
    print("EXACT RATIONAL THREE-AXIS SPLITTING: PASS")


def failure_gate_checks() -> None:
    source = SOURCE.read_text()
    cases = (
        ("invalid-budget", ["--max-depth", "0"]),
        ("box-budget", ["--max-boxes", "1"]),
        ("depth-budget", ["--max-depth", "1"]),
    )
    for case, argv in cases:
        failed, output = execute_source(source, argv)
        assert failed
        assert TERMINAL_PASS not in output
        print(f"FAILURE GATE {case}: REJECTED")


def precision_rerun() -> None:
    source = SOURCE.read_text()
    assert hashlib.sha256(source.encode()).hexdigest() == SOURCE_SHA256
    assert source.count("    ctx.prec = 192\n") == 1
    source = source.replace("    ctx.prec = 192\n", "    ctx.prec = 256\n")
    failed, output = execute_source(
        source,
        ["--max-depth", "42", "--max-boxes", "1000000"],
    )
    assert not failed
    assert TERMINAL_PASS in output
    digest = hashlib.sha256(output.encode()).hexdigest()
    print("256-BIT FULL RERUN: PASS")
    print("256-bit output SHA-256:", digest)
    for line in output.splitlines():
        if "processed-boxes" in line or "weakest-gap" in line:
            print(line.strip())


def mutation_control(case: str, needle: str, replacement: str) -> None:
    source = SOURCE.read_text()
    assert source.count(needle) == 1
    source = source.replace(needle, replacement)
    failed, output = execute_source(
        source,
        ["--max-depth", "42", "--max-boxes", "1000000"],
    )
    assert failed
    assert TERMINAL_PASS not in output
    print(f"INDEPENDENT MUTATION {case}: REJECTED")


def main() -> None:
    exact_reduction()
    boundary_and_split_checks()
    failure_gate_checks()
    precision_rerun()
    mutation_control(
        "full-cdf",
        "    cdf_upper = 1 - rational(p2_lo)\n",
        "    cdf_upper = arb(1)\n",
    )
    mutation_control(
        "zero-normal-payment",
        "        normal = normal_cdf(positive_phase - negative_phase).lower()\n",
        "        normal = arb(0)\n",
    )
    print("INDEPENDENT ROW10 CHECKS: PASS")


if __name__ == "__main__":
    main()
