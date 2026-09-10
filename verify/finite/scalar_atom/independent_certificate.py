#!/usr/bin/env python3
"""Separate evaluation of the scalar atom inequality.

The evaluator reconstructs the rational partition and the majorant over
all positive orders. It computes even moments by a direct coefficient sum.
It uses the same Arb library as compact_certificate.py, without importing
that program.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import platform
import sys
from fractions import Fraction
from math import ceil, comb
from pathlib import Path

from flint import arb, ctx


SOURCE = Path(__file__).resolve().parent
SOURCE_FILES = ("compact_certificate.py",)
EXPECTED_SOURCE_HASHES = {
    "compact_certificate.py": "35410ea2f44560967980d836834cf856c4ec088b9c318d2d17e7c207ca53b75f",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_source_hashes() -> None:
    for name in SOURCE_FILES:
        actual = sha256(SOURCE / name)
        expected = EXPECTED_SOURCE_HASHES[name]
        if actual != expected:
            raise AssertionError(
                f"source hash mismatch for {name}: {actual} != {expected}"
            )
        print(f"source-sha256 {actual} {name}")


def constants() -> tuple[arb, arb, arb, arb]:
    root10 = arb(10).sqrt()
    ke = (root10 + 3) / 6
    star = root10 - 3
    lower_skew = 1 / root10
    upper_skew = star + (1 / (5 * ke)).sqrt()
    return ke, star, lower_skew, upper_skew


def normalized_parameter_ball(index: int, cells: int) -> arb:
    """Return an outward ball for the exact interval [index/cells,(index+1)/cells]."""

    midpoint = Fraction(2 * index + 1, 2 * cells)
    radius = Fraction(1, 2 * cells)
    return arb(
        f"{midpoint.numerator}/{midpoint.denominator}",
        f"{radius.numerator}/{radius.denominator}",
    )


def exact_partition(cells: int) -> list[tuple[Fraction, Fraction]]:
    return [(Fraction(j, cells), Fraction(j + 1, cells)) for j in range(cells)]


def require_exact_cover(intervals: list[tuple[Fraction, Fraction]]) -> None:
    if not intervals:
        raise AssertionError("empty cover")
    ordered = sorted(intervals)
    if len(set(ordered)) != len(ordered):
        raise AssertionError("duplicate box")
    if ordered[0][0] != 0:
        raise AssertionError(f"left boundary missing: {ordered[0][0]}")
    if ordered[-1][1] != 1:
        raise AssertionError(f"right boundary missing: {ordered[-1][1]}")
    for previous, current in zip(ordered, ordered[1:]):
        if previous[1] != current[0]:
            raise AssertionError(f"gap or overlap: {previous} then {current}")


def coverage_and_deletion_controls(cells: int) -> None:
    cover = exact_partition(cells)
    require_exact_cover(cover)

    for label, damaged in (
        ("deleted-interior", cover[: cells // 2] + cover[cells // 2 + 1 :]),
        ("deleted-left-boundary", cover[1:]),
        ("deleted-right-boundary", cover[:-1]),
    ):
        try:
            require_exact_cover(damaged)
        except AssertionError:
            print(f"positive-control {label}: FIRED")
        else:
            raise AssertionError(f"coverage checker missed {label}")

    slabs = [(2 * m - 2, 2 * m) for m in range(2, 65)]
    assert len(slabs) == 63
    assert slabs[0][0] == 2 and slabs[-1][1] == 128
    assert all(left[1] == right[0] for left, right in zip(slabs, slabs[1:]))
    print(
        f"exact-cover cells={len(cover)} slabs={len(slabs)} box-slabs={len(cover)*len(slabs)}"
    )
    print("exact-cover parameter-endpoints=0,1 order-endpoints=2,128: PASS")


def even_moment_direct(m: int, t: arb) -> arb:
    """Evaluate q_m by a direct positive coefficient sum."""

    x = t * t
    denominator = 4**m
    answer = arb(0)
    power = arb(1)
    for k in range(m + 1):
        numerator = comb(2 * k, k) * comb(2 * (m - k), m - k)
        answer += (arb(numerator) / denominator) * power
        power *= x
    return answer


def even_moment_exact_fraction(m: int, t_squared: Fraction) -> Fraction:
    answer = Fraction(0)
    for k in range(m + 1):
        coefficient = Fraction(comb(2 * k, k) * comb(2 * (m - k), m - k), 4**m)
        answer += coefficient * t_squared**k
    return answer


def global_order_peak(m: int, t: arb) -> arb:
    """Independent implementation of the all-positive-order moment peak."""

    q_lower = even_moment_direct(m - 1, t)
    q_upper = even_moment_direct(m, t)
    ratio = q_upper / q_lower
    ell = (q_lower / q_upper).log()
    if not ell > 0:
        raise AssertionError(
            f"could not certify q_{m}<q_{m-1}: ratio={ratio}, ell={ell}"
        )

    # The maximum of sqrt(nu) exp(-ell*nu/2) over nu>0 occurs at nu=1/ell.
    # This algebraic form is inclusion-equivalent but was derived independently.
    c = 1 - t * t
    lower_order = 2 * m - 2
    log_peak = (
        ((arb.pi() * c) / (2 * arb(1).exp() * ell)).log() / 2
        + q_lower.log()
        + arb(lower_order) * ell / 2
    )
    return log_peak.exp()


def target_rhs(t: arb, baseline: arb | None = None) -> arb:
    ke, star, _, _ = constants()
    if baseline is None:
        baseline = arb(1)
    delta = t - star
    return baseline + ke * delta * delta


def direct_g(nu: arb, t: arb) -> arb:
    c = 1 - t * t
    return (arb.pi() * nu * c / 2).sqrt() * c.hypgeom_2f1(-nu / 2, arb(1) / 2, 1)


def mathematical_controls() -> None:
    ke, star, lower_skew, upper_skew = constants()

    # Exact coefficient anchors independent of special functions.
    for m in range(1, 13):
        assert even_moment_exact_fraction(m, Fraction(1)) == 1
        assert even_moment_exact_fraction(m, Fraction(0)) == Fraction(
            comb(2 * m, m), 4**m
        )

    # Direct polynomial versus Arb hypergeometric evaluation at unrelated anchors.
    for m, t_text in ((1, "0.37"), (7, "0.47"), (19, "0.52"), (64, "0.4")):
        t = arb(t_text)
        c = 1 - t * t
        polynomial = even_moment_direct(m, t)
        hypergeometric = c.hypgeom_2f1(-m, arb(1) / 2, 1)
        if not (polynomial - hypergeometric).contains(0):
            raise AssertionError(f"moment normalization mismatch at m={m}, t={t}")

    # Adjacent-moment inequality at noninteger orders, including points on both halves of slabs.
    for nu_text, t_text in (
        ("3.25", "0.33"),
        ("5.3", "0.47"),
        ("17.75", "0.58"),
        ("127.2", "0.4"),
    ):
        nu = arb(nu_text)
        t = arb(t_text)
        m = max(2, ceil(float(nu_text) / 2))
        lam = (nu - (2 * m - 2)) / 2
        qlo = even_moment_direct(m - 1, t)
        qhi = even_moment_direct(m, t)
        moment_envelope = qlo ** (1 - lam) * qhi**lam
        c = 1 - t * t
        g_envelope = (arb.pi() * nu * c / 2).sqrt() * moment_envelope
        actual = direct_g(nu, t)
        if not g_envelope >= actual:
            raise AssertionError(f"adjacent moment control failed at nu={nu}, t={t}")

    # A safe negative control and a known G>1 normalization control.
    known_over_one = direct_g(arb(5), arb("0.51"))
    assert known_over_one > 1
    lowered_gap = target_rhs(star, arb(99) / 100) - direct_g(arb(1000), star)
    assert lowered_gap < 0

    # The reversed comparison must not pass at a known safe compact point.
    sample_t = lower_skew
    sample_peak = global_order_peak(64, sample_t)
    sample_gap = target_rhs(sample_t) - sample_peak
    assert sample_gap > 0
    reversed_gap = sample_peak - target_rhs(sample_t)
    assert reversed_gap < 0

    # Interfaces and tangent constants.
    assert (1 - lower_skew * lower_skew - arb(9) / 10).contains(0)
    assert (target_rhs(upper_skew) - arb(6) / 5).contains(0)
    print("mathematical-controls: PASS")
    print("  known-G-over-one", known_over_one)
    print("  lowered-baseline-positive-control", lowered_gap)
    print("  reversed-comparison-positive-control", reversed_gap)
    print("  sample-safe-gap", sample_gap)
    print("  K_E", ke)


def bessel_tail_control() -> None:
    ke, star, lower_skew, upper_skew = constants()
    c_min = 1 - upper_skew * upper_skew
    w_min = 32 * c_min

    # Use unscaled I_0 times exp(-W), unlike the submitted scaled=True call.
    unscaled = (2 * arb.pi() * w_min).sqrt() * w_min.bessel_i(0) * (-w_min).exp()
    scaled_crosscheck = (2 * arb.pi() * w_min).sqrt() * w_min.bessel_i(0, scaled=True)
    if not (unscaled - scaled_crosscheck).contains(0):
        raise AssertionError("scaled and unscaled Bessel backends disagree")
    rhs_min = 1 + ke * (lower_skew - star) ** 2
    gap = rhs_min - unscaled
    assert w_min > 1 and gap > 0
    bad_gap = rhs_min - (unscaled + arb(1) / 20)
    assert bad_gap < 0
    print("bessel-tail-control: PASS")
    print("  W-min", w_min)
    print("  H-unscaled", unscaled)
    print("  H-scaled-difference", unscaled - scaled_crosscheck)
    print("  tail-gap", gap)
    print("  inflated-majorant-positive-control", bad_gap)


def error_function_controls() -> None:
    _, star, lower_skew, _ = constants()
    c_boundary = arb(9) / 10
    endpoint_defect = -(1 - c_boundary).log() - c_boundary * arb.pi() ** 2 / 4
    assert endpoint_defect > 0

    c_star = 1 - star * star
    for nu in (arb(1), arb(2), arb(128), arb(1000)):
        argument = arb.pi() * (nu * c_star / 8).sqrt()
        erf_bound = argument.erf()
        actual = direct_g(nu, star)
        # At large arguments erf rounds to a ball containing one, whereas erfc remains a
        # directly positive ball and certifies the strict inequality analytically used.
        assert argument.erfc() > 0
        assert erf_bound >= actual
    assert (1 - lower_skew * lower_skew - c_boundary).contains(0)
    print("error-function-controls: PASS")
    print("  endpoint-defect-at-c=9/10", endpoint_defect)


def direct_even_order_g(m: int, t: arb) -> arb:
    c = 1 - t * t
    return (arb.pi() * m * c).sqrt() * even_moment_direct(m, t)


def small_subcover_control(cells: int) -> None:
    """Certify an unmodified nine-box subcover as a negative control."""

    _, _, lower_skew, upper_skew = constants()
    checked = 0
    for index in (0, cells // 2, cells - 1):
        parameter = normalized_parameter_ball(index, cells)
        t = lower_skew + (upper_skew - lower_skew) * parameter
        for m in (2, 32, 64):
            gap = target_rhs(t) - global_order_peak(m, t)
            if not gap > 0:
                raise AssertionError(
                    f"safe subcover failed at index={index}, m={m}: {gap}"
                )
            checked += 1
    assert checked == 9
    print("unmodified-small-subcover-negative-control: PASS boxes=9")


def certify_compact(cells: int) -> None:
    _, _, lower_skew, upper_skew = constants()
    expected = cells * 63
    accepted = 0
    common_floor = arb(19) / 1000
    weakest: tuple[arb, int, int, arb] | None = None

    left_edge_min: arb | None = None
    right_edge_min: arb | None = None
    nu2_edge_min: arb | None = None
    nu128_edge_min: arb | None = None

    for index in range(cells):
        parameter = normalized_parameter_ball(index, cells)
        t = lower_skew + (upper_skew - lower_skew) * parameter
        right = target_rhs(t)

        # Direct checks of both order boundaries over every skew cell.
        gap_nu2 = right - direct_even_order_g(1, t)
        gap_nu128 = right - direct_even_order_g(64, t)
        if not gap_nu2 > 0 or not gap_nu128 > 0:
            raise AssertionError(f"order-edge failure at cell {index}")
        if nu2_edge_min is None or gap_nu2.lower() < nu2_edge_min.lower():
            nu2_edge_min = gap_nu2
        if nu128_edge_min is None or gap_nu128.lower() < nu128_edge_min.lower():
            nu128_edge_min = gap_nu128

        for m in range(2, 65):
            peak = global_order_peak(m, t)
            gap = right - peak
            if not gap > common_floor:
                raise AssertionError(
                    f"compact failure index={index} m={m} t={t} gap={gap}"
                )
            accepted += 1
            if weakest is None or gap.lower() < weakest[0].lower():
                weakest = (gap, index, m, t)

    if accepted != expected:
        raise AssertionError(f"wrong accepted count: {accepted} != {expected}")
    assert weakest is not None

    # Degenerate balls at both skew boundaries check every slab on the exact interfaces.
    for m in range(2, 65):
        left_gap = target_rhs(lower_skew) - global_order_peak(m, lower_skew)
        right_gap = target_rhs(upper_skew) - global_order_peak(m, upper_skew)
        assert left_gap > common_floor and right_gap > common_floor
        if left_edge_min is None or left_gap.lower() < left_edge_min.lower():
            left_edge_min = left_gap
        if right_edge_min is None or right_gap.lower() < right_edge_min.lower():
            right_edge_min = right_gap

    gap, index, m, t = weakest
    print("independent-compact-certificate: PASS")
    print("  cells", cells)
    print("  slabs", 63)
    print("  expected-box-slabs", expected)
    print("  accepted-box-slabs", accepted)
    print("  missing-box-slabs", expected - accepted)
    print("  duplicate-box-slabs", 0)
    print("  weakest-index", index)
    print("  weakest-m", m)
    print("  weakest-t-box", t)
    print("  weakest-gap", gap)
    print("  weakest-gap-lower", gap.lower().str(40))
    print("  left-skew-edge-min", left_edge_min)
    print("  right-skew-edge-min", right_edge_min)
    print("  nu=2-edge-min", nu2_edge_min)
    print("  nu=128-edge-min", nu128_edge_min)


def main() -> None:
    precision = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    cells = int(sys.argv[2]) if len(sys.argv) > 2 else 2048
    if precision < 128 or cells < 2:
        raise SystemExit(
            "usage: independent_certificate.py [precision>=128] [cells>=2]"
        )
    ctx.prec = precision
    print("python", platform.python_version())
    print("python-flint", importlib.metadata.version("python-flint"))
    print("arb-precision-bits", precision)
    verify_source_hashes()
    coverage_and_deletion_controls(cells)
    mathematical_controls()
    error_function_controls()
    bessel_tail_control()
    small_subcover_control(cells)
    certify_compact(cells)


if __name__ == "__main__":
    main()
