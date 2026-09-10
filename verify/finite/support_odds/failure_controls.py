#!/usr/bin/env python3
"""Check rejection of invalid sign patterns in the two odds programs."""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path


HERE = Path(__file__).resolve().parent


def run_mutation(
    path: Path,
    case: str,
    needle: str,
    injected_line: str,
    final_pass: str,
) -> None:
    source = path.read_text()
    if source.count(needle) != 1:
        raise AssertionError(f"{path.name}:{case}: mutation anchor changed")
    source = source.replace(needle, injected_line + needle, 1)

    stream = StringIO()
    caught = False
    try:
        with redirect_stdout(stream):
            exec(
                compile(source, str(path), "exec"),
                {"__name__": "__main__", "__file__": str(path)},
            )
    except AssertionError:
        caught = True

    output = stream.getvalue()
    if not caught:
        raise AssertionError(f"{path.name}:{case}: mutation was accepted")
    if final_pass in output:
        raise AssertionError(f"{path.name}:{case}: final PASS was printed")
    print(f"{path.name}:{case}: REJECTED")


def exact_controls() -> None:
    path = HERE / "exact_certificate.py"
    final_pass = "LIVE-STRIP ALL-DIMENSIONAL RESOURCE ENVELOPE: PASS"
    run_mutation(
        path,
        "small-negative",
        "    small_negative = [entry for entry in small_signs if entry[2] < 0]\n",
        "    small_signs.append((0, 0, -1))\n",
        final_pass,
    )
    run_mutation(
        path,
        "extra-small-zero",
        "    assert small_zero == [(60, 20, 0)]\n",
        "    small_zero.append((0, 0, 0))\n",
        final_pass,
    )
    run_mutation(
        path,
        "large-zero",
        "    if large_nonpositive:\n",
        "    large_nonpositive.append((0, 0, 0))\n",
        final_pass,
    )


def independent_controls() -> None:
    path = HERE / "independent_certificate.py"
    final_pass = "ALL-DIMENSIONAL PAIR COMPRESSION: PASS"
    run_mutation(
        path,
        "small-negative",
        "    small_negative = [entry for entry in small_signs if entry[2] < 0]\n",
        "    small_signs.append((0, 0, -1))\n",
        final_pass,
    )
    run_mutation(
        path,
        "extra-small-zero",
        "    assert small_zero == [(60, 20)], small_zero\n",
        "    small_zero.append((0, 0))\n",
        final_pass,
    )
    run_mutation(
        path,
        "large-zero",
        "    assert all(sign > 0 for sign in large_signs)\n",
        "    large_signs[0] = 0\n",
        final_pass,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "target",
        nargs="?",
        choices=("exact", "independent", "all"),
        default="all",
    )
    args = parser.parse_args()

    if args.target in {"exact", "all"}:
        exact_controls()
    if args.target in {"independent", "all"}:
        independent_controls()
    print("SUPPORT-ODDS FAILURE CONTROLS: PASS")


if __name__ == "__main__":
    main()
