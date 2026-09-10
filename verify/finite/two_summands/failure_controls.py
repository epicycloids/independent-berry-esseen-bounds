#!/usr/bin/env python3
"""Check rejection of invalid lower bounds in the two-summand programs."""

from __future__ import annotations

import sys
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path


HERE = Path(__file__).resolve().parent


def run_mutation(
    path: Path,
    case: str,
    needle: str,
    injected_line: str,
    argv: list[str],
    final_pass: str,
) -> None:
    source = path.read_text()
    if source.count(needle) != 1:
        raise AssertionError(f"{path.name}:{case}: mutation anchor changed")
    source = source.replace(needle, injected_line + needle, 1)

    stream = StringIO()
    caught = False
    old_argv = sys.argv
    try:
        sys.argv = [str(path), *argv]
        try:
            with redirect_stdout(stream):
                exec(
                    compile(source, str(path), "exec"),
                    {"__name__": "__main__", "__file__": str(path)},
                )
        except (AssertionError, SystemExit):
            caught = True
    finally:
        sys.argv = old_argv

    output = stream.getvalue()
    if not caught:
        raise AssertionError(f"{path.name}:{case}: mutation was accepted")
    if final_pass in output:
        raise AssertionError(f"{path.name}:{case}: final PASS was printed")
    print(f"{path.name}:{case}: REJECTED")


def row01_controls() -> None:
    path = HERE / "certificate.py"
    argv = ["--owner", "01", "--max-depth", "60", "--max-boxes", "100000"]
    final_pass = "GLOBAL CLOSED-01 GAP-RATIO CERTIFICATE: PASS"
    run_mutation(
        path,
        "negative-live-gap",
        "        if all(gaps[owner] > 0 for owner in owners):\n",
        "        gaps = (gaps[0], -arb(1))\n",
        argv,
        final_pass,
    )
    run_mutation(
        path,
        "wrong-owned-mass",
        "        cdf01 = (1 - p1 * p2).upper()\n",
        "        p2 = q2\n",
        argv,
        final_pass,
    )


def row10_controls() -> None:
    path = HERE / "row10_certificate.py"
    argv = ["--max-depth", "42", "--max-boxes", "1000000"]
    final_pass = "GLOBAL CLOSED-10 VARIANCE-SHARE CERTIFICATE: PASS"
    run_mutation(
        path,
        "negative-live-gap",
        "        if gap > 0:\n",
        "        gap = -arb(1)\n",
        argv,
        final_pass,
    )
    run_mutation(
        path,
        "lowered-constant",
        "    penalty = ce.lower() * resource\n",
        "    ce -= arb(1) / 20\n",
        argv,
        final_pass,
    )


def main() -> None:
    row01_controls()
    row10_controls()
    print("TWO-BIT FAILURE CONTROLS: PASS")


if __name__ == "__main__":
    main()
