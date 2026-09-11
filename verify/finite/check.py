"""Check the finite-results supplement; complete calculations are opt-in."""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys

from flint import ctx

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load(relative: str):
    path = HERE / relative
    name = "finite_" + relative.replace("/", "_").replace(".", "_")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def check_sources():
    record = json.loads((ROOT / "certificates/finite/results.json").read_text())
    for name, expected in record["source_sha256"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Finite certificate source differs: {name}")
    return len(record["source_sha256"])


def fixed_checks():
    """Exercise formulas, boundary conventions, and rejection of invalid inputs."""
    ctx.prec = 256
    stream = io.StringIO()
    with redirect_stdout(stream):
        load("common_diameter/singleton_certificate.py").fixed_controls()
        load("scalar_atom/compact_certificate.py").fixed_controls()
        scalar = load("scalar_atom/independent_certificate.py")
        scalar.verify_source_hashes()
        scalar.coverage_and_deletion_controls(2048)
        scalar.mathematical_controls()
        scalar.error_function_controls()
        scalar.bessel_tail_control()
        algebra = load("support_odds/exact_certificate.py")
        for value, sign in ((3-algebra.SQRT5, 1), (2-algebra.SQRT5, -1),
                            (algebra.SQRT5-algebra.SQRT5, 0)):
            if algebra.qsqrt5_sign(value) != sign:
                raise AssertionError("Incorrect algebraic sign")
    return stream.getvalue()


def complete_runs(independent: bool):
    commands = [
        ["common_diameter/singleton_certificate.py", "--max-depth", "38"],
        ["scalar_atom/compact_certificate.py"],
        ["support_odds/exact_certificate.py"],
    ]
    if independent:
        commands += [
            ["scalar_atom/independent_certificate.py", "256", "2048"],
            ["support_odds/independent_certificate.py"],
            ["support_odds/failure_controls.py"],
        ]
    for name, *arguments in commands:
        print(f"Running {name}", flush=True)
        subprocess.run([sys.executable, str(HERE / name), *arguments], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true",
                        help="Rerun all primary certificates, one at a time.")
    parser.add_argument("--independent", action="store_true",
                        help="With --full, also run the supplied independent checks.")
    parser.add_argument("--verbose", action="store_true",
                        help="With --full, print the individual fixed-check results.")
    args = parser.parse_args()
    if args.independent and not args.full:
        parser.error("--independent requires --full")
    count = check_sources()
    current = load("current/check.py")
    current.controls()
    current.verify()
    log = ""
    if args.verbose:
        print(log, end="")
    print(f"Finite supplement: {count} source hashes and exact current partitions passed; stored leaf enclosures were not recomputed.")
    if args.full:
        current.verify(arithmetic=True)
        log = fixed_checks()
        if args.verbose: print(log, end="")
        complete_runs(args.independent)


if __name__ == "__main__":
    main()
