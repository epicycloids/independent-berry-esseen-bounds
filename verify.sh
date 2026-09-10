#!/bin/sh
# Check supplied records and optionally build the papers.
set -eu

usage() {
    cat <<'USAGE'
Usage: ./verify.sh [--no-latex]
       ./verify.sh --help

Check file hashes, stored upper-bound records, and the default scalar and
finite certificate checks. The default also builds all PDFs under dist/.

  --no-latex    Run the checks without building PDFs.
  -h, --help    Show this help and exit without checking or building.

Requires uv and Python 3.12; PDF builds also require Tectonic.
See verify/README.md and certificates/finite/README.md for reevaluation
commands and the scope of the recorded checks.
USAGE
}

if [ "$#" -gt 1 ]; then
    usage >&2
    exit 2
fi
BUILD=yes
case "${1-}" in
    "") ;;
    --no-latex) BUILD=no ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
esac

cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
sha256sum --check --quiet SHA256SUMS
uv run --frozen python -m unittest discover -s verify -p 'test_*.py'
uv run --frozen python verify/check_saved.py
uv run --frozen python verify/check_replays.py
uv run --frozen python verify/check_scalar_certificates.py
uv run --frozen python verify/paper_values.py
uv run --frozen python verify/finite/check.py
if [ "$BUILD" = yes ]; then
    ./build.sh
fi
