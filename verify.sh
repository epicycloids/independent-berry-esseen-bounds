#!/bin/sh
# Check supplied records and optionally build the paper.
set -eu
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
BUILD=yes
case "${1-}" in
  "") ;;
  --no-latex) BUILD=no ;;
  *) echo "usage: ./verify.sh [--no-latex]" >&2; exit 2 ;;
esac
test "$#" -le 1
sha256sum --check --quiet SHA256SUMS
uv run --frozen python -m unittest discover -s verify -p 'test_*.py'
uv run --frozen python verify/check_saved.py
uv run --frozen python verify/check_replays.py
uv run --frozen python verify/check_scalar_certificates.py
uv run --frozen python verify/paper_values.py
if [ "$BUILD" = yes ]; then
  ./build.sh
fi
