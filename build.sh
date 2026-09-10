#!/bin/sh
# Build the manuscript in a temporary directory.
set -eu
cd "$(dirname "$0")"
DESTINATION="${1-dist}"
BUILD_DIR=$(mktemp -d)
trap 'rm -rf "$BUILD_DIR"' EXIT HUP INT TERM
mkdir -p "$DESTINATION"
cp paper/main.tex paper/values.tex paper/refs.bib "$BUILD_DIR/"
# Rerun LaTeX for the bibliography and cross-references.
tectonic --reruns 3 --keep-logs --outdir "$BUILD_DIR" "$BUILD_DIR/main.tex"
if grep -Eq 'undefined (references|citations)|Label\(s\) may have changed' "$BUILD_DIR/main.log"; then
    cat "$BUILD_DIR/main.log" >&2
    echo 'Unresolved manuscript references' >&2
    exit 1
fi
cp "$BUILD_DIR/main.pdf" "$DESTINATION/independent-berry-esseen.pdf"
printf 'Built %s/independent-berry-esseen.pdf\n' "$DESTINATION"
