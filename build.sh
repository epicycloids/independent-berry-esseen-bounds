#!/bin/sh
# Build each manuscript in a temporary directory.
set -eu

usage() {
    cat <<'USAGE'
Usage: ./build.sh [OUTPUT_DIRECTORY]
       ./build.sh --help

Build every manuscript with Tectonic. OUTPUT_DIRECTORY defaults to dist/;
relative paths are resolved from the repository root. The upper-bound PDF
is written directly to that directory, and each companion PDF to its own
subdirectory. Existing PDFs at those paths are replaced.

  -h, --help    Show this help and exit without building.
USAGE
}

if [ "$#" -gt 1 ]; then
    usage >&2
    exit 2
fi
case "${1-}" in
    -h|--help) usage; exit 0 ;;
    -*) usage >&2; exit 2 ;;
esac

cd "$(dirname "$0")"
DESTINATION="${1-dist}"
BUILD_DIR=$(mktemp -d)
trap 'rm -rf "$BUILD_DIR"' EXIT HUP INT TERM
mkdir -p "$DESTINATION"
for SOURCE_DIR in paper paper/*/; do
    test -f "$SOURCE_DIR/main.tex" || continue
    if [ "$SOURCE_DIR" = paper ]; then
        NAME=independent-berry-esseen
        OUTPUT_DIR="$DESTINATION"
    else
        NAME=$(basename "$SOURCE_DIR")
        OUTPUT_DIR="$DESTINATION/$NAME"
    fi
    WORK_DIR="$BUILD_DIR/$NAME"
    mkdir -p "$WORK_DIR" "$OUTPUT_DIR"
    cp "$SOURCE_DIR/"*.tex "$SOURCE_DIR/refs.bib" "$WORK_DIR/"
    tectonic --reruns 3 --keep-logs --outdir "$WORK_DIR" "$WORK_DIR/main.tex"
    if grep -Eq 'undefined (references|citations)|Label\(s\) may have changed' "$WORK_DIR/main.log"; then
        cat "$WORK_DIR/main.log" >&2
        echo "Unresolved references in $NAME" >&2
        exit 1
    fi
    cp "$WORK_DIR/main.pdf" "$OUTPUT_DIR/$NAME.pdf"
    printf 'Built %s/%s.pdf\n' "$OUTPUT_DIR" "$NAME"
done
