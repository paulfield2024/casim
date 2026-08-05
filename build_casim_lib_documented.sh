#!/usr/bin/env bash
# Build the CASIM microphysics library from the paper-documented sources
# (src_clean_documented/ + external_stubs/) into
# build_clean_documented/libcasim.a, to verify that adding paper-citation
# documentation to src_clean_documented/ has not changed any scientific
# behaviour (compare its unit-test/column-model output against the
# existing src_clean/ baseline).
#
# Usage: ./build_casim_lib_documented.sh [clean]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="$ROOT/build_clean_documented"
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -J "$BUILD" -I "$BUILD" -ffree-line-length-none -fPIC -DDEF_MODEL=MODEL_MONC -DMODEL_MONC=4}

if [[ "${1:-}" == "clean" ]]; then
  rm -rf "$BUILD"
  echo "Cleaned $BUILD"
  exit 0
fi

mkdir -p "$BUILD"

echo "Determining build order..."
mapfile -t FILES < <(python3 "$ROOT/build_order.py" "$ROOT/external_stubs" "$ROOT/src_clean_documented")

OBJS=()
for f in "${FILES[@]}"; do
  base=$(basename "$f" .F90)
  obj="$BUILD/${base}.o"
  echo "Compiling $f"
  $FC -c $FFLAGS -o "$obj" "$f"
  OBJS+=("$obj")
done

echo "Archiving build_clean_documented/libcasim.a"
ar rcs "$BUILD/libcasim.a" "${OBJS[@]}"

echo "Done. Modules and library are in $BUILD"
