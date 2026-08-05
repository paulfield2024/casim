#!/usr/bin/env bash
# Build the CASIM microphysics library from the UMDP3-cleaned sources
# (src_clean/ + external_stubs/) into build_clean/libcasim.a, to verify
# that the style cleanup in src_clean/ has not changed any scientific
# behaviour (compare its unit-test/column-model output against the
# baseline built from src/).
#
# Usage: ./build_casim_lib_clean.sh [clean]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="$ROOT/build_clean"
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -J "$BUILD" -I "$BUILD" -ffree-line-length-none -fPIC -DDEF_MODEL=MODEL_MONC -DMODEL_MONC=4}

if [[ "${1:-}" == "clean" ]]; then
  rm -rf "$BUILD"
  echo "Cleaned $BUILD"
  exit 0
fi

mkdir -p "$BUILD"

echo "Determining build order..."
mapfile -t FILES < <(python3 "$ROOT/build_order.py" "$ROOT/external_stubs" "$ROOT/src_clean")

OBJS=()
for f in "${FILES[@]}"; do
  base=$(basename "$f" .F90)
  obj="$BUILD/${base}.o"
  echo "Compiling $f"
  $FC -c $FFLAGS -o "$obj" "$f"
  OBJS+=("$obj")
done

echo "Archiving build_clean/libcasim.a"
ar rcs "$BUILD/libcasim.a" "${OBJS[@]}"

echo "Done. Modules and library are in $BUILD"
