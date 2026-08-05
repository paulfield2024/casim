#!/usr/bin/env bash
# Build the CASIM microphysics library (src/ + external_stubs/) into a
# static library build/libcasim.a, compiling files in dependency order as
# computed by build_order.py.
#
# Usage: ./build_casim_lib.sh [clean]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="$ROOT/build"
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -J "$BUILD" -I "$BUILD" -ffree-line-length-none -fPIC -DDEF_MODEL=MODEL_MONC -DMODEL_MONC=4}

if [[ "${1:-}" == "clean" ]]; then
  rm -rf "$BUILD"
  echo "Cleaned $BUILD"
  exit 0
fi

mkdir -p "$BUILD"

echo "Determining build order..."
mapfile -t FILES < <(python3 "$ROOT/build_order.py" "$ROOT/external_stubs" "$ROOT/src")

OBJS=()
for f in "${FILES[@]}"; do
  base=$(basename "$f" .F90)
  obj="$BUILD/${base}.o"
  echo "Compiling $f"
  $FC -c $FFLAGS -o "$obj" "$f"
  OBJS+=("$obj")
done

echo "Archiving build/libcasim.a"
ar rcs "$BUILD/libcasim.a" "${OBJS[@]}"

echo "Done. Modules and library are in $BUILD"
