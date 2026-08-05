#!/usr/bin/env bash
# Build the CASIM 1D column model driver, linking against build/libcasim.a
# (see ../build_casim_lib.sh). Produces column_model/casim_column_model.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HERE="$ROOT/column_model"
BUILD="$ROOT/build"
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -J "$BUILD" -I "$BUILD" -ffree-line-length-none}

if [[ ! -f "$BUILD/libcasim.a" ]]; then
  echo "libcasim.a not found - building CASIM library first..."
  "$ROOT/build_casim_lib.sh"
fi

echo "Compiling casim_column_model.F90"
$FC $FFLAGS -o "$HERE/casim_column_model" "$HERE/casim_column_model.F90" "$BUILD/libcasim.a"

echo "Built $HERE/casim_column_model"
