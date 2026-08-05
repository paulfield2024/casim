#!/usr/bin/env bash
# Build all CASIM unit test programs, linking each against build/libcasim.a
# (see ../build_casim_lib.sh). Each unit test is a standalone Fortran
# `program test_*` source file living in this directory; each is compiled
# into its own executable under unit_tests/bin/.
#
# Usage:
#   ./build.sh            # build all test_*.F90 files
#   ./build.sh test_foo    # build only unit_tests/test_foo.F90
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HERE="$ROOT/unit_tests"
BUILD="$ROOT/build"
BIN="$HERE/bin"
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -J "$BUILD" -I "$BUILD" -ffree-line-length-none}

if [[ ! -f "$BUILD/libcasim.a" ]]; then
  echo "libcasim.a not found - building CASIM library first..."
  "$ROOT/build_casim_lib.sh"
fi

mkdir -p "$BIN"

# Shared assertion-helper module used by every test_*.F90 (check_close,
# check_true, test_summary). Compile it once into $BUILD so its .mod file
# is available, and re-use the resulting .o when linking each test.
echo "Compiling test_utils.F90"
$FC $FFLAGS -c -o "$BUILD/test_utils.o" "$HERE/test_utils.F90"

if [[ $# -gt 0 ]]; then
  names=("$@")
else
  names=()
  for f in "$HERE"/test_*.F90; do
    [[ -e "$f" ]] || continue
    # test_utils.F90 is the shared helper module, not a standalone test
    # program - it is already compiled above, skip it here.
    [[ "$(basename "$f")" == "test_utils.F90" ]] && continue
    names+=("$(basename "${f%.F90}")")
  done
fi

if [[ ${#names[@]} -eq 0 ]]; then
  echo "No test_*.F90 files found in $HERE"
  exit 0
fi

for name in "${names[@]}"; do
  src="$HERE/$name.F90"
  if [[ ! -f "$src" ]]; then
    echo "error: $src not found" >&2
    exit 1
  fi
  echo "Compiling $name.F90"
  $FC $FFLAGS -o "$BIN/$name" "$src" "$BUILD/test_utils.o" "$BUILD/libcasim.a"
done

echo "Built $(printf '%s ' "${names[@]}")into $BIN"
