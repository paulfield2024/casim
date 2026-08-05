#!/usr/bin/env bash
# Run every built unit test executable under unit_tests/bin/ and save its
# stdout as the Known-Good-Output (KGO) baseline under kgo/<test_name>.txt.
#
# Usage:
#   ./run_all.sh          # run all built tests, write/refresh kgo/*.txt
#
# Each test program is expected to:
#   - print human-readable PASS/FAIL lines for its own internal checks
#   - print numeric results with enough precision to serve as a KGO baseline
#   - exit with a non-zero status if any of its internal checks fail
#
# This script itself exits non-zero if any test executable exits non-zero,
# so it can be used as a single pass/fail gate (e.g. in CI).
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/unit_tests/bin"
KGO="$ROOT/kgo"

mkdir -p "$KGO"

if [[ ! -d "$BIN" ]] || [[ -z "$(ls -A "$BIN" 2>/dev/null)" ]]; then
  echo "No built test executables found in $BIN - run unit_tests/build.sh first."
  exit 1
fi

overall_status=0

for exe in "$BIN"/*; do
  [[ -x "$exe" ]] || continue
  name="$(basename "$exe")"
  out="$KGO/$name.txt"
  echo "Running $name"
  if "$exe" > "$out" 2>&1; then
    echo "  PASS -> $out"
  else
    echo "  FAIL (exit $?) -> $out"
    overall_status=1
  fi
done

exit $overall_status
