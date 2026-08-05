# KGO verification results: src_clean_documented

This directory holds the verification artifacts produced when checking
that adding CASIM-paper documentation to `src_clean/` (producing
`src_clean_documented/`) made no scientific/numerical change.

## Method
1. `tools/verify_equivalence.py` was run for all 60 `.F90` files,
   comparing `src_clean/<f>` against `src_clean_documented/<f>` after
   stripping comments/whitespace - all 60 reported OK (token-equivalent).
2. `build_casim_lib_documented.sh` built `build_clean_documented/libcasim.a`
   from `src_clean_documented/` + `external_stubs/` - build succeeded
   with no errors.
3. All 9 `unit_tests/test_*.F90` programs were rebuilt against
   `build_clean_documented/libcasim.a` and run; stdout for each is saved
   here under `unit_tests/*.txt` and matches the existing `kgo/*.txt`
   baselines byte-for-byte (diff clean).
4. `column_model/casim_column_model.F90` was rebuilt against
   `build_clean_documented/libcasim.a`, run from an isolated scratch
   directory, and converted to netCDF via `column_model/make_netcdf.py`.
   All 113 variables in `column_model/casim_column_documented.nc` are
   numerically identical (`numpy.array_equal`) to the existing
   `column_model/output/casim_column_kgo.nc` baseline.

## Result
No scientific/numerical behaviour changed by the documentation pass:
all unit tests and all 113 column-model variables match exactly.

## Contents
- `unit_tests/*.txt` - stdout of each of the 9 unit test executables
  built against `src_clean_documented/`.
- `column_model/casim_column_documented.nc` - column-model netCDF output
  built against `src_clean_documented/` (compare to
  `../column_model/output/casim_column_kgo.nc`).
