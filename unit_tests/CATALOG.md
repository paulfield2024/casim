# CASIM micro_main.F90 Call-Graph Catalog

Catalog of the functions and subroutines reachable (directly or transitively)
from `shipway_microphysics` (module `micro_main`, [src/micro_main.F90](../src/micro_main.F90)),
excluding UM/parent-model glue code and the `casdiags` diagnostic-gathering
code in `generic_diagnostic_variables.F90` (already covered by
`column_model/casim_column_model.F90` / `/memories/repo/casim_build.md`).

This is used to drive which routines get standalone Fortran unit tests
under `unit_tests/`, with KGO output captured under `kgo/`.

## 1. Pure functions (return a value)

| File | Function | Signature | Returns | Description |
|---|---|---|---|---|
| [special.F90](../src/special.F90) | `gammalookup(x)` | `real(wp) x` | `real(wp)` | Linear-interpolation lookup table version of the Gamma function (this is what the generic `Gammafunc` interface actually resolves to). |
| special.F90 | `gammafunc1(x)` | `real(wp) x` | `real(wp)` | Gamma function via recursion + Stirling's asymptotic series (not currently bound into the `Gammafunc` interface, but usable directly). |
| special.F90 | `erfg(x, c)` | `real(wp) x`, `integer c` (0=erf,1=erfc) | `real(wp)` | Series/asymptotic error function or complementary error function. |
| special.F90 | `casim_erfc(x)` | `real(wp) x` | `real(wp)` | `erfc(x)`, implemented as `erfg(x, 1)`. |
| special.F90 | `erfinv1(x)` (interface `erfinv`) | `real(wp) x` | `real(wp)` | Inverse error function via truncated series (documented in-source as needing accuracy improvements). |
| special.F90 | `erfinv2(x)` | `real(wp) x` | `real(wp)` | Inverse error function via secant iteration (not bound to the `erfinv` generic interface). |
| special.F90 | `erfinv3(x, [tol])` | `real(wp) x`, optional `real(wp) tol` | `real(wp)` | Inverse error function via Newton-Raphson (falls back to `erfinv1` for \|x\|>0.95). |
| [qsat_casim_func.F90](../src/qsat_casim_func.F90) (module `qsat_funs`) | `Qsaturation(T, p)` | `real(wp) T` (K), `real(wp) p` (mb) | `real(wp)` | Tetens formula saturation mixing ratio over liquid water; returns 999.0 outside valid domain. |
| qsat_casim_func.F90 | `Qisaturation(T, p)` | same | `real(wp)` | Tetens formula saturation mixing ratio over ice. |
| qsat_casim_func.F90 | `dqwsatdt(qsat, T)` | `real(wp) qsat`, `real(wp) T` | `real(wp)` | d(qsat)/dT via Clausius-Clapeyron-derived closed form. |
| [gauss_4A_func.F90](../src/gauss_4A_func.F90) (module `gauss_casim_micro`) | `gauss_casim_func(a, b)` | `real(wp) a, b` | `real(wp)` | Brute-force double numerical integral used in the snow-aggregation collision kernel. |
| [lognormal_funcs.F90](../src/lognormal_funcs.F90) | `MNtoRm(M, N, density, sigma)` | 4x `real(wp)` | `real(wp)` | Inverts lognormal-mode mass/number to a mean radius; returns 0 below mass/number thresholds. |
| [lookup.F90](../src/lookup.F90) | `Gfunc(mu, p1, p2, p3)` | 4x `real(wp)` | `real(wp)` | Gamma-function combination used for 3rd-moment increments; calls `Gammafunc`. |
| lookup.F90 | `Hfunc(m1, m2, m3, p1, p2, p3)` | 6x `real(wp)` | `real(wp)` | Moment power-law combination for 3-moment inversion. **Not in module's public list** - not directly testable from outside `lookup.F90`. |
| lookup.F90 | `moment(n0, lam, mu, p)` | 4x `real(wp)` | `real(wp)` | p-th moment of a (modified) gamma size distribution. |
| [sweepout_rate.F90](../src/sweepout_rate.F90) | `sweepout(n0, lam, mu, params, rho, [mass_weight])` | `real(wp)` x3, `type(hydro_params)`, `real(wp)`, optional `logical` | `real(wp)` | Single-species collection/sweep-out rate (needs `mphys_parameters::hydro_params`). |
| sweepout_rate.F90 | `sweepout_1M2M(...)` | as above, no `mu` | `real(wp)` | Fast 1M/2M version using cached Gamma values in `params`. |
| sweepout_rate.F90 | `binary_collection(...)` | X/Y distribution params + `hydro_params` x2 + rho | `real(wp)` | Two-species collision kernel. |
| sweepout_rate.F90 | `binary_collection_1M2M(...)` | as above, no `mu` | `real(wp)` | Fast 1M/2M two-species collision kernel. |

`m3_incs.F90` (`m3_inc_type1..4`) and `which_mode_to_use.F90` (`which_mode`)
are algebraic **subroutines** (multiple `intent(out)` args) rather than pure
functions, listed under Subroutines below even though they contain no loops
over levels.

## 2. Subroutines (physics processes), in approximate call order

| File | Subroutine | Process | Direct callees |
|---|---|---|---|
| [micro_main.F90](../src/micro_main.F90) | `shipway_microphysics` | 3D grid entry point | `set_passive_fields`, `microphysics_common`, `casim_reflec`, diagnostics |
| micro_main.F90 | `microphysics_common` | per-column physics orchestration (substeps + sedimentation) | everything below |
| micro_main.F90 | `update_q` | time-integrate `qfields += tend*dt`, optional negative-clipping | - |
| micro_main.F90 | `gather_process_diagnostics` | copy `procs` rates into `casdiags` | - |
| [preconditioning.F90](../src/preconditioning.F90) | `preconditioner` | flag columns/levels with any work to do | - |
| [autoconversion.F90](../src/autoconversion.F90) | `raut` | cloud→rain autoconversion (KK2000 or Kogan2013) | - |
| [accretion.F90](../src/accretion.F90) | `racw` | rain accretes cloud | `sweepout` (if not using KK/Kogan empirical form) |
| [aggregation.F90](../src/aggregation.F90) | `racr` | rain self-collection | - |
| aggregation.F90 | `ice_aggregation` | snow self-collection | `gaussfunclookup` → `gauss_casim_func` |
| [snow_autoconversion.F90](../src/snow_autoconversion.F90) | `saut` | ice→snow autoconversion (Bergeron) | - |
| [evaporation.F90](../src/evaporation.F90) | `revp` | rain evaporation + CCN re-activation bookkeeping | `ventilation_3M`/`ventilation_1M_2M`, `which_mode` → `MNtoRm`, `Qsaturation` |
| [ice_nucleation.F90](../src/ice_nucleation.F90) | `inuc` | heterogeneous ice nucleation (5 schemes) | `Qsaturation`, `Qisaturation` |
| [ice_accretion.F90](../src/ice_accretion.F90) | `iacc` | binary accretion between hydrometeor types (→graupel/snow) | `binary_collection`/`binary_collection_1M2M` |
| [graupel_embryo.F90](../src/graupel_embryo.F90) | `graupel_embryos` | small snow + cloud → graupel embryos | - |
| [graupel_wetgrowth.F90](../src/graupel_wetgrowth.F90) | `wetgrowth` | graupel wet growth/shedding | - |
| [breakup.F90](../src/breakup.F90) | `ice_breakup` | snow mechanical breakup | - |
| [ice_multiplication.F90](../src/ice_multiplication.F90) | `hallet_mossop` | Hallet-Mossop secondary ice production | - |
| ice_multiplication.F90 | `droplet_shattering` | rime splintering | - |
| ice_multiplication.F90 | `ice_collision` | ice-ice collisional fragmentation | - |
| [ice_deposition.F90](../src/ice_deposition.F90) | `idep` | vapour deposition growth (ice/snow/graupel) | `ventilation_3M`/`ventilation_1M_2M` |
| [adjust_deposition.F90](../src/adjust_deposition.F90) | `adjust_dep` | rescale deposition to stay physically consistent | - |
| [ice_melting.F90](../src/ice_melting.F90) | `melting` | ice/snow/graupel melting | `ventilation_3M`/`ventilation_1M_2M` |
| [homogeneous_freezing.F90](../src/homogeneous_freezing.F90) | `ihom_rain` | homogeneous freezing of rain (T<-40C) | - |
| homogeneous_freezing.F90 | `ihom_droplets` | homogeneous freezing of cloud droplets | - |
| [ventfac.F90](../src/ventfac.F90) | `ventilation_3M` | ventilation factor, 3-moment | `Gammafunc` |
| ventfac.F90 | `ventilation_1M_2M` | ventilation factor, 1M/2M (cached Gamma values) | - |
| [distributions.F90](../src/distributions.F90) | `query_distributions` | derive (n0,lam,mu) from moments | `get_slope_generic` |
| [lookup.F90](../src/lookup.F90) | `get_slope_generic` | dispatch to 1M/2M/3M inversion | `get_lam_n0_1M/2M/3M` |
| lookup.F90 | `get_lam_n0_3M` | 3-moment inversion | `get_mu` |
| lookup.F90 | `get_lam_n0_2M` | 2-moment inversion (fixed mu) | `Gammafunc` |
| lookup.F90 | `get_lam_n0_1M` | 1-moment inversion (diagnostic) | - |
| lookup.F90 | `get_mu` | iterative solve for shape parameter mu given 3 moments | - |
| lookup.F90 | `get_n0` | back out intercept n0 given a moment, power, mu, lam | `Gammafunc` |
| [m3_incs.F90](../src/m3_incs.F90) | `m3_inc_type4` | 3rd-moment increment on phase change (algebraic) | - |
| m3_incs.F90 | `m3_inc_type3` | 3rd-moment increment from dm1,dm2 assuming fixed mu | `Gfunc` |
| m3_incs.F90 | `m3_inc_type2` | 3rd-moment increment given existing m1,m2,m3 | `m3_inc_type3` (fallback if m1<=0) |
| m3_incs.F90 | `m3_inc_type1` | simplest number increment scaling: `dm2 = dm1*m1/m2` | - |
| [which_mode_to_use.F90](../src/which_mode_to_use.F90) | `which_mode` | partition evaporated aerosol between CCN modes | `MNtoRm` |
| [sum_procs.F90](../src/sum_procs.F90) (module `sum_process`) | `sum_procs` | accumulate hydrometeor process rates into tendencies (+ latent heating) | `m3_inc_type*` (if 3rd moment requested) |
| sum_procs.F90 | `sum_aprocs` | accumulate aerosol process rates into tendencies | - |
| [mphys_tidy.F90](../src/mphys_tidy.F90) | `tidy_qin` | enforce hydrometeor moment positivity/consistency | - |
| mphys_tidy.F90 | `tidy_ain` | enforce aerosol moment positivity/consistency | - |
| mphys_tidy.F90 | `ensure_positive` | rescale a group of processes to avoid negative moments | - |
| mphys_tidy.F90 | `ensure_positive_aerosol` | aerosol version of `ensure_positive` | - |
| mphys_tidy.F90 | `ensure_saturated` | rescale deposition/sublimation to not over-remove vapour | - |
| [condensation.F90](../src/condensation.F90) | `condevp` | condensation/evaporation + droplet activation | `Qsaturation`, `Qisaturation`, (Shipway activation tables if `iopt_act` selects that scheme) |
| [sedimentation.F90](../src/sedimentation.F90) | `sedr` | Eulerian flux sedimentation, full 3-moment | - |
| sedimentation.F90 | `sedr_1M_2M` | fast sedimentation, 1M/2M (cached Gamma values) | - |
| sedimentation.F90 | `terminal_velocity_CFL` | choose sedimentation substeps to satisfy CFL | - |
| [process_routines.F90](../src/process_routines.F90) | `allocate_procs`/`deallocate_procs`/`zero_procs` | process-rate bookkeeping array management | - |
| [aerosol_routines.F90](../src/aerosol_routines.F90) | `examine_aerosol` | derive aerosol mode radii/activation partitioning (only exercised if `l_process`, i.e. `aerosol_option>0`) | `MNtoRm` |
| aerosol_routines.F90 | `AbdulRazzakGhan2000[_dust]` | CCN/dust activation supersaturation scheme (only if `iopt_act` selects it) | `casim_erfc` |

Note: `aerosol_routines.F90` and the STASH aerosol diagnostics are only
exercised when `l_process = .true.` (i.e. `aerosol_option > 0`); with this
repo's default column-model config (`aerosol_option = 0`) that code path is
dead at runtime (see `/memories/repo/casim_build.md` for the confirmed
SIGSEGV if the corresponding STASH diagnostic flags are force-enabled
without `aerosol_option > 0`).

## Testing priority

Pure functions with simple scalar signatures (no derived types) are the
easiest and highest-value first targets for `unit_tests/`:
`Gammafunc`/`gammafunc1`, `casim_erfc`, `erfinv`, `Qsaturation`,
`Qisaturation`, `dqwsatdt`, `gauss_casim_func`, `MNtoRm`, and the `m3_incs`
family (`m3_inc_type1/3/4` are pure algebra; `m3_inc_type2` calls
`m3_inc_type3` as a fallback). `Gfunc`/`Hfunc`/`moment` and the
`sweepout_rate.F90`/`lookup.F90` distribution-inversion routines need a
`type(hydro_params)` (from `mphys_parameters`, populated by
`set_mphys_switches`) and are next in line. The full per-process subroutines
(`raut`, `racw`, `revp`, `inuc`, ...) need populated `qfields`/`cffields`/
`procs` arrays and are best tested via `microphysics_common`/
`shipway_microphysics` integration (already covered by the column model) or
with hand-built minimal column arrays.

## Subroutine tests written (in addition to the pure-function tests above)

Beyond the pure functions, three subroutines with plain real(wp)/simple-array
signatures (no `type(hydro_params)` needed) were also given standalone tests:

- `which_mode` (`which_mode_to_use.F90`, `unit_tests/test_which_mode.F90`):
  with the repo's default `l_aeroproc_midway=.FALSE.`, `imethod` is a
  compile-time parameter fixed to `isimple_method`, so only the simple
  rm-vs-0.5-micron-threshold branch is reachable at runtime; tested both
  sides of the threshold plus the `dm*dn<=0` no-op case.
- `m3_inc_type2` (`m3_incs.F90`, `unit_tests/test_m3_inc_type2.F90`): tests
  both the `m1>0` main branch (using the hardwired `p3=6` shortcut formula)
  and the `m1<=0` fallback to `m3_inc_type3` (reusing the already-verified
  `Gfunc(0,0,3,6)=20` case from `test_m3_incs.F90`).
- `get_n0` and `get_mu` (`lookup.F90`, added to `unit_tests/test_lookup.F90`):
  `get_n0` is plain algebra. `get_mu` needs a (mu -> Gfunc value) lookup
  table; built locally via the public `set_mu_lookup` helper (same call
  signature used by `initialize.F90`) rather than depending on module-level
  `mu_i`/`mu_g` state, then solved for a known target mu by picking
  `m1=1, m2=1` so `Hfunc(m1,m2,m3,p1,p2,p3)=m3` exactly equals the target
  Gfunc value.

The remaining physics-process subroutines (`raut`, `racw`, `revp`, `inuc`,
sedimentation, etc.) still require populated `qfields`/`cffields`/`procs`
arrays or a full `type(hydro_params)` and are left to the column-model
integration test rather than isolated unit tests (see
`column_model/casim_column_model.F90` and `column_model/output/casim_column_kgo.nc`).
