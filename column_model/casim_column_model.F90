! =============================================================================
! casim_column_model.F90
!
! A minimal 1D (single column) vertical atmospheric wrapper around CASIM's
! `shipway_microphysics` entry point (src/micro_main.F90). This program
! builds an idealised vertical profile, then repeatedly calls CASIM to
! integrate the microphysics forward in time, writing the prognostic and
! diagnostic fields out to plain-text files (one row per timestep, one
! column per vertical level). `make_netcdf.py` subsequently converts these
! text files into a single NetCDF file with (time, height) dimensions,
! which serves as the Known-Good-Output (KGO) for this benchmark.
!
! Default case (as specified in prompts/make_benchmark_tests):
!   - 0-10 km column, dz = 100 m (100 levels, cell centres)
!   - surface temperature 15 degC, lapse rate 6.5 K/km (dry hydrostatic
!     background atmosphere, fixed in time - CASIM is a microphysics-only
!     scheme with no dynamical core, so the thermodynamic background state
!     is held fixed while theta/qv/hydrometeors evolve under microphysics)
!   - water vapour density 10 g/m3 up to 1 km, decreasing linearly by
!     1 g/m3 per km above that (converted to specific humidity using the
!     local air density)
!   - parabolic vertical velocity, max 5 m/s at 5 km, zero at 0 and 10 km
!   - dt = 10 s for 500 s (50 steps)
!   - full double-moment configuration (option 22222), no aerosol tracers
!     (aerosol_option = 0)
! =============================================================================
program casim_column_model

  use variable_precision, only: wp
  use casim_parent_mod,    only: casim_parent, parent_monc
  use mphys_switches,      only: set_mphys_switches
  use mphys_constants,     only: Rd, cp, g
  use initialize,          only: mphys_init, mphys_finalise
  use micro_main,          only: shipway_microphysics
  use generic_diagnostic_variables, only: casdiags, allocate_diagnostic_space, &
                                           deallocate_diagnostic_space

  implicit none

  ! -------------------------------------------------------------------------
  ! Column / time configuration
  ! -------------------------------------------------------------------------
  integer,  parameter :: nz          = 100         ! number of vertical levels
  real(wp), parameter :: layer_dz    = 100.0_wp    ! m
  real(wp), parameter :: dt          = 10.0_wp     ! s
  real(wp), parameter :: run_length  = 500.0_wp    ! s
  integer,  parameter :: nt          = nint(run_length/dt)

  ! Idealised background atmosphere parameters
  real(wp), parameter :: Tsurf_C  = 15.0_wp        ! degC
  real(wp), parameter :: lapse    = 6.5e-3_wp       ! K / m
  real(wp), parameter :: p0       = 101325.0_wp     ! Pa
  real(wp), parameter :: T0K      = Tsurf_C + 273.15_wp

  ! Moisture profile parameters (interpreted as vapour density, kg/m3)
  real(wp), parameter :: rhov_low   = 10.0e-3_wp    ! kg/m3 below 1 km
  real(wp), parameter :: rhov_lapse = 1.0e-3_wp     ! kg/m3 per km above 1 km
  real(wp), parameter :: z_moist_top = 1000.0_wp     ! m

  ! Vertical velocity profile parameters
  real(wp), parameter :: w_max = 5.0_wp             ! m/s
  real(wp), parameter :: z_wmax = 5000.0_wp         ! m
  real(wp), parameter :: z_top  = 10000.0_wp        ! m

  ! -------------------------------------------------------------------------
  ! Column state arrays. CASIM's shipway_microphysics expects (k,i,j) shaped
  ! arrays; here i=j=1 (single column).
  ! -------------------------------------------------------------------------
  real(wp) :: z(nz)
  real(wp) :: exner(nz,1,1), pressure(nz,1,1), rho(nz,1,1), dz(nz,1,1)
  real(wp) :: w(nz,1,1), tke(nz,1,1)
  real(wp) :: theta(nz,1,1), qv(nz,1,1)

  ! Hydrometeor prognostics: q1=ql, q2=qr, q3=nl, q4=nr, q5=m3r,
  !                          q6=qi, q7=qs, q8=qg, q9=ni, q10=ns, q11=ng,
  !                          q12=m3s, q13=m3g
  real(wp) :: q1(nz,1,1), q2(nz,1,1), q3(nz,1,1), q4(nz,1,1), q5(nz,1,1)
  real(wp) :: q6(nz,1,1), q7(nz,1,1), q8(nz,1,1), q9(nz,1,1), q10(nz,1,1)
  real(wp) :: q11(nz,1,1), q12(nz,1,1), q13(nz,1,1)

  ! Aerosol fields in (unused: aerosol_option = 0)
  real(wp) :: a1(nz,1,1), a2(nz,1,1), a3(nz,1,1), a4(nz,1,1), a5(nz,1,1)
  real(wp) :: a6(nz,1,1), a7(nz,1,1), a8(nz,1,1), a9(nz,1,1), a10(nz,1,1)
  real(wp) :: a11(nz,1,1), a12(nz,1,1), a13(nz,1,1), a14(nz,1,1), a15(nz,1,1)
  real(wp) :: a16(nz,1,1), a17(nz,1,1), a18(nz,1,1), a19(nz,1,1), a20(nz,1,1)

  ! Cloud fractions (in/out)
  real(wp) :: cfliq(nz,1,1), cfice(nz,1,1), cfsnow(nz,1,1), cfrain(nz,1,1), cfgr(nz,1,1)

  ! Tendency work arrays: zeroed before each call (no external forcing),
  ! returned by CASIM as the per-timestep increment (l_tendency = .false.)
  real(wp) :: dqv(nz,1,1), dth(nz,1,1)
  real(wp) :: dq1(nz,1,1), dq2(nz,1,1), dq3(nz,1,1), dq4(nz,1,1), dq5(nz,1,1)
  real(wp) :: dq6(nz,1,1), dq7(nz,1,1), dq8(nz,1,1), dq9(nz,1,1), dq10(nz,1,1)
  real(wp) :: dq11(nz,1,1), dq12(nz,1,1), dq13(nz,1,1)
  real(wp) :: da1(nz,1,1), da2(nz,1,1), da3(nz,1,1), da4(nz,1,1), da5(nz,1,1)
  real(wp) :: da6(nz,1,1), da7(nz,1,1), da8(nz,1,1), da9(nz,1,1), da10(nz,1,1)
  real(wp) :: da11(nz,1,1), da12(nz,1,1), da13(nz,1,1), da14(nz,1,1), da15(nz,1,1)
  real(wp) :: da16(nz,1,1), da17(nz,1,1)

  real(wp) :: Tk, p_here, rhov
  integer  :: k, n

  ! Output file units, organised as (name, unit)
  integer, parameter :: nvars_prog = 15  ! qv,theta,q1..q13
  integer, parameter :: nvars_diag = 87  ! profile diagnostics: see list below
  integer, parameter :: nvars_scalar = 9 ! surface precip rates + water paths + surface_cloud
  character(len=32) :: prog_names(nvars_prog)
  character(len=32) :: diag_names(nvars_diag)
  character(len=32) :: scalar_names(nvars_scalar)
  integer :: prog_units(nvars_prog), diag_units(nvars_diag), scalar_units(nvars_scalar)
  integer :: time_unit, height_unit
  character(len=*), parameter :: outdir = 'output/'

  ! ---------------------------------------------------------------------
  ! 1. Build the idealised background atmosphere and initial conditions
  ! ---------------------------------------------------------------------
  do k = 1, nz
    z(k) = (real(k,wp) - 0.5_wp) * layer_dz

    Tk = T0K - lapse*z(k)
    ! Exact hydrostatic solution for a constant lapse-rate dry atmosphere:
    ! p(z) = p0 * (T(z)/T0)^(g/(Rd*lapse))
    p_here = p0 * (Tk/T0K) ** (g/(Rd*lapse))

    pressure(k,1,1) = p_here
    exner(k,1,1)    = (p_here/p0) ** (Rd/cp)
    theta(k,1,1)    = Tk / exner(k,1,1)
    rho(k,1,1)      = p_here / (Rd*Tk)
    dz(k,1,1)       = layer_dz

    ! Vapour density profile -> specific humidity via local air density
    if (z(k) <= z_moist_top) then
      rhov = rhov_low
    else
      rhov = max(0.0_wp, rhov_low - rhov_lapse*(z(k)-z_moist_top)/1000.0_wp)
    end if
    qv(k,1,1) = rhov / rho(k,1,1)

    ! Parabolic vertical velocity, max w_max at z_wmax, zero at 0 and z_top
    w(k,1,1) = w_max * (1.0_wp - ((z(k)-z_wmax)/z_wmax)**2)

    tke(k,1,1) = 0.0_wp
  end do

  ! No initial hydrometeors.
  q1=0.0_wp; q2=0.0_wp; q3=0.0_wp; q4=0.0_wp; q5=0.0_wp
  q6=0.0_wp; q7=0.0_wp; q8=0.0_wp; q9=0.0_wp; q10=0.0_wp
  q11=0.0_wp; q12=0.0_wp; q13=0.0_wp

  ! Cloud fraction fields. This column model has no prognostic cloud
  ! fraction scheme, so assume the grid box is fully overcast (cf=1)
  ! wherever hydrometeors occur - this matches the "assumed cloud
  ! fraction is 1" behaviour that CASIM's own radar reflectivity code
  ! (casim_reflec_mod) warns about when no cloud-fraction scheme is
  ! active. NOTE: leaving these at 0 (as originally written) makes
  ! casim_reflec_mod's dBZ output multiply by zero unconditionally
  ! (it has no l_prf_cfrac floor like most other process routines do),
  ! so radar reflectivity would be stuck at the ref_lim=-35dBZ floor
  ! regardless of actual hydrometeor content.
  cfliq=1.0_wp; cfice=1.0_wp; cfsnow=1.0_wp; cfrain=1.0_wp; cfgr=1.0_wp

  ! No aerosol tracers (aerosol_option = 0)
  a1=0.0_wp; a2=0.0_wp; a3=0.0_wp; a4=0.0_wp; a5=0.0_wp
  a6=0.0_wp; a7=0.0_wp; a8=0.0_wp; a9=0.0_wp; a10=0.0_wp
  a11=0.0_wp; a12=0.0_wp; a13=0.0_wp; a14=0.0_wp; a15=0.0_wp
  a16=0.0_wp; a17=0.0_wp; a18=0.0_wp; a19=0.0_wp; a20=0.0_wp

  ! ---------------------------------------------------------------------
  ! 2. Configure and initialise CASIM
  ! ---------------------------------------------------------------------
  casim_parent = parent_monc

  ! option = 22222 -> double-moment cloud, rain, ice, snow, graupel
  ! aerosol_option = 0 -> no aerosol tracers (fixed CCN/IN assumptions)
  call set_mphys_switches(22222, 0)

  ! Request the diagnostics we want to write out
  casdiags % l_surface_rain  = .true.
  casdiags % l_surface_snow  = .true.
  casdiags % l_surface_graup = .true.
  casdiags % l_rainfall_3d   = .true.
  casdiags % l_snowfall_3d   = .true.
  casdiags % l_snowonly_3d   = .true.
  casdiags % l_graupfall_3d  = .true.
  casdiags % l_dth = .true.
  casdiags % l_dqv = .true.
  casdiags % l_dqc = .true.
  casdiags % l_dqr = .true.
  casdiags % l_dqi = .true.
  casdiags % l_dqs = .true.
  casdiags % l_dqg = .true.
  casdiags % l_lwp = .true.
  casdiags % l_rwp = .true.
  casdiags % l_iwp = .true.
  casdiags % l_swp = .true.
  casdiags % l_gwp = .true.
  casdiags % l_surface_cloud = .true.
  casdiags % l_radar = .true.

  ! Process-rate and number-tendency diagnostics (procs structure), see
  ! src/generic_diagnostic_variables.F90 for full descriptions.
  casdiags % l_phomc = .true.
  casdiags % l_pinuc = .true.
  casdiags % l_pidep = .true.
  casdiags % l_psdep = .true.
  casdiags % l_piacw = .true.
  casdiags % l_psacw = .true.
  casdiags % l_psacr = .true.
  casdiags % l_pisub = .true.
  casdiags % l_pssub = .true.
  casdiags % l_pimlt = .true.
  casdiags % l_psmlt = .true.
  casdiags % l_psaut = .true.
  casdiags % l_psaci = .true.
  casdiags % l_praut = .true.
  casdiags % l_pracw = .true.
  casdiags % l_prevp = .true.
  casdiags % l_pgacw = .true.
  casdiags % l_pgacs = .true.
  casdiags % l_pgmlt = .true.
  casdiags % l_pgsub = .true.
  casdiags % l_psedi = .true.
  casdiags % l_pseds = .true.
  casdiags % l_psedr = .true.
  casdiags % l_psedg = .true.
  casdiags % l_psedl = .true.
  casdiags % l_pcond = .true.
  casdiags % l_phomr = .true.
  casdiags % l_nhomc = .true.
  casdiags % l_nhomr = .true.
  casdiags % l_nihal = .true.
  casdiags % l_ninuc = .true.
  casdiags % l_nsedi = .true.
  casdiags % l_nseds = .true.
  casdiags % l_nsedg = .true.
  casdiags % l_nraut = .true.
  casdiags % l_nsedl = .true.
  casdiags % l_nracw = .true.
  casdiags % l_nracr = .true.
  casdiags % l_nsedr = .true.
  casdiags % l_nrevp = .true.
  casdiags % l_nisub = .true.
  casdiags % l_nssub = .true.
  casdiags % l_nsaut = .true.
  casdiags % l_nsaci = .true.
  casdiags % l_ngacs = .true.
  casdiags % l_ngsub = .true.
  casdiags % l_niacw = .true.
  casdiags % l_nsacw = .true.
  casdiags % l_nsacr = .true.
  casdiags % l_nimlt = .true.
  casdiags % l_nsmlt = .true.
  casdiags % l_ngacw = .true.
  casdiags % l_ngmlt = .true.
  casdiags % l_pihal = .true.
  casdiags % l_praci_g = .true.
  casdiags % l_praci_r = .true.
  casdiags % l_praci_i = .true.
  casdiags % l_nraci_g = .true.
  casdiags % l_nraci_r = .true.
  casdiags % l_nraci_i = .true.
  casdiags % l_pidps = .true.
  casdiags % l_nidps = .true.
  casdiags % l_pgaci = .true.
  casdiags % l_ngaci = .true.
  casdiags % l_niics_s = .true.
  casdiags % l_niics_i = .true.

  ! Remaining diagnostics not otherwise needed for the KGO text/NetCDF
  ! output, enabled anyway per request ("all diagnostic logicals on
  ! unless they cause a crash"). l_precip is allocated but never
  ! populated by shipway_microphysics (harmless, stays zero).
  casdiags % l_precip = .true.

  ! l_mphys_pts allocates a LOGICAL (not REAL) array - safe to enable,
  ! but intentionally not wired into write_diagnostics_* since it is
  ! not a real-valued field.
  casdiags % l_mphys_pts = .true.

  ! NOTE: the aerosol-processing STASH diagnostics (l_aact_*, l_aaut,
  ! l_asedr_*, l_dnuc_*, etc., codes 601-676) are intentionally left
  ! .FALSE. here. They DO cause a crash with this configuration: their
  ! population code in micro_main.F90's gather_process_diagnostics is
  ! gated on the aswitch%l_XXX process-enable switches (which default
  ! to .true. regardless of aerosol_option), not on l_process. With
  ! aerosol_option = 0, l_process is .false., so the aerosol process
  ! indices (e.g. i_aact%id) are never allocated (stay 0), and
  ! gather_process_diagnostics then indexes aerosol_procs(:, 0, :),
  ! which is out of bounds -> SIGSEGV. Confirmed via a debug build
  ! (-fcheck=all): "Index '0' of dimension 2 of array 'aerosol_procs'
  ! below lower bound of 1" at micro_main.F90:3360.

  call allocate_diagnostic_space(1, 1, 1, 1, 1, nz)

  call mphys_init(1, 1, 1, 1, 1, nz,                                     &
       is_in=1, ie_in=1, js_in=1, je_in=1, ks_in=1, ke_in=nz,             &
       l_tendency=.false.)

  ! ---------------------------------------------------------------------
  ! 3. Open output files and write the initial state (t = 0)
  ! ---------------------------------------------------------------------
  call execute_command_line('mkdir -p '//outdir)

  prog_names(1)  = 'qv'
  prog_names(2)  = 'theta'
  prog_names(3)  = 'q1_ql'
  prog_names(4)  = 'q2_qr'
  prog_names(5)  = 'q3_nl'
  prog_names(6)  = 'q4_nr'
  prog_names(7)  = 'q5_m3r'
  prog_names(8)  = 'q6_qi'
  prog_names(9)  = 'q7_qs'
  prog_names(10) = 'q8_qg'
  prog_names(11) = 'q9_ni'
  prog_names(12) = 'q10_ns'
  prog_names(13) = 'q11_ng'
  prog_names(14) = 'q12_m3s'
  prog_names(15) = 'q13_m3g'

  diag_names(1)  = 'rainfall_3d'
  diag_names(2)  = 'snowfall_3d'
  diag_names(3)  = 'snowonly_3d'
  diag_names(4)  = 'graupfall_3d'
  diag_names(5)  = 'dth_total'
  diag_names(6)  = 'dqv_total'
  diag_names(7)  = 'dqc'
  diag_names(8)  = 'dqr'
  diag_names(9)  = 'dqi'
  diag_names(10) = 'dqs'
  diag_names(11) = 'dqg'
  diag_names(12) = 'w'
  diag_names(13) = 'rho'
  diag_names(14) = 'pressure'
  diag_names(15) = 'exner'
  diag_names(16) = 'dbz_tot'
  diag_names(17) = 'dbz_g'
  diag_names(18) = 'dbz_i'
  diag_names(19) = 'dbz_s'
  diag_names(20) = 'dbz_l'
  diag_names(21) = 'dbz_r'
  diag_names(22) = 'phomc'
  diag_names(23) = 'pinuc'
  diag_names(24) = 'pidep'
  diag_names(25) = 'psdep'
  diag_names(26) = 'piacw'
  diag_names(27) = 'psacw'
  diag_names(28) = 'psacr'
  diag_names(29) = 'pisub'
  diag_names(30) = 'pssub'
  diag_names(31) = 'pimlt'
  diag_names(32) = 'psmlt'
  diag_names(33) = 'psaut'
  diag_names(34) = 'psaci'
  diag_names(35) = 'praut'
  diag_names(36) = 'pracw'
  diag_names(37) = 'prevp'
  diag_names(38) = 'pgacw'
  diag_names(39) = 'pgacs'
  diag_names(40) = 'pgmlt'
  diag_names(41) = 'pgsub'
  diag_names(42) = 'psedi'
  diag_names(43) = 'pseds'
  diag_names(44) = 'psedr'
  diag_names(45) = 'psedg'
  diag_names(46) = 'psedl'
  diag_names(47) = 'pcond'
  diag_names(48) = 'phomr'
  diag_names(49) = 'nhomc'
  diag_names(50) = 'nhomr'
  diag_names(51) = 'nihal'
  diag_names(52) = 'ninuc'
  diag_names(53) = 'nsedi'
  diag_names(54) = 'nseds'
  diag_names(55) = 'nsedg'
  diag_names(56) = 'nraut'
  diag_names(57) = 'nsedl'
  diag_names(58) = 'nracw'
  diag_names(59) = 'nracr'
  diag_names(60) = 'nsedr'
  diag_names(61) = 'nrevp'
  diag_names(62) = 'nisub'
  diag_names(63) = 'nssub'
  diag_names(64) = 'nsaut'
  diag_names(65) = 'nsaci'
  diag_names(66) = 'ngacs'
  diag_names(67) = 'ngsub'
  diag_names(68) = 'niacw'
  diag_names(69) = 'nsacw'
  diag_names(70) = 'nsacr'
  diag_names(71) = 'nimlt'
  diag_names(72) = 'nsmlt'
  diag_names(73) = 'ngacw'
  diag_names(74) = 'ngmlt'
  diag_names(75) = 'pihal'
  diag_names(76) = 'praci_g'
  diag_names(77) = 'praci_r'
  diag_names(78) = 'praci_i'
  diag_names(79) = 'nraci_g'
  diag_names(80) = 'nraci_r'
  diag_names(81) = 'nraci_i'
  diag_names(82) = 'pidps'
  diag_names(83) = 'nidps'
  diag_names(84) = 'pgaci'
  diag_names(85) = 'ngaci'
  diag_names(86) = 'niics_s'
  diag_names(87) = 'niics_i'

  scalar_names(1) = 'surface_rain'
  scalar_names(2) = 'surface_snow'
  scalar_names(3) = 'surface_graup'
  scalar_names(4) = 'lwp'
  scalar_names(5) = 'rwp'
  scalar_names(6) = 'iwp'
  scalar_names(7) = 'swp'
  scalar_names(8) = 'gwp'
  scalar_names(9) = 'surface_cloud'

  open(newunit=height_unit, file=outdir//'height.txt', status='replace')
  do k = 1, nz
    write(height_unit,'(F12.4)') z(k)
  end do
  close(height_unit)

  open(newunit=time_unit, file=outdir//'time.txt', status='replace')

  do n = 1, nvars_prog
    open(newunit=prog_units(n), file=outdir//trim(prog_names(n))//'.txt', status='replace')
  end do
  do n = 1, nvars_diag
    open(newunit=diag_units(n), file=outdir//trim(diag_names(n))//'.txt', status='replace')
  end do
  do n = 1, nvars_scalar
    open(newunit=scalar_units(n), file=outdir//trim(scalar_names(n))//'.txt', status='replace')
  end do

  call write_time(time_unit, 0.0_wp)
  call write_prognostics()
  call write_diagnostics_initial()
  call write_scalars_initial()

  ! ---------------------------------------------------------------------
  ! 4. Time integration loop
  ! ---------------------------------------------------------------------
  do n = 1, nt

    dqv=0.0_wp; dth=0.0_wp
    dq1=0.0_wp; dq2=0.0_wp; dq3=0.0_wp; dq4=0.0_wp; dq5=0.0_wp
    dq6=0.0_wp; dq7=0.0_wp; dq8=0.0_wp; dq9=0.0_wp; dq10=0.0_wp
    dq11=0.0_wp; dq12=0.0_wp; dq13=0.0_wp
    da1=0.0_wp; da2=0.0_wp; da3=0.0_wp; da4=0.0_wp; da5=0.0_wp
    da6=0.0_wp; da7=0.0_wp; da8=0.0_wp; da9=0.0_wp; da10=0.0_wp
    da11=0.0_wp; da12=0.0_wp; da13=0.0_wp; da14=0.0_wp; da15=0.0_wp
    da16=0.0_wp; da17=0.0_wp

    call shipway_microphysics(1, 1, 1, 1, 1, nz, dt,                     &
         qv, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13,      &
         theta, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13,   &
         a14, a15, a16, a17, a18, a19, a20,                               &
         exner, pressure, rho, w, tke, dz,                                &
         cfliq, cfice, cfsnow, cfrain, cfgr,                              &
         dqv, dq1, dq2, dq3, dq4, dq5, dq6, dq7, dq8, dq9, dq10, dq11,    &
         dq12, dq13, dth, da1, da2, da3, da4, da5, da6, da7, da8, da9,    &
         da10, da11, da12, da13, da14, da15, da16, da17,                  &
         is_in=1, ie_in=1, js_in=1, je_in=1)

    ! l_tendency = .false. => dq*/dth/dqv are already increments over dt
    qv    = qv    + dqv
    theta = theta + dth
    q1=q1+dq1; q2=q2+dq2; q3=q3+dq3; q4=q4+dq4; q5=q5+dq5
    q6=q6+dq6; q7=q7+dq7; q8=q8+dq8; q9=q9+dq9; q10=q10+dq10
    q11=q11+dq11; q12=q12+dq12; q13=q13+dq13

    call write_time(time_unit, real(n,wp)*dt)
    call write_prognostics()
    call write_diagnostics_step()
    call write_scalars_step()

  end do

  close(time_unit)
  do n = 1, nvars_prog
    close(prog_units(n))
  end do
  do n = 1, nvars_diag
    close(diag_units(n))
  end do
  do n = 1, nvars_scalar
    close(scalar_units(n))
  end do

  call mphys_finalise()
  call deallocate_diagnostic_space()

  print *, 'CASIM 1D column model complete: ', nt, ' steps of ', dt, ' s written to ', outdir

contains

  subroutine write_time(unit_no, t)
    integer,  intent(in) :: unit_no
    real(wp), intent(in) :: t
    write(unit_no,'(F12.4)') t
  end subroutine write_time

  subroutine write_row(unit_no, vec)
    integer,  intent(in) :: unit_no
    real(wp), intent(in) :: vec(nz)
    character(len=32) :: fmt
    write(fmt,'(A,I0,A)') '(', nz, 'ES16.7)'
    write(unit_no, fmt) vec(:)
  end subroutine write_row

  subroutine write_prognostics()
    call write_row(prog_units(1),  qv(:,1,1))
    call write_row(prog_units(2),  theta(:,1,1))
    call write_row(prog_units(3),  q1(:,1,1))
    call write_row(prog_units(4),  q2(:,1,1))
    call write_row(prog_units(5),  q3(:,1,1))
    call write_row(prog_units(6),  q4(:,1,1))
    call write_row(prog_units(7),  q5(:,1,1))
    call write_row(prog_units(8),  q6(:,1,1))
    call write_row(prog_units(9),  q7(:,1,1))
    call write_row(prog_units(10), q8(:,1,1))
    call write_row(prog_units(11), q9(:,1,1))
    call write_row(prog_units(12), q10(:,1,1))
    call write_row(prog_units(13), q11(:,1,1))
    call write_row(prog_units(14), q12(:,1,1))
    call write_row(prog_units(15), q13(:,1,1))
  end subroutine write_prognostics

  ! At t=0, no microphysics has run yet: write zeros for the
  ! process/tendency diagnostics, but the background w/rho/pressure/exner.
  subroutine write_diagnostics_initial()
    real(wp) :: zero_vec(nz)
    integer  :: m
    zero_vec = 0.0_wp
    call write_row(diag_units(1),  zero_vec)          ! rainfall_3d
    call write_row(diag_units(2),  zero_vec)          ! snowfall_3d
    call write_row(diag_units(3),  zero_vec)          ! snowonly_3d
    call write_row(diag_units(4),  zero_vec)          ! graupfall_3d
    call write_row(diag_units(5),  zero_vec)          ! dth_total
    call write_row(diag_units(6),  zero_vec)          ! dqv_total
    call write_row(diag_units(7),  zero_vec)          ! dqc
    call write_row(diag_units(8),  zero_vec)          ! dqr
    call write_row(diag_units(9),  zero_vec)          ! dqi
    call write_row(diag_units(10), zero_vec)          ! dqs
    call write_row(diag_units(11), zero_vec)          ! dqg
    call write_row(diag_units(12), w(:,1,1))          ! w
    call write_row(diag_units(13), rho(:,1,1))        ! rho
    call write_row(diag_units(14), pressure(:,1,1))   ! pressure
    call write_row(diag_units(15), exner(:,1,1))      ! exner
    do m = 16, nvars_diag
      call write_row(diag_units(m), zero_vec)
    end do
  end subroutine write_diagnostics_initial

  subroutine write_diagnostics_step()
    call write_row(diag_units(1),  real(casdiags % rainfall_3d(1,1,:), wp))
    call write_row(diag_units(2),  real(casdiags % snowfall_3d(1,1,:), wp))
    call write_row(diag_units(3),  real(casdiags % snowonly_3d(1,1,:), wp))
    call write_row(diag_units(4),  real(casdiags % graupfall_3d(1,1,:), wp))
    call write_row(diag_units(5),  real(casdiags % dth_total(1,1,:), wp))
    call write_row(diag_units(6),  real(casdiags % dqv_total(1,1,:), wp))
    call write_row(diag_units(7),  real(casdiags % dqc(1,1,:), wp))
    call write_row(diag_units(8),  real(casdiags % dqr(1,1,:), wp))
    call write_row(diag_units(9),  real(casdiags % dqi(1,1,:), wp))
    call write_row(diag_units(10), real(casdiags % dqs(1,1,:), wp))
    call write_row(diag_units(11), real(casdiags % dqg(1,1,:), wp))
    call write_row(diag_units(12), w(:,1,1))
    call write_row(diag_units(13), rho(:,1,1))
    call write_row(diag_units(14), pressure(:,1,1))
    call write_row(diag_units(15), exner(:,1,1))
    call write_row(diag_units(16), real(casdiags % dbz_tot(1,1,:), wp))
    call write_row(diag_units(17), real(casdiags % dbz_g(1,1,:), wp))
    call write_row(diag_units(18), real(casdiags % dbz_i(1,1,:), wp))
    call write_row(diag_units(19), real(casdiags % dbz_s(1,1,:), wp))
    call write_row(diag_units(20), real(casdiags % dbz_l(1,1,:), wp))
    call write_row(diag_units(21), real(casdiags % dbz_r(1,1,:), wp))
    call write_row(diag_units(22), real(casdiags % phomc(1,1,:), wp))
    call write_row(diag_units(23), real(casdiags % pinuc(1,1,:), wp))
    call write_row(diag_units(24), real(casdiags % pidep(1,1,:), wp))
    call write_row(diag_units(25), real(casdiags % psdep(1,1,:), wp))
    call write_row(diag_units(26), real(casdiags % piacw(1,1,:), wp))
    call write_row(diag_units(27), real(casdiags % psacw(1,1,:), wp))
    call write_row(diag_units(28), real(casdiags % psacr(1,1,:), wp))
    call write_row(diag_units(29), real(casdiags % pisub(1,1,:), wp))
    call write_row(diag_units(30), real(casdiags % pssub(1,1,:), wp))
    call write_row(diag_units(31), real(casdiags % pimlt(1,1,:), wp))
    call write_row(diag_units(32), real(casdiags % psmlt(1,1,:), wp))
    call write_row(diag_units(33), real(casdiags % psaut(1,1,:), wp))
    call write_row(diag_units(34), real(casdiags % psaci(1,1,:), wp))
    call write_row(diag_units(35), real(casdiags % praut(1,1,:), wp))
    call write_row(diag_units(36), real(casdiags % pracw(1,1,:), wp))
    call write_row(diag_units(37), real(casdiags % prevp(1,1,:), wp))
    call write_row(diag_units(38), real(casdiags % pgacw(1,1,:), wp))
    call write_row(diag_units(39), real(casdiags % pgacs(1,1,:), wp))
    call write_row(diag_units(40), real(casdiags % pgmlt(1,1,:), wp))
    call write_row(diag_units(41), real(casdiags % pgsub(1,1,:), wp))
    call write_row(diag_units(42), real(casdiags % psedi(1,1,:), wp))
    call write_row(diag_units(43), real(casdiags % pseds(1,1,:), wp))
    call write_row(diag_units(44), real(casdiags % psedr(1,1,:), wp))
    call write_row(diag_units(45), real(casdiags % psedg(1,1,:), wp))
    call write_row(diag_units(46), real(casdiags % psedl(1,1,:), wp))
    call write_row(diag_units(47), real(casdiags % pcond(1,1,:), wp))
    call write_row(diag_units(48), real(casdiags % phomr(1,1,:), wp))
    call write_row(diag_units(49), real(casdiags % nhomc(1,1,:), wp))
    call write_row(diag_units(50), real(casdiags % nhomr(1,1,:), wp))
    call write_row(diag_units(51), real(casdiags % nihal(1,1,:), wp))
    call write_row(diag_units(52), real(casdiags % ninuc(1,1,:), wp))
    call write_row(diag_units(53), real(casdiags % nsedi(1,1,:), wp))
    call write_row(diag_units(54), real(casdiags % nseds(1,1,:), wp))
    call write_row(diag_units(55), real(casdiags % nsedg(1,1,:), wp))
    call write_row(diag_units(56), real(casdiags % nraut(1,1,:), wp))
    call write_row(diag_units(57), real(casdiags % nsedl(1,1,:), wp))
    call write_row(diag_units(58), real(casdiags % nracw(1,1,:), wp))
    call write_row(diag_units(59), real(casdiags % nracr(1,1,:), wp))
    call write_row(diag_units(60), real(casdiags % nsedr(1,1,:), wp))
    call write_row(diag_units(61), real(casdiags % nrevp(1,1,:), wp))
    call write_row(diag_units(62), real(casdiags % nisub(1,1,:), wp))
    call write_row(diag_units(63), real(casdiags % nssub(1,1,:), wp))
    call write_row(diag_units(64), real(casdiags % nsaut(1,1,:), wp))
    call write_row(diag_units(65), real(casdiags % nsaci(1,1,:), wp))
    call write_row(diag_units(66), real(casdiags % ngacs(1,1,:), wp))
    call write_row(diag_units(67), real(casdiags % ngsub(1,1,:), wp))
    call write_row(diag_units(68), real(casdiags % niacw(1,1,:), wp))
    call write_row(diag_units(69), real(casdiags % nsacw(1,1,:), wp))
    call write_row(diag_units(70), real(casdiags % nsacr(1,1,:), wp))
    call write_row(diag_units(71), real(casdiags % nimlt(1,1,:), wp))
    call write_row(diag_units(72), real(casdiags % nsmlt(1,1,:), wp))
    call write_row(diag_units(73), real(casdiags % ngacw(1,1,:), wp))
    call write_row(diag_units(74), real(casdiags % ngmlt(1,1,:), wp))
    call write_row(diag_units(75), real(casdiags % pihal(1,1,:), wp))
    call write_row(diag_units(76), real(casdiags % praci_g(1,1,:), wp))
    call write_row(diag_units(77), real(casdiags % praci_r(1,1,:), wp))
    call write_row(diag_units(78), real(casdiags % praci_i(1,1,:), wp))
    call write_row(diag_units(79), real(casdiags % nraci_g(1,1,:), wp))
    call write_row(diag_units(80), real(casdiags % nraci_r(1,1,:), wp))
    call write_row(diag_units(81), real(casdiags % nraci_i(1,1,:), wp))
    call write_row(diag_units(82), real(casdiags % pidps(1,1,:), wp))
    call write_row(diag_units(83), real(casdiags % nidps(1,1,:), wp))
    call write_row(diag_units(84), real(casdiags % pgaci(1,1,:), wp))
    call write_row(diag_units(85), real(casdiags % ngaci(1,1,:), wp))
    call write_row(diag_units(86), real(casdiags % niics_s(1,1,:), wp))
    call write_row(diag_units(87), real(casdiags % niics_i(1,1,:), wp))
  end subroutine write_diagnostics_step

  subroutine write_scalars_initial()
    integer :: m
    do m = 1, nvars_scalar
      write(scalar_units(m), '(ES16.7)') 0.0_wp
    end do
  end subroutine write_scalars_initial

  subroutine write_scalars_step()
    write(scalar_units(1), '(ES16.7)') real(casdiags % SurfaceRainR(1,1), wp)
    write(scalar_units(2), '(ES16.7)') real(casdiags % SurfaceSnowR(1,1), wp)
    write(scalar_units(3), '(ES16.7)') real(casdiags % SurfaceGraupR(1,1), wp)
    write(scalar_units(4), '(ES16.7)') real(casdiags % lwp(1,1), wp)
    write(scalar_units(5), '(ES16.7)') real(casdiags % rwp(1,1), wp)
    write(scalar_units(6), '(ES16.7)') real(casdiags % iwp(1,1), wp)
    write(scalar_units(7), '(ES16.7)') real(casdiags % swp(1,1), wp)
    write(scalar_units(8), '(ES16.7)') real(casdiags % gwp(1,1), wp)
    write(scalar_units(9), '(ES16.7)') real(casdiags % SurfaceCloudR(1,1), wp)
  end subroutine write_scalars_step

end program casim_column_model
