#!/usr/bin/env python3
"""
Convert the plain-text output of the CASIM 1D column model into a single
NetCDF file with (time, height) and (time,) dimensioned variables.

Usage:
    python3 make_netcdf.py [output_dir] [output_file]

Defaults:
    output_dir  = 'output'
    output_file = 'output/casim_column_kgo.nc'
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset

# Profile variables: (file stem, netcdf var name, long_name, units)
PROG_VARS = [
    ('qv',      'qv',      'water vapour specific humidity', 'kg kg-1'),
    ('theta',   'theta',   'potential temperature', 'K'),
    ('q1_ql',   'ql',      'cloud liquid water specific humidity', 'kg kg-1'),
    ('q2_qr',   'qr',      'rain water specific humidity', 'kg kg-1'),
    ('q3_nl',   'nl',      'cloud droplet number mixing ratio', 'kg-1'),
    ('q4_nr',   'nr',      'rain drop number mixing ratio', 'kg-1'),
    ('q5_m3r',  'm3r',     'rain 3rd moment mixing ratio', 'm3 kg-1'),
    ('q6_qi',   'qi',      'ice specific humidity', 'kg kg-1'),
    ('q7_qs',   'qs',      'snow specific humidity', 'kg kg-1'),
    ('q8_qg',   'qg',      'graupel specific humidity', 'kg kg-1'),
    ('q9_ni',   'ni',      'ice crystal number mixing ratio', 'kg-1'),
    ('q10_ns',  'ns',      'snow number mixing ratio', 'kg-1'),
    ('q11_ng',  'ng',      'graupel number mixing ratio', 'kg-1'),
    ('q12_m3s', 'm3s',     'snow 3rd moment mixing ratio', 'm3 kg-1'),
    ('q13_m3g', 'm3g',     'graupel 3rd moment mixing ratio', 'm3 kg-1'),
]

DIAG_VARS = [
    ('rainfall_3d',  'rainfall_3d',  'rain sedimentation flux', 'kg m-2 s-1'),
    ('snowfall_3d',  'snowfall_3d',  'total snow+graupel sedimentation flux', 'kg m-2 s-1'),
    ('snowonly_3d',  'snowonly_3d',  'snow only sedimentation flux', 'kg m-2 s-1'),
    ('graupfall_3d', 'graupfall_3d', 'graupel sedimentation flux', 'kg m-2 s-1'),
    ('dth_total',    'dth_total',    'total potential temperature tendency', 'K s-1'),
    ('dqv_total',    'dqv_total',    'total water vapour tendency', 'kg kg-1 s-1'),
    ('dqc',          'dqc',          'cloud water tendency', 'kg kg-1 s-1'),
    ('dqr',          'dqr',          'rain water tendency', 'kg kg-1 s-1'),
    ('dqi',          'dqi',          'ice tendency', 'kg kg-1 s-1'),
    ('dqs',          'dqs',          'snow tendency', 'kg kg-1 s-1'),
    ('dqg',          'dqg',          'graupel tendency', 'kg kg-1 s-1'),
    ('w',            'w',            'vertical velocity', 'm s-1'),
    ('rho',          'rho',          'air density', 'kg m-3'),
    ('pressure',     'pressure',     'air pressure', 'Pa'),
    ('exner',        'exner',        'Exner function', '1'),
    # Radar reflectivity diagnostics (casdiags % l_radar)
    ('dbz_tot', 'dbz_tot', 'total radar reflectivity', 'dBZ'),
    ('dbz_g',   'dbz_g',   'graupel radar reflectivity', 'dBZ'),
    ('dbz_i',   'dbz_i',   'ice radar reflectivity', 'dBZ'),
    ('dbz_s',   'dbz_s',   'snow radar reflectivity', 'dBZ'),
    ('dbz_l',   'dbz_l',   'liquid radar reflectivity', 'dBZ'),
    ('dbz_r',   'dbz_r',   'rain radar reflectivity', 'dBZ'),
    # Process-rate and number-tendency diagnostics (procs structure); see
    # src/generic_diagnostic_variables.F90 for full descriptions.
    ('phomc', 'phomc', 'homogeneous nucleation rate of cloud', 'kg kg-1 s-1'),
    ('pinuc', 'pinuc', 'heterogeneous ice nucleation rate', 'kg kg-1 s-1'),
    ('pidep', 'pidep', 'deposition rate of ice crystals', 'kg kg-1 s-1'),
    ('psdep', 'psdep', 'deposition rate for snow', 'kg kg-1 s-1'),
    ('piacw', 'piacw', 'ice-water accretion (riming) rate', 'kg kg-1 s-1'),
    ('psacw', 'psacw', 'snow-cloud water accretion (riming) rate', 'kg kg-1 s-1'),
    ('psacr', 'psacr', 'snow-rain accretion rate', 'kg kg-1 s-1'),
    ('pisub', 'pisub', 'sublimation rate of ice crystals', 'kg kg-1 s-1'),
    ('pssub', 'pssub', 'sublimation rate of snow', 'kg kg-1 s-1'),
    ('pimlt', 'pimlt', 'melting rate of ice crystals', 'kg kg-1 s-1'),
    ('psmlt', 'psmlt', 'melting rate of snow', 'kg kg-1 s-1'),
    ('psaut', 'psaut', 'snow autoconversion rate (from ice)', 'kg kg-1 s-1'),
    ('psaci', 'psaci', 'snow-ice accretion rate', 'kg kg-1 s-1'),
    ('praut', 'praut', 'rain autoconversion rate', 'kg kg-1 s-1'),
    ('pracw', 'pracw', 'rain-cloud water accretion rate', 'kg kg-1 s-1'),
    ('prevp', 'prevp', 'rain evaporation rate', 'kg kg-1 s-1'),
    ('pgacw', 'pgacw', 'graupel-cloud water accretion rate', 'kg kg-1 s-1'),
    ('pgacs', 'pgacs', 'graupel-snow accretion rate', 'kg kg-1 s-1'),
    ('pgmlt', 'pgmlt', 'graupel melting rate', 'kg kg-1 s-1'),
    ('pgsub', 'pgsub', 'graupel sublimation rate', 'kg kg-1 s-1'),
    ('psedi', 'psedi', 'ice crystal sedimentation rate', 'kg kg-1 s-1'),
    ('pseds', 'pseds', 'snow sedimentation rate', 'kg kg-1 s-1'),
    ('psedr', 'psedr', 'rain sedimentation rate', 'kg kg-1 s-1'),
    ('psedg', 'psedg', 'graupel sedimentation rate', 'kg kg-1 s-1'),
    ('psedl', 'psedl', 'cloud liquid sedimentation rate', 'kg kg-1 s-1'),
    ('pcond', 'pcond', 'condensation/evaporation rate', 'kg kg-1 s-1'),
    ('phomr', 'phomr', 'homogeneous freezing rate of rain', 'kg kg-1 s-1'),
    ('nhomc', 'nhomc', 'homogeneous nucleation number tendency (cloud)', 'kg-1 s-1'),
    ('nhomr', 'nhomr', 'homogeneous freezing number tendency (rain)', 'kg-1 s-1'),
    ('nihal', 'nihal', 'Hallett-Mossop secondary ice number tendency', 'kg-1 s-1'),
    ('ninuc', 'ninuc', 'heterogeneous ice nucleation number tendency', 'kg-1 s-1'),
    ('nsedi', 'nsedi', 'ice crystal sedimentation number tendency', 'kg-1 s-1'),
    ('nseds', 'nseds', 'snow sedimentation number tendency', 'kg-1 s-1'),
    ('nsedg', 'nsedg', 'graupel sedimentation number tendency', 'kg-1 s-1'),
    ('nraut', 'nraut', 'rain autoconversion number tendency', 'kg-1 s-1'),
    ('nsedl', 'nsedl', 'cloud liquid sedimentation number tendency', 'kg-1 s-1'),
    ('nracw', 'nracw', 'rain-cloud water accretion number tendency', 'kg-1 s-1'),
    ('nracr', 'nracr', 'rain-rain accretion number tendency', 'kg-1 s-1'),
    ('nsedr', 'nsedr', 'rain sedimentation number tendency', 'kg-1 s-1'),
    ('nrevp', 'nrevp', 'rain evaporation number tendency', 'kg-1 s-1'),
    ('nisub', 'nisub', 'ice sublimation number tendency', 'kg-1 s-1'),
    ('nssub', 'nssub', 'snow sublimation number tendency', 'kg-1 s-1'),
    ('nsaut', 'nsaut', 'snow autoconversion number tendency', 'kg-1 s-1'),
    ('nsaci', 'nsaci', 'snow-ice accretion number tendency', 'kg-1 s-1'),
    ('ngacs', 'ngacs', 'graupel-snow accretion number tendency', 'kg-1 s-1'),
    ('ngsub', 'ngsub', 'graupel sublimation number tendency', 'kg-1 s-1'),
    ('niacw', 'niacw', 'ice-water accretion number tendency', 'kg-1 s-1'),
    ('nsacw', 'nsacw', 'snow-cloud water accretion number tendency', 'kg-1 s-1'),
    ('nsacr', 'nsacr', 'snow-rain accretion number tendency', 'kg-1 s-1'),
    ('nimlt', 'nimlt', 'ice melting number tendency', 'kg-1 s-1'),
    ('nsmlt', 'nsmlt', 'snow melting number tendency', 'kg-1 s-1'),
    ('ngacw', 'ngacw', 'graupel-cloud water accretion number tendency', 'kg-1 s-1'),
    ('ngmlt', 'ngmlt', 'graupel melting number tendency', 'kg-1 s-1'),
    ('pihal', 'pihal', 'Hallett-Mossop secondary ice production rate', 'kg kg-1 s-1'),
    ('praci_g', 'praci_g', 'rain-ice accretion rate onto graupel', 'kg kg-1 s-1'),
    ('praci_r', 'praci_r', 'rain-ice accretion rate onto rain', 'kg kg-1 s-1'),
    ('praci_i', 'praci_i', 'rain-ice accretion rate onto ice', 'kg kg-1 s-1'),
    ('nraci_g', 'nraci_g', 'rain-ice accretion number tendency (graupel)', 'kg-1 s-1'),
    ('nraci_r', 'nraci_r', 'rain-ice accretion number tendency (rain)', 'kg-1 s-1'),
    ('nraci_i', 'nraci_i', 'rain-ice accretion number tendency (ice)', 'kg-1 s-1'),
    ('pidps', 'pidps', 'droplet shattering secondary ice production rate', 'kg kg-1 s-1'),
    ('nidps', 'nidps', 'droplet shattering secondary ice number tendency', 'kg-1 s-1'),
    ('pgaci', 'pgaci', 'graupel-ice accretion rate', 'kg kg-1 s-1'),
    ('ngaci', 'ngaci', 'graupel-ice accretion number tendency', 'kg-1 s-1'),
    ('niics_s', 'niics_s', 'ice collision secondary ice number tendency (snow)', 'kg-1 s-1'),
    ('niics_i', 'niics_i', 'ice collision secondary ice number tendency (ice)', 'kg-1 s-1'),
]

SCALAR_VARS = [
    ('surface_rain',  'surface_rain_rate',  'surface rainfall rate', 'kg m-2 s-1'),
    ('surface_snow',  'surface_snow_rate',  'surface snowfall rate', 'kg m-2 s-1'),
    ('surface_graup', 'surface_graup_rate', 'surface graupel fall rate', 'kg m-2 s-1'),
    ('lwp',           'lwp',                'liquid water path', 'kg m-2'),
    ('rwp',           'rwp',                'rain water path', 'kg m-2'),
    ('iwp',           'iwp',                'ice water path', 'kg m-2'),
    ('swp',           'swp',                'snow water path', 'kg m-2'),
    ('gwp',           'gwp',                'graupel water path', 'kg m-2'),
    ('surface_cloud', 'surface_cloud_rate', 'surface cloud liquid deposition rate', 'kg m-2 s-1'),
]


def load_column(path):
    """Load a whitespace-delimited text file of shape (ntime, nz)."""
    return np.loadtxt(path, ndmin=2)


def load_scalar(path):
    """Load a text file of one value per line, shape (ntime,)."""
    return np.loadtxt(path, ndmin=1)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else 'output'
    outfile = sys.argv[2] if len(sys.argv) > 2 else os.path.join(outdir, 'casim_column_kgo.nc')

    height = np.loadtxt(os.path.join(outdir, 'height.txt'), ndmin=1)
    time = np.loadtxt(os.path.join(outdir, 'time.txt'), ndmin=1)
    nz = height.shape[0]
    nt = time.shape[0]

    with Dataset(outfile, 'w', format='NETCDF4') as nc:
        nc.createDimension('time', nt)
        nc.createDimension('height', nz)

        nc.title = 'CASIM 1D column model output'
        nc.source = 'casim_column_model (shipway_microphysics driver)'

        v_height = nc.createVariable('height', 'f8', ('height',))
        v_height.units = 'm'
        v_height.long_name = 'height above surface'
        v_height[:] = height

        v_time = nc.createVariable('time', 'f8', ('time',))
        v_time.units = 's'
        v_time.long_name = 'time since start of simulation'
        v_time[:] = time

        for stem, name, long_name, units in PROG_VARS + DIAG_VARS:
            path = os.path.join(outdir, stem + '.txt')
            if not os.path.exists(path):
                print(f'Warning: missing {path}, skipping')
                continue
            data = load_column(path)
            var = nc.createVariable(name, 'f8', ('time', 'height'), zlib=True)
            var.units = units
            var.long_name = long_name
            var[:, :] = data

        for stem, name, long_name, units in SCALAR_VARS:
            path = os.path.join(outdir, stem + '.txt')
            if not os.path.exists(path):
                print(f'Warning: missing {path}, skipping')
                continue
            data = load_scalar(path)
            var = nc.createVariable(name, 'f8', ('time',), zlib=True)
            var.units = units
            var.long_name = long_name
            var[:] = data

    print(f'Wrote {outfile}')


if __name__ == '__main__':
    main()
