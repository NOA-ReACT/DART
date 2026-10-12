#!/usr/bin/env python3
"""Test DART indexing and ICON locations using small sloping-column fixtures."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import sys
sys.dont_write_bytecode = True
from test_geometry import GRID, VERTICAL

FIELDS = [
    ('u', 'QTY_U_WIND_COMPONENT', 'm s-1'),
    ('v', 'QTY_V_WIND_COMPONENT', 'm s-1'),
    ('w', 'QTY_VERTICAL_VELOCITY', 'm s-1'),
    ('temp', 'QTY_TEMPERATURE', 'K'),
    ('qv', 'QTY_SPECIFIC_HUMIDITY', 'kg kg-1'),
    ('qc', 'QTY_CLOUDWATER_MIXING_RATIO', 'kg kg-1'),
    ('qi', 'QTY_ICE_MIXING_RATIO', 'kg kg-1'),
    ('qr', 'QTY_RAINWATER_MIXING_RATIO', 'kg kg-1'),
    ('qs', 'QTY_SNOW_MIXING_RATIO', 'kg kg-1'),
    ('qg', 'QTY_GRAUPEL_MIXING_RATIO', 'kg kg-1'),
    ('qh', 'QTY_HAIL_MIXING_RATIO', 'kg kg-1'),
    ('qnc', 'QTY_DROPLET_NUMBER_CONCENTR', 'kg-1'),
    ('qni', 'QTY_ICE_NUMBER_CONCENTRATION', 'kg-1'),
    ('qnr', 'QTY_RAIN_NUMBER_CONCENTR', 'kg-1'),
    ('qns', 'QTY_SNOW_NUMBER_CONCENTR', 'kg-1'),
    ('qng', 'QTY_GRAUPEL_NUMBER_CONCENTR', 'kg-1'),
    ('qnh', 'QTY_HAIL_NUMBER_CONCENTR', 'kg-1'),
]

def history(records=1):
    declarations, data = [], ['time=' + ','.join(str(60*r) for r in range(records)) + ';']
    for name, qty, unit in FIELDS:
        axis = 'height' if name == 'w' else 'height_2'
        declarations.append(f'float {name}(time,{axis},ncells); {name}:units=" {unit}";')
        data.append(name + '=' + ','.join(['1'] * ((12 if name == 'w' else 8) * records)) + ';')
    return '\n'.join(['netcdf history {',
        'dimensions: time=UNLIMITED; ncells=4; height=3; height_2=2;',
        'variables: double time(time); time:units="minutes since 2024-11-26 00:00:00";',
        'time:calendar="gregorian";', *declarations,
        ':uuidOfHGrid="test-horizontal"; :uuidOfVGrid="test-vertical";',
        'data:', *data, '}'])


def table(fields):
    return 'state_variables=' + ',\n'.join(f"'{name}','{qty}','UPDATE'" for name, qty, unit in fields)


def run_tests(executable, root):
    source = history()
    cases = [
        ('default', source, '', None, 60),
        ('temperature_only', source, table([FIELDS[3]]), None, 8),
        ('w_only', source, table([FIELDS[2]]), None, 12),
        ('reordered', source, table([FIELDS[2], FIELDS[0]]), None, 20),
        ('all_hydrometeors', source, table(FIELDS), None, 140),
        ('no_copy_back', source, table([FIELDS[0]]).replace('UPDATE', 'NO_COPY_BACK'), None, 8),
        ('pressure', source, "state_variables='pres','QTY_PRESSURE','UPDATE'", 'pressure deferred', None),
        ('duplicate', source, table([FIELDS[0], FIELDS[0]]), 'Duplicate state variable', None),
        ('wrong_quantity', source, table([FIELDS[0]]).replace('QTY_U_WIND_COMPONENT','QTY_TEMPERATURE'),
         'Incorrect quantity', None),
        ('incomplete', source, "state_variables='u','QTY_U_WIND_COMPONENT'", 'complete name/quantity', None),
        ('missing_field', source.replace('float u(', 'float absent(').replace('u:units','absent:units')
         .replace('\nu=', '\nabsent='), '', 'NetCDF: Variable not found', None),
        ('wrong_stagger', source.replace('w(time,height,ncells)', 'w(time,height_2,ncells)')
         .replace('w=1,1,1,1,1,1,1,1,1,1,1,1;', 'w=1,1,1,1,1,1,1,1;'), '', 'Incorrect vertical staggering', None),
        ('wrong_units', source.replace('temp:units=" K"','temp:units=" degC"'), '', 'Incorrect units', None),
        ('horizontal_uuid', source.replace('test-horizontal','other-horizontal'), '', 'horizontal UUID mismatch', None),
        ('vertical_uuid', source.replace('test-vertical','other-vertical'), '', 'vertical UUID mismatch', None),
        ('two_times', history(2), '', 'one unlimited time record', None),
        ('zero_index', source, '', 'State index out of range', None),
        ('large_index', source, '', 'State index out of range', None),
    ]
    for name, state, selection, error, size in cases:
        directory = root / name
        directory.mkdir()
        for stem, content in [('grid', GRID), ('geometry', VERTICAL), ('history', state)]:
            cdl = directory / (stem + '.cdl')
            cdl.write_text(content)
            subprocess.run(['ncgen', '-o', str(directory / (stem + '.nc')), str(cdl)], check=True)
        (directory / 'input.nml').write_text(
            "&utilities_nml\n/\n&obs_kind_nml\n/\n&location_nml\n/\n&icon_geometry_nml\n horizontal_grid_file='grid.nc',\n"
            " vertical_geometry_file='geometry.nc', nudge_zone_width=1\n/\n"
            "&model_nml\n template_file='history.nc',\n" + selection + '\n/\n')
        args = [str(executable)]
        if name == 'zero_index': args.append('0')
        if name == 'large_index': args.append('61')
        result = subprocess.run(args, cwd=directory, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=60)
        (directory / 'test.log').write_text(result.stdout)
        if error:
            assert error in result.stdout and 'state_metadata_check PASS' not in result.stdout, (name, result.stdout)
        else:
            assert result.returncode == 0 and 'state_metadata_check PASS' in result.stdout, (name, result.stdout)
            assert f'STATE SIZE {size}\n' in result.stdout, (name, result.stdout)
        print('PASS', name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--work-dir', type=Path)
    args = parser.parse_args()
    if args.work_dir:
        args.work_dir.mkdir(parents=True, exist_ok=False)
        run_tests(args.executable.resolve(), args.work_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='icon-state-') as temporary:
            run_tests(args.executable.resolve(), Path(temporary))
