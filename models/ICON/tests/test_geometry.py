#!/usr/bin/env python3
"""Integration tests; requires ncgen and the compiled geometry_check executable.

Uses temporary synthetic NetCDF inputs, never edits scientific data. Pass a
--work-dir under the experiment directory to retain fixtures and individual logs.
"""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

GRID = '''netcdf grid {
dimensions: cell=4; vertex=3; nv=3;
variables:
double clon(cell); clon:units="radian";
double clat(cell); clat:units="radian";
double vlon(vertex); vlon:units="radian";
double vlat(vertex); vlat:units="radian";
int vertex_of_cell(nv,cell); int neighbor_cell_index(nv,cell);
int refin_c_ctrl(cell);
:uuidOfHGrid="test-horizontal"; :global_grid=0; :semi_major_axis=6371229.;
data:
clon=0,0.1,0.2,0.3; clat=0,0.1,0.2,0.3;
vlon=0,0.2,0.3; vlat=0,0.2,0.3;
vertex_of_cell=1,1,1,1,2,2,2,2,3,3,3,3;
neighbor_cell_index=-1,1,2,3,2,3,4,1,3,4,1,2;
refin_c_ctrl=1,4,5,0;
}'''
VERTICAL = '''netcdf geometry {
dimensions: ncells=4; height=2; height_2=3;
variables:
double z_mc(height,ncells); z_mc:units="m";
double z_ifc(height_2,ncells); z_ifc:units="m";
:uuidOfHGrid="test-horizontal"; :uuidOfVGrid="test-vertical";
data:
z_mc=150,160,170,180,50,60,70,80;
z_ifc=200,210,220,230,100,110,120,130,0,10,20,30;
}'''


def run_tests(executable, root):
    # Positive cases distinguish the specified zone, nudging zone, and global case.
    cases = [
        ('lam', GRID, VERTICAL, 1, 'excluded cells=3'),
        ('no_nudging', GRID, VERTICAL, 0, 'excluded cells=2'),
        ('global', GRID.replace(':global_grid=0', ':global_grid=1'), VERTICAL, -1, 'excluded cells=0'),
        ('float_heights', GRID, VERTICAL.replace('double z_', 'float z_'), 1, 'excluded cells=3'),
        ('axis_names', GRID, VERTICAL.replace('height_2', 'interfaces').replace('height', 'full'), 1,
         'excluded cells=3'),
    ]
    failures = [
        ('missing_width', GRID, VERTICAL, -1, 'LAM requires explicit'),
        ('width_coverage', GRID, VERTICAL, 2, 'Insufficient refin_c_ctrl coverage'),
        ('uuid', GRID, VERTICAL.replace('test-horizontal', 'other-grid'), 1, 'Horizontal UUID mismatch'),
        ('units', GRID, VERTICAL.replace('units="m"', 'units="km"'), 1, 'Expected metres'),
        ('latitude', GRID.replace('clat=0,0.1', 'clat=2,0.1'), VERTICAL, 1, 'Latitude outside'),
        ('neighbor', GRID.replace('neighbor_cell_index=-1', 'neighbor_cell_index=99'), VERTICAL, 1,
         'Invalid neighbor_cell_index'),
        ('vertex', GRID.replace('vertex_of_cell=1', 'vertex_of_cell=0'), VERTICAL, 1, 'Invalid vertex_of_cell'),
        ('bracketing', GRID, VERTICAL.replace('z_mc=150', 'z_mc=250'), 1, 'strictly between'),
        ('fill', GRID, VERTICAL.replace('z_mc:units="m";', 'z_mc:units="m"; z_mc:_FillValue=150.;'), 1,
         'Missing values'),
        ('nan', GRID, VERTICAL.replace('z_mc=150', 'z_mc=NaN'), 1, 'Nonfinite values'),
        ('cell_count', GRID, VERTICAL.replace('ncells=4', 'ncells=5'), 1, 'Cell count/order mismatch'),
        ('axis_order', GRID, VERTICAL.replace('z_mc(height,ncells)', 'z_mc(ncells,height)'), 1,
         'Cell count/order mismatch'),
    ]
    for name, grid, vertical, width, expected in cases + failures:
        directory = root / name
        directory.mkdir()
        for stem, content in [('grid', grid), ('vertical', vertical)]:
            cdl = directory / (stem + '.cdl')
            cdl.write_text(content)
            subprocess.run(['ncgen', '-o', str(directory / (stem + '.nc')), str(cdl)], check=True)
        (directory / 'input.nml').write_text(
            '&utilities_nml\n/\n&icon_geometry_nml\n'
            " horizontal_grid_file='grid.nc', vertical_geometry_file='vertical.nc',\n"
            f' nudge_zone_width={width}\n/\n')
        result = subprocess.run([str(executable)], cwd=directory, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
        (directory / 'test.log').write_text(result.stdout)
        # Some Fortran STOP/error-handler implementations return zero on fatal errors.
        passed = 'geometry_check PASS' in result.stdout
        assert expected in result.stdout, (name, result.stdout)
        if name in {c[0] for c in cases}:
            assert result.returncode == 0 and passed, (name, result.stdout)
        else:
            assert not passed, (name, result.stdout)
        print('PASS', name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--work-dir', type=Path)
    args = parser.parse_args()
    assert shutil.which('ncgen'), 'Load the NetCDF environment first'
    if args.work_dir:
        args.work_dir.mkdir(parents=True, exist_ok=False)
        run_tests(args.executable.resolve(), args.work_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='icon-geometry-') as temporary:
            run_tests(args.executable.resolve(), Path(temporary))
