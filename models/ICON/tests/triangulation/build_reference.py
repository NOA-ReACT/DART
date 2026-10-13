#!/usr/bin/env python3
"""Build the pinned ICON serial oracle; no adapted DART numerical modules are used.

Run in a compute allocation with the Intel/NetCDF environment loaded. The two
Delaunay sources are copied verbatim. Only infrastructure is supplied locally;
the coordinate routines are extracted verbatim from pinned iconmath.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

ICON_COMMIT = '8597da45ef4b86323f3fb844caedc4ae5e1ffc01'
MATH_COMMIT = '622c01ba78e4fbb19fdf5f3c9f18a90b1e551389'


def routine(text, name):
    match = re.search(r'^  ELEMENTAL FUNCTION ' + name + r'\b.*?^  END FUNCTION ' + name + r'\b',
                      text, re.M | re.S)
    assert match, name
    return match[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('icon', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    icon = args.icon.resolve()
    root = args.output.resolve()
    for path, revision in [(icon, ICON_COMMIT), (icon/'externals/iconmath', MATH_COMMIT)]:
        actual = subprocess.check_output(['git','-C',str(path),'rev-parse','HEAD'],text=True).strip()
        assert actual == revision, (path, actual)
    root.mkdir(parents=True, exist_ok=False)
    hashes = {}
    sources = ['src/grid/mo_delaunay_types.f90', 'src/grid/mo_delaunay.f90',
               'externals/iconmath/src/support/mo_math_utilities.F90',
               'src/shr_horizontal/mo_intp_lonlat_baryctr.f90']
    for source in sources:
        content = (icon/source).read_bytes()
        # Verify the reference source itself has no local edits.
        repo = icon/'externals/iconmath' if source.startswith('externals/iconmath/') else icon
        relative = source.removeprefix('externals/iconmath/') if repo != icon else source
        original = subprocess.check_output(['git','-C',str(repo),'show',f'HEAD:{relative}'])
        assert content == original, f'Modified reference: {source}'
        hashes[source] = hashlib.sha256(content).hexdigest()
        (root/Path(source).name).write_bytes(content)
    (root/'adapters.f90').write_text('''module mo_kind
integer, parameter :: wp=selected_real_kind(13), qp=selected_real_kind(33,4931)
end module
module mo_exception
contains
subroutine finish(routine,message)
character(len=*), intent(in) :: routine,message
write(*,*) trim(routine), ': ', trim(message)
error stop 1
end subroutine
end module
module mo_impl_constants
integer, parameter :: SUCCESS=0
end module
module mo_mpi
integer, parameter :: p_comm_work=0, p_real=0
end module
module mo_util_file
contains
logical function util_file_is_writable(filename)
character(len=*), intent(in) :: filename
! Reference never uses ICON's NetCDF writer; do not allow accidental use.
error stop 'Reference file-writing adapter is not implemented'
util_file_is_writable=.false.
end function
end module
module mo_netcdf
use netcdf
end module
module mo_netcdf_errhandler
contains
subroutine nf(status,routine)
integer, intent(in) :: status
character(len=*), intent(in) :: routine
if (status/=0) then
write(*,*) routine,status
error stop 1
endif
end subroutine
end module
''')
    math = (icon/sources[2]).read_text()
    (root/'reference_coordinates.f90').write_text('''! Verbatim ICON BSD-3-Clause coordinate routines; see reference provenance.
module reference_coordinates
use mo_kind, only : wp
implicit none
real(wp), parameter :: dbl_eps=epsilon(1._wp)
real(wp), parameter :: pi_2=1.57079632679489661923132169163975144_wp
type t_cartesian_coordinates
real(wp) :: x(3)
end type
type t_geographical_coordinates
real(wp) :: lon,lat
end type
contains
''' + routine(math,'cc2gc_sphere')+'\n'+routine(math,'gc2cc_sphere')+'\nend module\n')
    driver = Path(__file__).with_name('reference_driver.f90.in')
    shutil.copyfile(driver, root/'reference_driver.f90')
    flags = ['-O0','-g','-traceback','-check','bounds','-fp-model','precise',
             '-fpp','-DNOMPI','-D__HAVE_QUAD_PRECISION']
    includes = subprocess.check_output(['nf-config','--fflags'],text=True).split()
    libs = (subprocess.check_output(['nf-config','--flibs'],text=True).split() +
            subprocess.check_output(['nc-config','--libs'],text=True).split())
    command = ['ifx',*flags,*includes,'adapters.f90','reference_coordinates.f90',
               'mo_delaunay_types.f90','mo_delaunay.f90','reference_driver.f90',*libs,'-o','icon_reference']
    (root/'provenance.json').write_text(json.dumps(dict(icon=ICON_COMMIT,iconmath=MATH_COMMIT,
                                        hashes=hashes,command=command),indent=2)+'\n')
    with (root/'build.log').open('w') as log:
        subprocess.run(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True)

if __name__ == '__main__':
    main()
