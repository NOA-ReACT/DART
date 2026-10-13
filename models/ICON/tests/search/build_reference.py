#!/usr/bin/env python3
"""Build pinned ICON search oracle, independently of all adapted DART modules.

Numerical routines are extracted from verified source without algebraic changes.
Two diagnostic assignments expose the selected triangle and scaling iteration.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

ICON = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2]).resolve()
base = Path(__file__).resolve().parents[1] / 'triangulation/build_reference.py'
subprocess.run([sys.executable, str(base), str(ICON), str(OUT)], check=True)
support = ICON / 'externals/fortran-support'
assert subprocess.check_output(['git','-C',str(support),'rev-parse','HEAD'],text=True).strip() == \
    '34246610b17db29f214fb2f95ca1c9087f09b89c'
hashes = {}
for source in ['src/mo_octree.F90', 'src/mo_util_sort.F90']:
    content = (support/source).read_bytes()
    assert content == subprocess.check_output(['git','-C',str(support),'show','HEAD:'+source])
    (OUT/Path(source).name).write_bytes(content)
    hashes[source] = hashlib.sha256(content).hexdigest()

def routine(text,name,kind='SUBROUTINE'):
    return re.search(r'^ *[^!\n]*\b'+kind+' '+name+r'\b.*?^ *END '+kind+' '+name+r'\b[^\n]*',
                     text,re.M|re.S).group(0)

math=(OUT/'mo_math_utilities.F90').read_text()
coords=(OUT/'reference_coordinates.f90').read_text()
funcs='\n'.join(routine(math,n,'FUNCTION') for n in ['project_point_to_plane','cc_dot_product','vector_product'])
# Expand exactly the source macros; keep the native projection evaluation order.
funcs=funcs.replace('d_normalize(norm_to_plane)',
                  'norm_to_plane%x=norm_to_plane%x/SQRT(DOT_PRODUCT(norm_to_plane%x,norm_to_plane%x))')
(OUT/'reference_coordinates.f90').write_text(coords.replace('end module\n',funcs+'\nend module\n'))
(OUT/'support_kind.f90').write_text('''module mo_iconlib_kind
integer,parameter :: wp=selected_real_kind(13),dp=selected_real_kind(13),sp=kind(1.)
end module
''')
src=(OUT/'mo_intp_lonlat_baryctr.f90').read_text()
funcs='\n'.join([routine(src,'compute_barycentric_coords'),routine(src,'inside_triangle','FUNCTION'),
                routine(src,'inside_scaled_triangle','FUNCTION'),routine(src,'compute_triangle_bboxes'),
                routine(src,'compute_barycentric_coordinates')])
funcs=funcs.replace('            IF (dbg_level > 5) THEN',
                    '            reference_retry=i_scale\n            reference_triangle=idx0+1\n'
                    '            IF (dbg_level > 5) THEN')
normalize_source=ICON/'src/shr_horizontal/mo_intp_lonlat.f90'
raw=normalize_source.read_bytes()
assert raw==subprocess.check_output(['git','-C',str(ICON),'show','HEAD:src/shr_horizontal/mo_intp_lonlat.f90'])
hashes['src/shr_horizontal/mo_intp_lonlat.f90']=hashlib.sha256(raw).hexdigest()
funcs+='\n'+routine(raw.decode(),'normalize_weights')
(OUT/'reference_search.F90').write_text('''! ICON BSD-3-Clause numerical routines with test-only infrastructure.
module reference_search
use mo_kind,only:wp
use mo_exception,only:finish
use mo_delaunay_types
use mo_octree
use reference_coordinates, gc2cc=>gc2cc_sphere
implicit none
character(len=*),parameter :: modname='reference_search'
integer,parameter :: dbg_level=0,nproma=1,SUCCESS=0
integer :: reference_retry=-1,reference_triangle=0
type t_intp_scalar_coeff
real(wp),allocatable :: coeff(:,:,:),v(:,:,:,:)
integer,allocatable :: stencil(:,:),idx(:,:,:),blk(:,:,:)
end type
type t_lon_lat_intp
type(t_intp_scalar_coeff) :: baryctr
type(t_geographical_coordinates),allocatable :: ll_coord(:,:)
contains
procedure :: nblks_lonlat
procedure :: npromz_lonlat
end type
contains
integer function nblks_lonlat(this,n)
class(t_lon_lat_intp),intent(in)::this
integer,intent(in)::n
nblks_lonlat=size(this%ll_coord,2)
end function
integer function npromz_lonlat(this,n)
class(t_lon_lat_intp),intent(in)::this
integer,intent(in)::n
npromz_lonlat=1
end function
elemental integer function idx_no(i)
integer,intent(in)::i
idx_no=i
end function
elemental integer function blk_no(i)
integer,intent(in)::i
blk_no=1
end function
integer function get_my_mpi_work_id()
get_my_mpi_work_id=0
end function
function int2string(i) result(s)
integer,intent(in)::i
character(len=32)::s
write(s,'(I0)') i
end function
subroutine message(a,b)
character(len=*),intent(in)::a,b
write(*,*) a,b
end subroutine
'''+funcs+'\nend module\n')
(OUT/'search_driver.f90').write_text(Path(__file__).with_name('reference_driver.f90.in').read_text())
flags=['-O0','-g','-traceback','-check','bounds','-fp-model','precise','-fpp','-DNOMPI','-D__HAVE_QUAD_PRECISION']
inc=subprocess.check_output(['nf-config','--fflags'],text=True).split()
libs=(subprocess.check_output(['nf-config','--flibs'],text=True).split()+
      subprocess.check_output(['nc-config','--libs'],text=True).split())
cmd=['ifx',*flags,*inc,'adapters.f90','support_kind.f90','reference_coordinates.f90','mo_delaunay_types.f90',
     'mo_util_sort.F90','mo_octree.F90','reference_search.F90','search_driver.f90',*libs,'-o','icon_search_reference']
(OUT/'search_provenance.json').write_text(json.dumps(dict(hashes=hashes,command=cmd),indent=2)+'\n')
with (OUT/'search_build.log').open('w') as log:
    subprocess.run(cmd,cwd=OUT,stdout=log,stderr=subprocess.STDOUT,check=True)
