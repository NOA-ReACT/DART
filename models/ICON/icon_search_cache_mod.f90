! DART software - Copyright UCAR. See the DART license for terms of use.
module icon_search_cache_mod
! Portable numeric NetCDF representation, not a compiler-dependent type dump.
! Thin structural/identity validation only; no production Delaunay validation.
use types_mod, only : r8
use utilities_mod, only : error_handler, E_ERR, E_WARN, E_MSG
use netcdf_utilities_mod, only : nc_check
use netcdf
use icon_horizontal_search_mod, only : horizontal_search_type, build_horizontal_search
use icon_octree_mod, only : OCTREE_DEPTH
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use, intrinsic :: iso_c_binding, only : c_int, c_char, c_null_char
implicit none
private
public :: load_or_build_search, read_search_cache, write_search_cache
integer, parameter :: schema_version=1, algorithm_version=1
character(len=*), parameter :: routine='icon_search_cache_mod'
interface
   function c_mkstemp(path) bind(C,name='mkstemp') result(fd)
      import c_int,c_char
      character(c_char), intent(inout) :: path(*)
      integer(c_int) :: fd
   end function
   function c_close(fd) bind(C,name='close') result(status)
      import c_int
      integer(c_int), value :: fd
      integer(c_int) :: status
   end function
   function c_link(old,new) bind(C,name='link') result(status)
      import c_int,c_char
      character(c_char), intent(in) :: old(*),new(*)
      integer(c_int) :: status
   end function
   function c_unlink(path) bind(C,name='unlink') result(status)
      import c_int,c_char
      character(c_char), intent(in) :: path(*)
      integer(c_int) :: status
   end function
end interface
contains
subroutine require(ok,message)
logical, intent(in) :: ok
character(len=*), intent(in) :: message
if (.not. ok) call error_handler(E_ERR,routine,message)
end subroutine

subroutine load_or_build_search(path,uuid,lon,lat,s)
character(len=*), intent(in) :: path,uuid
real(r8), intent(in) :: lon(:),lat(:)
type(horizontal_search_type), intent(out) :: s
logical :: exists
integer :: status
call require(len_trim(path)>0,'Empty horizontal_search_cache_file')
inquire(file=trim(path),exist=exists)
if (exists) then
   call read_search_cache(path,uuid,lon,lat,s)
   call error_handler(E_MSG,routine,'Loaded horizontal search: '//trim(path))
else
   call error_handler(E_MSG,routine,'Constructing horizontal triangulation and octree')
   call build_horizontal_search(lon,lat,s)
   call write_search_cache(path,uuid,lon,lat,s,status)
   if (status==1) then
      ! Another invocation won publication. Never overwrite its product.
      call read_search_cache(path,uuid,lon,lat,s)
      call error_handler(E_MSG,routine,'Loaded concurrently published horizontal search: '//trim(path))
   else if (status/=0) then
      call error_handler(E_WARN,routine,'Could not publish cache; using in-memory search: '//trim(path))
   else
      call error_handler(E_MSG,routine,'Saved horizontal search: '//trim(path))
   endif
endif
end subroutine

subroutine read_search_cache(path,uuid,lon,lat,s)
character(len=*), intent(in) :: path,uuid
real(r8), intent(in) :: lon(:),lat(:)
type(horizontal_search_type), intent(out) :: s
integer :: nc,dc,dt,db,d3,d2,n,t,b,value,id,i
real(r8), allocatable :: saved_lon(:),saved_lat(:)
logical, allocatable :: used(:)
character(len=128) :: saved_uuid
call nc_check(nf90_open(trim(path),nf90_nowrite,nc),routine,'open '//trim(path))
call nc_check(nf90_get_att(nc,nf90_global,'schema_version',value),routine)
call require(value==schema_version,'Unsupported search cache schema; remove cache or select another path')
call nc_check(nf90_get_att(nc,nf90_global,'algorithm_version',value),routine)
call require(value==algorithm_version,'Unsupported search algorithm version')
call nc_check(nf90_get_att(nc,nf90_global,'octree_depth',value),routine)
call require(value==OCTREE_DEPTH,'Incompatible octree depth')
call nc_check(nf90_get_att(nc,nf90_global,'uuidOfHGrid',saved_uuid),routine)
call require(trim(saved_uuid)==trim(uuid),'Search cache horizontal UUID mismatch')
call dimension(nc,'cell',dc,n)
call require(n==size(lon).and.n==size(lat),'Search cache cell count mismatch')
call dimension(nc,'triangle',dt,t)
call require(t>=1.and.t<=2*n-4,'Invalid cached triangle count')
call dimension(nc,'box_entry',db,b)
call require(b>=2.and.b<=t+1,'Invalid cached box count')
call dimension(nc,'xyz_axis',d3,value)
call require(value==3,'Invalid xyz_axis dimension')
call dimension(nc,'pair',d2,value)
call require(value==2,'Invalid pair dimension')
allocate(saved_lon(n),saved_lat(n),s%triangles%point_cell(n),s%triangles%vertex(3,t),s%triangles%xyz(3,n))
allocate(s%tree%object_idx(t),s%tree%object_data(t),s%tree%obj_min(t,3),s%tree%obj_max(t,3),s%tree%box(b,2))
call variable(nc,'cell_lon',nf90_double,[dc],id)
call nc_check(nf90_get_var(nc,id,saved_lon),routine,'cell_lon')
call variable(nc,'cell_lat',nf90_double,[dc],id)
call nc_check(nf90_get_var(nc,id,saved_lat),routine,'cell_lat')
call variable(nc,'point_cell',nf90_int,[dc],id)
call nc_check(nf90_get_var(nc,id,s%triangles%point_cell),routine,'point_cell')
call variable(nc,'vertex',nf90_int,[d3,dt],id)
call nc_check(nf90_get_var(nc,id,s%triangles%vertex),routine,'vertex')
call variable(nc,'xyz',nf90_double,[d3,dc],id)
call nc_check(nf90_get_var(nc,id,s%triangles%xyz),routine,'xyz')
call variable(nc,'object_idx',nf90_int,[dt],id)
call nc_check(nf90_get_var(nc,id,s%tree%object_idx),routine,'object_idx')
call variable(nc,'object_data',nf90_int,[dt],id)
call nc_check(nf90_get_var(nc,id,s%tree%object_data),routine,'object_data')
call variable(nc,'obj_min',nf90_double,[dt,d3],id)
call nc_check(nf90_get_var(nc,id,s%tree%obj_min),routine,'obj_min')
call variable(nc,'obj_max',nf90_double,[dt,d3],id)
call nc_check(nf90_get_var(nc,id,s%tree%obj_max),routine,'obj_max')
call variable(nc,'box',nf90_int,[db,d2],id)
call nc_check(nf90_get_var(nc,id,s%tree%box),routine,'box')
call variable(nc,'root_range',nf90_double,[d2,d3],id)
call nc_check(nf90_get_var(nc,id,s%tree%brange0),routine,'root_range')

call nc_check(nf90_get_att(nc,nf90_global,'delta',s%tree%delta),routine)
call nc_check(nf90_close(nc),routine)
call require(all(saved_lon==lon).and.all(saved_lat==lat),'Search cache ordered cell coordinates mismatch')
call require(all(ieee_is_finite(s%triangles%xyz)),'Nonfinite cached xyz')
call require(all(s%triangles%point_cell>=1.and.s%triangles%point_cell<=n),'Invalid cached cell mapping')
allocate(used(n))
used=.false.
used(s%triangles%point_cell)=.true.
call require(all(used),'Cached cell mapping is not a permutation')
call require(all(s%triangles%vertex>=1.and.s%triangles%vertex<=n),'Invalid cached triangle index')
do i=1,t
   call require(s%triangles%vertex(1,i)/=s%triangles%vertex(2,i).and. &
                s%triangles%vertex(1,i)/=s%triangles%vertex(3,i).and. &
                s%triangles%vertex(2,i)/=s%triangles%vertex(3,i),'Repeated cached triangle vertex')
enddo
call require(all(s%tree%object_idx>=1.and.s%tree%object_idx<=t),'Invalid cached object index')
deallocate(used)
allocate(used(t))
used=.false.
used(s%tree%object_idx)=.true.
call require(all(used),'Cached object indices are not a permutation')
call require(all(s%tree%object_data==[(i,i=1,t)]),'Invalid cached triangle mapping')
call require(all(ieee_is_finite(s%tree%obj_min)).and.all(ieee_is_finite(s%tree%obj_max)), &
             'Nonfinite cached object bounds')
call require(all(s%tree%obj_min<=s%tree%obj_max),'Reversed cached object bounds')
call require(all(s%tree%brange0(1,:)==-1._r8).and.all(s%tree%brange0(2,:)==1._r8),'Invalid cached root range')
call require(s%tree%delta==2._r8/(2._r8**(OCTREE_DEPTH-1)),'Invalid cached octree resolution')
call require(all(s%tree%box(1:b-1,1)>=1).and.all(s%tree%box(1:b-1,1)<=(8**OCTREE_DEPTH-1)/7), &
             'Invalid cached box ID')
call require(s%tree%box(b,1)==(8**OCTREE_DEPTH-1)/7+1,'Invalid cached box sentinel')
call require(all(s%tree%box(2:b,1)>s%tree%box(1:b-1,1)),'Unsorted cached box IDs')
call require(s%tree%box(1,2)==1.and.s%tree%box(b,2)==t+1,'Invalid cached offset endpoints')
call require(all(s%tree%box(2:b,2)>s%tree%box(1:b-1,2)),'Invalid cached box offsets')
end subroutine

subroutine dimension(nc,name,id,n)
integer,intent(in) :: nc
character(len=*),intent(in) :: name
integer,intent(out) :: id,n
call nc_check(nf90_inq_dimid(nc,name,id),routine,name)
call nc_check(nf90_inquire_dimension(nc,id,len=n),routine,name)
call require(n>0,'Empty cache dimension '//name)
end subroutine

subroutine variable(nc,name,xtype,dims,id)
integer,intent(in) :: nc,xtype,dims(:)
character(len=*),intent(in) :: name
integer,intent(out) :: id
integer :: actual_type,rank,actual_dims(nf90_max_var_dims)
call nc_check(nf90_inq_varid(nc,name,id),routine,name)
call nc_check(nf90_inquire_variable(nc,id,xtype=actual_type,ndims=rank,dimids=actual_dims),routine,name)
call require(actual_type==xtype.and.rank==size(dims),'Wrong cache variable type or rank: '//name)
call require(all(actual_dims(:rank)==dims),'Wrong cache variable dimensions: '//name)
end subroutine

subroutine write_search_cache(path,uuid,lon,lat,s,status)
character(len=*), intent(in) :: path,uuid
real(r8), intent(in) :: lon(:),lat(:)
type(horizontal_search_type), intent(in) :: s
! 0 published, 1 another publisher won, 2 I/O failure (caller may keep memory).
integer, intent(out) :: status
integer :: nc,rc,dc,dt,db,d3,d2,ids(11)
integer(c_int) :: fd,ignored
character(kind=c_char,len=:),allocatable :: temporary
logical :: exists
status=2
temporary=trim(path)//'.tmp.XXXXXX'//c_null_char
fd=c_mkstemp(temporary)
if (fd<0) return
ignored=c_close(fd)
rc=nf90_create(temporary(:len(temporary)-1),ior(nf90_clobber,nf90_64bit_offset),nc)
if (rc==nf90_noerr) then
   writing: block
      rc=nf90_put_att(nc,nf90_global,'schema_version',schema_version)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_att(nc,nf90_global,'algorithm_version',algorithm_version)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_att(nc,nf90_global,'octree_depth',OCTREE_DEPTH)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_att(nc,nf90_global,'uuidOfHGrid',trim(uuid))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_att(nc,nf90_global,'delta',s%tree%delta)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_att(nc,nf90_global,'ICON_revision','8597da45ef4b86323f3fb844caedc4ae5e1ffc01')
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_att(nc,nf90_global,'fortran_support_revision','34246610b17db29f214fb2f95ca1c9087f09b89c')
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_att(nc,nf90_global,'iconmath_revision','622c01ba78e4fbb19fdf5f3c9f18a90b1e551389')
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_dim(nc,'cell',size(lon),dc)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_dim(nc,'triangle',size(s%triangles%vertex,2),dt)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_dim(nc,'box_entry',size(s%tree%box,1),db)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_dim(nc,'xyz_axis',3,d3)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_dim(nc,'pair',2,d2)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'cell_lon',nf90_double,[dc],ids(1))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'cell_lat',nf90_double,[dc],ids(2))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'point_cell',nf90_int,[dc],ids(3))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'vertex',nf90_int,[d3,dt],ids(4))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'xyz',nf90_double,[d3,dc],ids(5))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'object_idx',nf90_int,[dt],ids(6))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'object_data',nf90_int,[dt],ids(7))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'obj_min',nf90_double,[dt,d3],ids(8))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'obj_max',nf90_double,[dt,d3],ids(9))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'box',nf90_int,[db,d2],ids(10))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_def_var(nc,'root_range',nf90_double,[d2,d3],ids(11))
      if (rc/=nf90_noerr) exit writing
      rc=nf90_enddef(nc)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(1),lon)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(2),lat)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(3),s%triangles%point_cell)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(4),s%triangles%vertex)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(5),s%triangles%xyz)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(6),s%tree%object_idx)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(7),s%tree%object_data)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(8),s%tree%obj_min)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(9),s%tree%obj_max)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(10),s%tree%box)
      if (rc/=nf90_noerr) exit writing
      rc=nf90_put_var(nc,ids(11),s%tree%brange0)
      if (rc/=nf90_noerr) exit writing
   end block writing
   if (rc==nf90_noerr) then
      rc=nf90_close(nc)
   else
      ignored=nf90_close(nc)
   endif
endif
if (rc==nf90_noerr) then
   ignored=c_link(temporary,trim(path)//c_null_char)
   if (ignored==0) then
      status=0
   else
      inquire(file=trim(path),exist=exists)
      if (exists) status=1
   endif
endif
ignored=c_unlink(temporary)
end subroutine
end module
