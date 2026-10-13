! DART software - Copyright UCAR. See the DART license for terms of use.
module icon_search_init_mod
! Collective model lifecycle; numerical modules have no MPI dependencies.
use types_mod, only : r8
use mpi_utilities_mod, only : my_task_id, task_count, array_broadcast
use icon_geometry_mod, only : geometry
use icon_horizontal_search_mod, only : horizontal_search_type
use icon_search_cache_mod, only : load_or_build_search
implicit none
private
public :: init_icon_search, end_icon_search, horizontal_search
type(horizontal_search_type), protected :: horizontal_search
logical :: initialized=.false.
interface distribute
   module procedure distribute_integer_1, distribute_integer_2, distribute_real_2
end interface
contains
subroutine init_icon_search(path)
character(len=*), intent(in) :: path
integer :: n,t,b
real(r8) :: lengths(3),resolution(1)
if (initialized) return
if (my_task_id()==0) then
   call load_or_build_search(path,geometry%horizontal_uuid,geometry%clon,geometry%clat,horizontal_search)
endif
if (task_count()>1) then
   if (my_task_id()==0) lengths=real([size(horizontal_search%triangles%point_cell), &
      size(horizontal_search%triangles%vertex,2),size(horizontal_search%tree%box,1)],r8)
   call array_broadcast(lengths,0)
   n=nint(lengths(1)); t=nint(lengths(2)); b=nint(lengths(3))
   if (my_task_id()/=0) then
      allocate(horizontal_search%triangles%point_cell(n),horizontal_search%triangles%vertex(3,t), &
               horizontal_search%triangles%xyz(3,n))
      allocate(horizontal_search%tree%object_idx(t),horizontal_search%tree%object_data(t), &
               horizontal_search%tree%obj_min(t,3),horizontal_search%tree%obj_max(t,3), &
               horizontal_search%tree%box(b,2))
   endif
   call distribute(horizontal_search%triangles%point_cell)
   call distribute(horizontal_search%triangles%vertex)
   call distribute(horizontal_search%triangles%xyz)
   call distribute(horizontal_search%tree%object_idx)
   call distribute(horizontal_search%tree%object_data)
   call distribute(horizontal_search%tree%obj_min)
   call distribute(horizontal_search%tree%obj_max)
   call distribute(horizontal_search%tree%box)
   call distribute(horizontal_search%tree%brange0)
   if (my_task_id()==0) resolution(1)=horizontal_search%tree%delta
   call array_broadcast(resolution,0)
   horizontal_search%tree%delta=resolution(1)
endif
initialized=.true.
end subroutine

! DART broadcasts r8 arrays. Every 32-bit integer is represented exactly in
! r8; convert in bounded chunks rather than duplicating whole large arrays.
subroutine distribute_integer_1(a)
integer, intent(inout) :: a(:)
real(r8) :: buffer(65536)
integer :: first,last,n
 do first=1,size(a),size(buffer)
   last=min(size(a),first+size(buffer)-1)
   n=last-first+1
   if (my_task_id()==0) buffer(:n)=real(a(first:last),r8)
   call array_broadcast(buffer(:n),0)
   if (my_task_id()/=0) a(first:last)=nint(buffer(:n))
 enddo
end subroutine
subroutine distribute_integer_2(a)
integer,intent(inout),target,contiguous :: a(:,:)
integer,pointer :: flat(:)
flat(1:size(a))=>a
call distribute_integer_1(flat)
end subroutine
subroutine distribute_real_2(a)
real(r8),intent(inout),target,contiguous :: a(:,:)
real(r8),pointer :: flat(:)
flat(1:size(a))=>a
call array_broadcast(flat,0)
end subroutine
subroutine end_icon_search()
call horizontal_search%clear()
initialized=.false.
end subroutine
end module
