! DART software - Copyright UCAR. See the DART license for terms of use.
program triangulation_check
use types_mod, only : r8
use utilities_mod, only : initialize_utilities, finalize_utilities, error_handler, E_ERR
use icon_triangulation_mod, only : auxiliary_triangulation_type, build_auxiliary_triangulation
use icon_geometry_mod, only : init_icon_geometry, end_icon_geometry, geometry
implicit none
type(auxiliary_triangulation_type) :: tri, repeated
real(r8), allocatable :: lon(:),lat(:)
real(r8) :: begin,finish
integer :: n,u,i
character(len=512) :: input,output
call initialize_utilities('triangulation_check')
call get_command_argument(1,input)
call get_command_argument(2,output)
if (len_trim(output)==0) call error_handler(E_ERR,'triangulation_check','Usage: points.txt|--geometry output.txt')
if (trim(input)=='--geometry') then
   call init_icon_geometry()
   lon=geometry%clon
   lat=geometry%clat
   n=size(lon)
else
   open(newunit=u,file=trim(input),status='old')
   read(u,*) n
   allocate(lon(n),lat(n))
   do i=1,n
      read(u,*) lon(i),lat(i)
   enddo
   close(u)
endif
call cpu_time(begin)
call build_auxiliary_triangulation(lon,lat,tri)
call cpu_time(finish)
write(*,'(A,I0,A,I0,A,F12.4)') 'TRIANGULATION cells=',n,', triangles=',size(tri%vertex,2), &
                            ', construction CPU seconds=',finish-begin
! Separate reconstruction verifies repeatability; no static/module cache can mask it.
call build_auxiliary_triangulation(lon,lat,repeated)
if (size(tri%vertex,2)/=size(repeated%vertex,2)) &
   call error_handler(E_ERR,'triangulation_check','Repeat construction changed triangle count')
if (any(tri%vertex/=repeated%vertex) .or. any(tri%point_cell/=repeated%point_cell) .or. &
    any(tri%xyz/=repeated%xyz)) call error_handler(E_ERR,'triangulation_check','Repeat construction changed product')
call repeated%clear()
open(newunit=u,file=trim(output),status='new')
write(u,*) n,size(tri%vertex,2)
do i=1,n
   write(u,'(I0,3(1X,ES26.17E3))') tri%point_cell(i),tri%xyz(:,i)
enddo
do i=1,size(tri%vertex,2)
   write(u,*) tri%vertex(:,i)
enddo
close(u)
call tri%clear()
call tri%clear()
if (trim(input)=='--geometry') call end_icon_geometry()
write(*,'(A)') 'triangulation_check PASS'
call finalize_utilities()
end program
