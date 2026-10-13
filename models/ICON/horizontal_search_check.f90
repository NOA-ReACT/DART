! DART software - Copyright UCAR. See the DART license for terms of use.
program horizontal_search_check
! points.txt cache.nc queries.txt results.txt [roundtrip.nc]
! or --geometry cache.nc queries.txt results.txt (input.nml geometry required).
! Empty query file (count=0) is the standalone preparation mode.
use types_mod, only : r8
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities, my_task_id
use utilities_mod, only : error_handler, E_ERR
use icon_horizontal_search_mod, only : horizontal_search_type, search_candidates, query_horizontal_search
use icon_search_cache_mod, only : load_or_build_search, write_search_cache
use icon_search_init_mod, only : init_icon_search, end_icon_search, horizontal_search
use icon_geometry_mod, only : init_icon_geometry, end_icon_geometry
implicit none
type(horizontal_search_type) :: local
real(r8),allocatable :: lon(:),lat(:)
real(r8) :: qlon,qlat,weights(3),start,finish
integer,allocatable :: candidates(:)
integer :: u,v,i,n,status,t,cells(3),retry,ios
character(len=1024) :: points_path,cache_path,queries_path,output_path,roundtrip
character(len=32) :: rank_string
logical :: from_geometry
call initialize_mpi_utilities('horizontal_search_check')
call get_command_argument(1,points_path)
call get_command_argument(2,cache_path)
call get_command_argument(3,queries_path)
call get_command_argument(4,output_path)
call get_command_argument(5,roundtrip)
if (len_trim(output_path)==0) call error_handler(E_ERR,'horizontal_search_check', &
   'Usage: points.txt|--geometry cache.nc queries.txt results.txt [roundtrip.nc]')
from_geometry=trim(points_path)=='--geometry'
call cpu_time(start)
if (from_geometry) then
   call init_icon_geometry()
   call init_icon_search(trim(cache_path))
   call init_icon_search(trim(cache_path))
else
   open(newunit=u,file=trim(points_path),status='old')
   read(u,*) n
   allocate(lon(n),lat(n))
   do i=1,n
      read(u,*) lon(i),lat(i)
   enddo
   close(u)
   call load_or_build_search(trim(cache_path),'synthetic-search-fixture',lon,lat,local)
   if (len_trim(roundtrip)>0) then
      call write_search_cache(trim(roundtrip),'synthetic-search-fixture',lon,lat,local,status)
      if (status/=0) call error_handler(E_ERR,'horizontal_search_check','Roundtrip publication failed')
   endif
endif
call cpu_time(finish)
write(*,'(A,F12.6)') 'SEARCH INITIALIZATION CPU SECONDS ',finish-start
write(rank_string,'(I0)') my_task_id()
open(newunit=v,file=trim(output_path)//'.'//trim(rank_string),status='replace')
open(newunit=u,file=trim(queries_path),status='old')
read(u,*) n
do i=1,n
   read(u,*,iostat=ios) qlon,qlat
   if (ios/=0) call error_handler(E_ERR,'horizontal_search_check','Invalid query record')
   if (from_geometry) then
      call search_candidates(horizontal_search,qlon,qlat,candidates,status)
      call query_horizontal_search(horizontal_search,qlon,qlat,t,cells,weights,status,retry)
   else
      call search_candidates(local,qlon,qlat,candidates,status)
      call query_horizontal_search(local,qlon,qlat,t,cells,weights,status,retry)
   endif
   write(v,'(*(I0,1X))') i,status,t,cells,retry,size(candidates),candidates
   write(v,'(3(ES26.17E3,1X))') weights
enddo
close(u)
close(v)
call local%clear()
call local%clear()
if (from_geometry) then
   call end_icon_search()
   call end_icon_search()
   call end_icon_geometry()
endif
write(*,'(A)') 'horizontal_search_check PASS'
call finalize_mpi_utilities()
end program
