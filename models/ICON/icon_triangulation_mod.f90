! DART software - Copyright UCAR. See the DART license for terms of use.
module icon_triangulation_mod

! Serial adaptation of ICON's repartition setup, global spherical-cap branch.
! No MPI ownership, octree, observation weights or persistent cache in this step.
use types_mod, only : r8, PI
use utilities_mod, only : error_handler, E_ERR
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use icon_delaunay_types_mod, only : t_point_list, t_triangulation, point, spherical_cap
use icon_delaunay_mod, only : triangulate_serial
use icon_spherical_coordinates_mod, only : t_cartesian_coordinates, t_geographical_coordinates, &
                                          gc2cc_sphere, cc2gc_sphere
implicit none
private
public :: auxiliary_triangulation_type, build_auxiliary_triangulation

type auxiliary_triangulation_type
   ! One-based indices into sorted points; triangle row/corner order is native.
   integer, allocatable :: vertex(:,:) ! (3,ntriangles)
   integer, allocatable :: point_cell(:) ! sorted point -> original cell ID
   real(r8), allocatable :: xyz(:,:) ! original, UNPERTURBED centers in sorted order
contains
   procedure :: clear => clear_triangulation
end type

contains

subroutine build_auxiliary_triangulation(lon, lat, result)
real(r8), intent(in) :: lon(:), lat(:) ! radians, original native cell order
! intent(out) releases a previously allocated product automatically.
type(auxiliary_triangulation_type), intent(out) :: result
type(t_point_list) :: points
type(t_triangulation) :: triangles
type(t_geographical_coordinates) :: gc
type(t_cartesian_coordinates) :: cc
real(r8), allocatable :: original(:,:)
integer :: i, j, n
logical, allocatable :: used(:)
character(len=*), parameter :: routine='build_auxiliary_triangulation'

n=size(lon)
if (size(lat)/=n .or. n<4) call error_handler(E_ERR,routine,'Require at least four cell centers with matching lon/lat sizes')
if (.not. all(ieee_is_finite(lon)) .or. .not. all(ieee_is_finite(lat))) &
   call error_handler(E_ERR,routine,'Nonfinite cell-center coordinates')
if (any(abs(lat)>PI/2._r8) .or. any(abs(lon)>2._r8*PI)) &
   call error_handler(E_ERR,routine,'Cell-center coordinates must be geographic radians')
allocate(original(3,n),result%point_cell(n),result%xyz(3,n))
call points%initialize()
call points%reserve(n)
do i=1,n
   gc=t_geographical_coordinates(lon(i),lat(i))
   cc=gc2cc_sphere(gc)
   original(:,i)=cc%x
   ! Reproduce ICON's Cartesian/geographic round trip, not just a direct
   ! addition to the supplied longitude. It matters near symmetry and poles.
   gc=cc2gc_sphere(cc)
   gc%lon=gc%lon+cos(gc%lon)*sin(gc%lat)*1.e-10_r8
   cc=gc2cc_sphere(gc)
   call points%push_back(point(cc%x(1),cc%x(2),cc%x(3),iindex=i-1))
enddo
! Preserve the native equal-key quicksort, including its permutation of ties.
call points%quicksort()
do i=0,n-1
   result%point_cell(i+1)=points%a(i)%gindex+1
   points%a(i)%gindex=i
enddo
call triangles%initialize()
call triangulate_serial(points,triangles,spherical_cap(point(0._r8,0._r8,1._r8),-1._r8), &
                        ignore_completeness=.false.)
! A triangulation on the sphere cannot exceed 2*N-4 triangles. These are
! structural checks, not new geometric tolerances or Delaunay validation.
if (triangles%nentries<1 .or. triangles%nentries>2*n-4) &
   call error_handler(E_ERR,routine,'ICON construction returned an invalid triangle count')
allocate(result%vertex(3,triangles%nentries))
do i=0,triangles%nentries-1
   ! Kernel normalizes corner rotation in sorted-point space, not native-cell space.
   result%vertex(:,i+1)=triangles%a(i)%p+1
enddo
if (any(result%vertex<1) .or. any(result%vertex>n)) &
   call error_handler(E_ERR,routine,'Invalid point index returned by ICON construction')
allocate(used(n))
used=.false.
do j=1,size(result%vertex,2)
   used(result%vertex(:,j))=.true.
   if (result%vertex(1,j)==result%vertex(2,j) .or. result%vertex(1,j)==result%vertex(3,j) .or. &
       result%vertex(2,j)==result%vertex(3,j)) call error_handler(E_ERR,routine,'Repeated vertex in triangle')
enddo
if (.not. all(used)) call error_handler(E_ERR,routine,'ICON construction omitted input cell centers')
! Restore the original physical positions. Perturbation chooses connectivity only.
result%xyz=original(:,result%point_cell)
call triangles%destructor()
call points%destructor()
end subroutine

subroutine clear_triangulation(this)
class(auxiliary_triangulation_type), intent(inout) :: this
if (allocated(this%vertex)) deallocate(this%vertex)
if (allocated(this%point_cell)) deallocate(this%point_cell)
if (allocated(this%xyz)) deallocate(this%xyz)
end subroutine
end module icon_triangulation_mod
