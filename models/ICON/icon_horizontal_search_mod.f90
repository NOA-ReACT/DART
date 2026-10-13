! DART software - Copyright UCAR. See the DART license for terms of use.
module icon_horizontal_search_mod
! Serial search over all auxiliary triangles. A successful query is NOT a
! native-domain or boundary-mask acceptance decision.
use types_mod, only : r8, PI
use icon_triangulation_mod, only : auxiliary_triangulation_type, build_auxiliary_triangulation
use icon_octree_mod, only : t_range_octree, octree_init, octree_finalize, octree_count_point, octree_query_point
use icon_delaunay_types_mod, only : t_point, point, sagitta_on_unit_sphere
use icon_spherical_coordinates_mod, only : t_geographical_coordinates, t_cartesian_coordinates, gc2cc_sphere
use icon_search_math_mod, only : compute_barycentric_coords, inside_triangle, inside_scaled_triangle
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
implicit none
private
public :: horizontal_search_type, build_horizontal_search, query_horizontal_search, search_candidates
public :: SEARCH_OK, SEARCH_NOT_FOUND, SEARCH_INVALID_LOCATION, SEARCH_TOO_MANY_HITS, SEARCH_BAD_WEIGHTS
integer, parameter :: SEARCH_OK=0, SEARCH_NOT_FOUND=1, SEARCH_INVALID_LOCATION=2
integer, parameter :: SEARCH_TOO_MANY_HITS=3, SEARCH_BAD_WEIGHTS=4
integer, parameter :: max_hits=99, max_scale=13
real(r8), parameter :: inside_tol=1.e-4_r8, scale_factor=1.e6_r8*epsilon(1._r8)
type horizontal_search_type
   type(auxiliary_triangulation_type) :: triangles
   type(t_range_octree) :: tree
contains
   procedure :: clear => clear_search
end type
contains
subroutine build_horizontal_search(lon,lat,search)
real(r8), intent(in) :: lon(:),lat(:)
type(horizontal_search_type), intent(out) :: search
real(r8), allocatable :: pmin(:,:),pmax(:,:)
real(r8) :: brange(2,3),sagitta,xyz(3)
type(t_point) :: previous,current
integer :: i,j,k,n
call build_auxiliary_triangulation(lon,lat,search%triangles)
n=size(search%triangles%vertex,2)
allocate(pmin(n,3),pmax(n,3))
do i=1,n
   pmin(i,:)=99._r8
   pmax(i,:)=-99._r8
   sagitta=0._r8
   xyz=search%triangles%xyz(:,search%triangles%vertex(3,i))
   current=point(xyz(1),xyz(2),xyz(3))
   do j=1,3
      previous=current
      xyz=search%triangles%xyz(:,search%triangles%vertex(j,i))
      current=point(xyz(1),xyz(2),xyz(3))
      do k=1,3
         pmin(i,k)=min(pmin(i,k),xyz(k))
         pmax(i,k)=max(pmax(i,k),xyz(k))
      enddo
      sagitta=max(sagitta,sagitta_on_unit_sphere(current,previous))
   enddo
   pmin(i,:)=pmin(i,:)-2._r8*sagitta
   pmax(i,:)=pmax(i,:)+2._r8*sagitta
enddo
brange(1,:)=-1._r8
brange(2,:)=1._r8
! Explicit one-based triangle mapping, in native insertion order.
call octree_init(search%tree,brange,pmin,pmax,[(i,i=1,n)])
end subroutine

subroutine search_candidates(search,lon,lat,candidates,status)
type(horizontal_search_type), intent(in) :: search
real(r8), intent(in) :: lon,lat
integer, allocatable, intent(out) :: candidates(:)
integer, intent(out) :: status
type(t_cartesian_coordinates) :: cc
integer :: n
status=SEARCH_INVALID_LOCATION
allocate(candidates(0))
if (.not. ieee_is_finite(lon) .or. .not. ieee_is_finite(lat)) return
if (abs(lat)>PI/2._r8 .or. abs(lon)>2._r8*PI) return
cc=gc2cc_sphere(t_geographical_coordinates(lon,lat))
n=octree_count_point(search%tree,cc%x)
deallocate(candidates)
allocate(candidates(n))
call octree_query_point(search%tree,cc%x,candidates)
status=SEARCH_OK
if (n>max_hits) status=SEARCH_TOO_MANY_HITS
end subroutine

subroutine query_horizontal_search(search,lon,lat,triangle,cells,weights,status,retry)
type(horizontal_search_type), intent(in) :: search
real(r8), intent(in) :: lon,lat
integer, intent(out) :: triangle,cells(3),status
real(r8), intent(out) :: weights(3)
integer, intent(out), optional :: retry
integer, allocatable :: candidates(:)
integer :: i,j,k,scale,classification
real(r8) :: v(3,3),raw(3),checksum
type(t_point) :: p(3)
type(t_geographical_coordinates) :: gc
type(t_cartesian_coordinates) :: cc
triangle=0
cells=0
weights=0._r8
if (present(retry)) retry=-1
call search_candidates(search,lon,lat,candidates,status)
if (status/=SEARCH_OK) return
gc=t_geographical_coordinates(lon,lat)
cc=gc2cc_sphere(gc)
status=SEARCH_NOT_FOUND
do scale=0,max_scale
   do i=1,size(candidates)
      j=candidates(i)
      if (j<1) cycle
      v=search%triangles%xyz(:,search%triangles%vertex(:,j))
      do k=1,3
         p(k)=point(v(1,k),v(2,k),v(3,k))
      enddo
      if (scale==0) then
         call compute_barycentric_coords(gc,v,raw)
         if (.not. all(ieee_is_finite(raw))) then
            candidates(i)=-1
            cycle
         endif
         if (.not.(all(raw>=-inside_tol).and.all(raw<=1._r8+inside_tol))) then
            candidates(i)=-1
            cycle
         endif
         classification=inside_triangle(cc%x,p(1),p(2),p(3))
      else
         classification=inside_scaled_triangle(cc%x,p(1),p(2),p(3),2**scale*scale_factor)
      endif
      if (classification==-1) then
         candidates(i)=-1
         cycle
      endif
      if (classification==0 .and. scale/=max_scale) cycle
      ! ICON exits on the first accepted candidate, including the final
      ! uncertain retry. Final weights always use original physical vertices.
      call compute_barycentric_coords(gc,v,raw)
      if (.not. all(ieee_is_finite(raw))) then
         status=SEARCH_BAD_WEIGHTS
         return
      endif
      do k=1,3
         raw(k)=max(min(1._r8,raw(k)),0._r8)
      enddo
      checksum=0._r8
      do k=1,3
         checksum=checksum+raw(k)
      enddo
      if (checksum<=0._r8) then
         status=SEARCH_BAD_WEIGHTS
         return
      endif
      weights=raw/checksum
      triangle=j
      cells=search%triangles%point_cell(search%triangles%vertex(:,j))
      status=SEARCH_OK
      if (present(retry)) retry=scale
      return
   enddo
enddo
end subroutine

subroutine clear_search(this)
class(horizontal_search_type), intent(inout) :: this
call this%triangles%clear()
if (allocated(this%tree%object_idx)) call octree_finalize(this%tree)
end subroutine
end module
