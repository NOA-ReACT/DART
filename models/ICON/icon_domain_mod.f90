! DART software - Copyright UCAR. See the DART license for terms of use.
module icon_domain_mod
! Native-domain membership and nearest-cell service. This index is ephemeral;
! the persistent auxiliary interpolation cache is deliberately independent.
use types_mod, only : r8
use icon_geometry_mod, only : geometry
use xyz_location_mod, only : xyz_location_type, xyz_get_close_type, xyz_set_location, &
   xyz_get_close_init, xyz_get_close_destroy, xyz_find_nearest
use icon_octree_mod, only : t_range_octree, octree_finalize, octree_count_point, octree_query_point
use icon_horizontal_search_mod, only : build_triangle_index
use icon_delaunay_types_mod, only : t_point, point, ccw_spherical
use icon_spherical_coordinates_mod, only : t_cartesian_coordinates, t_geographical_coordinates, gc2cc_sphere
implicit none
private
public :: init_icon_domain, end_icon_domain, native_contains, nearest_icon_cell, native_triangle_contains
type(t_range_octree) :: native_tree
type(xyz_get_close_type) :: nearest_index
type(xyz_location_type), allocatable :: cell_locations(:)
real(r8), allocatable :: vertex_xyz(:,:)
logical :: initialized=.false.
contains
subroutine init_icon_domain()
integer :: i,j,neighbor
real(r8) :: padding,xyz(3)
type(t_cartesian_coordinates) :: cc
if (initialized) return
allocate(cell_locations(geometry%ncells))
padding=0._r8
do i=1,geometry%ncells
   xyz=geometry%sphere_radius*geometry%xyz(:,i)
   cell_locations(i)=xyz_set_location(xyz(1),xyz(2),xyz(3))
   do j=1,3
      neighbor=geometry%neighbors(i,j)
      if (neighbor<1) cycle
      padding=max(padding,norm2(geometry%xyz(:,i)-geometry%xyz(:,neighbor))*geometry%sphere_radius)
   enddo
enddo
! Padding only sets the spatial bins; it is not an observation cutoff.
call xyz_get_close_init(nearest_index,padding,geometry%ncells,cell_locations)
if (.not. geometry%global_grid) then
   allocate(vertex_xyz(3,geometry%nvertices))
   do i=1,geometry%nvertices
      cc=gc2cc_sphere(t_geographical_coordinates(geometry%vlon(i),geometry%vlat(i)))
      vertex_xyz(:,i)=cc%x
   enddo
   call build_triangle_index(vertex_xyz,transpose(geometry%vertices),native_tree)
endif
initialized=.true.
end subroutine

subroutine nearest_icon_cell(xyz,cell,status)
real(r8), intent(in) :: xyz(3) ! unit sphere
integer, intent(out) :: cell,status
real(r8) :: p(3)
p=geometry%sphere_radius*xyz
call xyz_find_nearest(nearest_index,xyz_set_location(p(1),p(2),p(3)),cell_locations,cell,status)
end subroutine

logical function native_contains(xyz) result(inside)
real(r8), intent(in) :: xyz(3)
integer, allocatable :: candidates(:)
integer :: i,n
inside=.true.
if (geometry%global_grid) return
inside=.false.
n=octree_count_point(native_tree,xyz)
allocate(candidates(n))
call octree_query_point(native_tree,xyz,candidates)
do i=1,n
   if (native_triangle_contains(vertex_xyz(:,geometry%vertices(candidates(i),:)),xyz)) then
      inside=.true.
      return
   endif
enddo
end function

pure logical function native_triangle_contains(vertices,xyz) result(inside)
real(r8), intent(in) :: vertices(3,3),xyz(3)
type(t_point) :: v(3),p
integer :: i,orientation
p=point(xyz(1),xyz(2),xyz(3))
do i=1,3
   v(i)=point(vertices(1,i),vertices(2,i),vertices(3,i))
enddo
orientation=ccw_spherical(v(1),v(2),v(3))
inside=.false.
! An unresolved orientation cannot establish an inward hemisphere.
if (orientation==0) return
do i=1,3
   if (orientation*ccw_spherical(v(i),v(mod(i,3)+1),p)<0) return
enddo
inside=.true.
end function

subroutine end_icon_domain()
if (.not. initialized) return
call xyz_get_close_destroy(nearest_index)
deallocate(cell_locations)
if (allocated(vertex_xyz)) then
   deallocate(vertex_xyz)
   call octree_finalize(native_tree)
endif
initialized=.false.
end subroutine
end module
