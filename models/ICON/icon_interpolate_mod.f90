! DART software - Copyright UCAR. See the DART license for terms of use.
module icon_interpolate_mod
use types_mod, only : r8, MISSING_R8, deg2rad
use location_mod, only : location_type, get_location, query_location, VERTISHEIGHT
use obs_kind_mod, only : QTY_VERTICAL_VELOCITY
use ensemble_manager_mod, only : ensemble_type
use distributed_state_mod, only : get_state
use state_structure_mod, only : get_varid_from_kind, get_dart_vector_index
use icon_geometry_mod, only : geometry
use icon_search_init_mod, only : horizontal_search
use icon_horizontal_search_mod, only : query_horizontal_search, SEARCH_OK
use icon_spherical_coordinates_mod, only : t_cartesian_coordinates, t_geographical_coordinates, gc2cc_sphere
use icon_domain_mod, only : native_contains, nearest_icon_cell
use icon_vertical_mod, only : height_bracket
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
implicit none
private
public :: interpolate_icon, observation_stencil_type, build_observation_stencil, evaluate_icon_stencil
public :: INTERP_OK, INTERP_UNSUPPORTED, INTERP_UNAVAILABLE, INTERP_DOMAIN, INTERP_BOUNDARY
public :: INTERP_SEARCH, INTERP_VERTICAL, INTERP_INVALID_DATA, INTERP_LOCATION
integer, parameter :: INTERP_OK=0, INTERP_UNSUPPORTED=1, INTERP_UNAVAILABLE=2, INTERP_DOMAIN=3
integer, parameter :: INTERP_BOUNDARY=4, INTERP_SEARCH=5, INTERP_VERTICAL=6, INTERP_INVALID_DATA=7
integer, parameter :: INTERP_LOCATION=8
type observation_stencil_type
   integer :: cells(3)=0, lower(3)=0, upper(3)=0
   real(r8) :: weights(3)=0._r8, fraction(3)=0._r8
end type
contains
subroutine interpolate_icon(state_handle,ens_size,location,qty,domain_id,expected_obs,istatus)
type(ensemble_type), intent(in) :: state_handle
integer, intent(in) :: ens_size,qty,domain_id
type(location_type), intent(in) :: location
real(r8), intent(out) :: expected_obs(ens_size)
integer, intent(out) :: istatus(ens_size)
type(observation_stencil_type) :: stencil
integer :: varid,status
expected_obs=MISSING_R8
istatus=INTERP_UNSUPPORTED
if (nint(query_location(location))/=VERTISHEIGHT .or. qty==QTY_VERTICAL_VELOCITY .or. qty<0) return
varid=get_varid_from_kind(domain_id,qty)
istatus=INTERP_UNAVAILABLE
if (varid<1) return
call build_observation_stencil(location,stencil,status)
istatus=status
if (status/=INTERP_OK) return
call evaluate_icon_stencil(state_handle,ens_size,domain_id,varid,stencil,expected_obs,istatus)
end subroutine

subroutine build_observation_stencil(location,stencil,status)
type(location_type), intent(in) :: location
type(observation_stencil_type), intent(out) :: stencil
integer, intent(out) :: status
real(r8) :: loc(3)
type(t_cartesian_coordinates) :: cc
integer :: nearest,rc,triangle,j
logical :: valid
status=INTERP_UNSUPPORTED
if (nint(query_location(location))/=VERTISHEIGHT) return
loc=get_location(location)
status=INTERP_LOCATION
if (.not. all(ieee_is_finite(loc))) return
loc(1:2)=loc(1:2)*deg2rad
cc=gc2cc_sphere(t_geographical_coordinates(loc(1),loc(2)))
status=INTERP_DOMAIN
if (.not. native_contains(cc%x)) return
call nearest_icon_cell(cc%x,nearest,rc)
status=INTERP_SEARCH
if (rc/=0 .or. nearest<1 .or. nearest>geometry%ncells) return
status=INTERP_BOUNDARY
if (geometry%excluded(nearest)) return
call query_horizontal_search(horizontal_search,loc(1),loc(2),triangle,stencil%cells,stencil%weights,rc)
status=INTERP_SEARCH
if (rc/=SEARCH_OK) return
status=INTERP_BOUNDARY
if (any(geometry%excluded(stencil%cells))) return
! All three columns must be admissible, including corners with zero weight.
status=INTERP_VERTICAL
do j=1,3
   call height_bracket(geometry%z_mc(stencil%cells(j),:),loc(3), &
      stencil%lower(j),stencil%upper(j),stencil%fraction(j),valid)
   if (.not. valid) return
enddo
status=INTERP_OK
end subroutine

subroutine evaluate_icon_stencil(state_handle,ens_size,domain_id,varid,stencil,expected_obs,istatus)
! Caller supplies a valid stencil. No collectives or direct ensemble storage:
! get_state uses the state window established by DART's caller.
type(ensemble_type), intent(in) :: state_handle
integer, intent(in) :: ens_size,domain_id,varid
type(observation_stencil_type), intent(in) :: stencil
real(r8), intent(out) :: expected_obs(ens_size)
integer, intent(out) :: istatus(ens_size)
real(r8) :: values(ens_size),column(ens_size),vertical_weight(2),coefficient
integer :: j,k,m,levels(2)
expected_obs=0._r8
istatus=INTERP_OK
do j=1,3
   column=0._r8
   levels=[stencil%lower(j),stencil%upper(j)]
   vertical_weight=[1._r8-stencil%fraction(j),stencil%fraction(j)]
   do k=1,2
      coefficient=stencil%weights(j)*vertical_weight(k)
      if (coefficient==0._r8) cycle
      values=get_state(get_dart_vector_index(stencil%cells(j),levels(k),1,domain_id,varid),state_handle)
      do m=1,ens_size
         if (istatus(m)/=INTERP_OK) cycle
         if (.not. ieee_is_finite(values(m)) .or. values(m)==MISSING_R8) then
            istatus(m)=INTERP_INVALID_DATA
         else
            column(m)=column(m)+vertical_weight(k)*values(m)
         endif
      enddo
   enddo
   do m=1,ens_size
      if (istatus(m)==INTERP_OK) expected_obs(m)=expected_obs(m)+stencil%weights(j)*column(m)
   enddo
enddo
do m=1,ens_size
   if (.not. ieee_is_finite(expected_obs(m))) istatus(m)=INTERP_INVALID_DATA
   if (istatus(m)/=INTERP_OK) expected_obs(m)=MISSING_R8
enddo
end subroutine
end module
