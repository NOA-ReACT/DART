! DART software - Copyright UCAR. See the DART license for terms of use.
program geometry_check

! Exercise geometry independently of state registration and observation operators.
use types_mod, only : r8
use utilities_mod, only : initialize_utilities, finalize_utilities, error_handler, E_ERR
use icon_geometry_mod, only : init_icon_geometry, end_icon_geometry, geometry
implicit none
integer :: cells, excluded

call initialize_utilities('geometry_check')
call init_icon_geometry()
cells = geometry%ncells
excluded = count(geometry%excluded)
! Repeated initialization must be harmless (no reallocation or file reread).
call init_icon_geometry()
if (any(abs(sum(geometry%xyz**2,dim=1)-1.0_r8)>1.0e-12_r8)) &
   call error_handler(E_ERR,'geometry_check','Cartesian centers are not unit vectors')
write(*,'(A,2(1X,ES24.16))') 'z_mc range:', minval(geometry%z_mc), maxval(geometry%z_mc)
write(*,'(A,2(1X,ES24.16))') 'ground range:', minval(geometry%z_ifc(:,geometry%nlevels+1)), &
                                           maxval(geometry%z_ifc(:,geometry%nlevels+1))
call end_icon_geometry()
if (allocated(geometry%z_mc) .or. geometry%ncells/=0) &
   call error_handler(E_ERR,'geometry_check','Geometry cleanup failed')
call init_icon_geometry()
if (geometry%ncells/=cells .or. count(geometry%excluded)/=excluded) &
   call error_handler(E_ERR,'geometry_check','Geometry reinitialization changed dimensions or mask')
call end_icon_geometry()
write(*,'(A)') 'geometry_check PASS'
call finalize_utilities()
end program
