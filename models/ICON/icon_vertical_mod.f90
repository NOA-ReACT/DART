! DART software - Copyright UCAR. See the DART license for terms of use.
module icon_vertical_mod
! MPAS find_height_bounds default policy, with ICON top-down storage.
use types_mod, only : r8
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
implicit none
private
public :: height_bracket
contains
pure subroutine height_bracket(z, height, lower, upper, fraction, valid)
real(r8), intent(in) :: z(:), height
integer, intent(out) :: lower, upper
real(r8), intent(out) :: fraction
logical, intent(out) :: valid
integer :: k
lower=0
upper=0
fraction=0._r8
valid=.false.
if (size(z)<2) return
if (.not. ieee_is_finite(height)) return
if (height<z(size(z)) .or. height>z(1)) return
! Scan upwards physically: an exact interior level belongs to the interval
! below it, just as in MPAS. Geometry initialization validates monotonicity.
do k=size(z),2,-1
   if (height<=z(k-1)) then
      lower=k
      upper=k-1
      fraction=(height-z(k))/(z(k-1)-z(k))
      valid=.true.
      return
   endif
enddo
end subroutine
end module
