! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See licenses/ICON-AUTHORS.TXT for a list of authors
! See licenses/ICON-BSD-3-Clause.txt for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Extracted from iconmath 622c01ba78e4fbb19fdf5f3c9f18a90b1e551389.
! Keep the Cartesian/geographic round trip used in ICON point perturbation.
module icon_spherical_coordinates_mod
use types_mod, only : wp => r8
implicit none
private
public :: cc2gc_sphere, gc2cc_sphere, t_cartesian_coordinates, t_geographical_coordinates
public :: project_point_to_plane, cc_dot_product
real(wp), parameter :: dbl_eps=epsilon(1._wp)
real(wp), parameter :: pi_2=1.57079632679489661923132169163975144_wp
type t_cartesian_coordinates
real(wp) :: x(3)
end type
type t_geographical_coordinates
real(wp) :: lon,lat
end type
contains
  ELEMENTAL FUNCTION cc2gc_sphere(p_x) RESULT(p_pos)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x ! Cart. coordinates
    TYPE(t_geographical_coordinates)          :: p_pos ! geo. coordinates

    REAL(wp)                                :: z_x, z_y, z_z, z_r
    !-----------------------------------------------------------------------

    z_x = p_x%x(1)
    z_y = p_x%x(2)
    z_z = p_x%x(3)

    z_r = z_x*z_x + z_y*z_y
    z_r = SQRT(z_r)

    IF (ABS(z_r) < dbl_eps) THEN ! one of the poles

      IF (z_z > 0.0_wp) THEN
        p_pos%lat = pi_2
      ELSE
        p_pos%lat = -1._wp*pi_2
      END IF
      p_pos%lon = 0._wp

    ELSE

      p_pos%lat = ATAN2(z_z, z_r)

      IF (ABS(z_x) < dbl_eps) THEN ! z_x == 0 ?

        IF (z_y >= 0.0_wp) THEN
          p_pos%lon = pi_2
        ELSE
          p_pos%lon = -1._wp*pi_2
        END IF

      ELSE
        p_pos%lon = ATAN2(z_y, z_x)
      END IF

    END IF

  END FUNCTION cc2gc_sphere
  ELEMENTAL FUNCTION gc2cc_sphere(p_pos) RESULT(p_x)
    TYPE(t_geographical_coordinates), INTENT(in) :: p_pos ! geo. coordinates
    TYPE(t_cartesian_coordinates)                :: p_x ! Cart. coordinates
    REAL(wp) :: z_cln, z_sln, z_clt, z_slt

    z_sln = SIN(p_pos%lon)
    z_cln = COS(p_pos%lon)
    z_slt = SIN(p_pos%lat)
    z_clt = COS(p_pos%lat)

    p_x%x(1) = z_cln*z_clt
    p_x%x(2) = z_sln*z_clt
    p_x%x(3) = z_slt

  END FUNCTION gc2cc_sphere
  ELEMENTAL FUNCTION project_point_to_plane(in_point, to_v0, to_v1, to_v2) RESULT(projected_point)
    TYPE(t_cartesian_coordinates), INTENT(in) :: in_point, to_v0, to_v1, to_v2
    TYPE(t_cartesian_coordinates) :: projected_point

    TYPE(t_cartesian_coordinates) :: in_p, v1, v2, norm_to_plane

    ! translate everybody to the center
    v1%x = to_v1%x - to_v0%x
    v2%x = to_v2%x - to_v0%x
    in_p%x = in_point%x - to_v0%x

    norm_to_plane = vector_product(v1, v2)
    norm_to_plane%x=norm_to_plane%x/SQRT(DOT_PRODUCT(norm_to_plane%x,norm_to_plane%x))
    norm_to_plane%x = norm_to_plane%x* &
      & DOT_PRODUCT(norm_to_plane%x, (v1%x - in_p%x))

    projected_point%x = in_p%x + norm_to_plane%x + to_v0%x

  END FUNCTION project_point_to_plane
  ELEMENTAL FUNCTION cc_dot_product(cc_x1, cc_x2) RESULT(p_prod)
    !$ACC ROUTINE SEQ
    TYPE(t_cartesian_coordinates), INTENT(in):: cc_x1, cc_x2 ! cart. coordinates

    REAL(wp) :: p_prod ! scalar product of cart.coordinates

    p_prod = DOT_PRODUCT(cc_x1%x, cc_x2%x)

  END FUNCTION cc_dot_product
  ELEMENTAL FUNCTION vector_product(x0, x1) RESULT(x2)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1

    TYPE(t_cartesian_coordinates) :: x2

    !-----------------------------------------------------------------------
    x2%x(1) = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
    x2%x(2) = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
    x2%x(3) = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)

  END FUNCTION vector_product
end module icon_spherical_coordinates_mod
