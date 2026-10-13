! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See licenses/ICON-AUTHORS.TXT for a list of authors
! See licenses/ICON-BSD-3-Clause.txt for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Numerical routines from ICON mo_intp_lonlat_baryctr, release 2026.04.
module icon_search_math_mod
use types_mod, only : wp => r8
use icon_spherical_coordinates_mod, only : t_geographical_coordinates, t_cartesian_coordinates, &
   gc2cc => gc2cc_sphere, project_point_to_plane, cc_dot_product
use icon_delaunay_types_mod, only : t_point, ccw_spherical, operator(+), operator(*), operator(/)
implicit none
private
public :: compute_barycentric_coords, inside_triangle, inside_scaled_triangle
contains
    SUBROUTINE compute_barycentric_coords(pt, v, u)
      TYPE(t_geographical_coordinates), INTENT(IN)  :: pt                  !< query point (longitude/latitude)
      !> triangle vertices, first index is x,y,z,
      !! second index is 1st,2nd or 3rd vertex
      REAL(wp),                         INTENT(IN)  :: v(3,3)
      REAL(wp),                         INTENT(OUT) :: u(3)                !< barycentric coordinates (dim: 3)
      ! local variables
      TYPE(t_cartesian_coordinates) :: a,b,c,p
      REAL(wp) :: B1(2), B2(2), r(2), det
      INTEGER :: i

      ! Cartesian coordinates of triangle vertices
      a%x(:) = v(:, 1)
      b%x(:) = v(:, 2)
      c%x(:) = v(:, 3)
      ! query point and project in onto the plane (p1,p2,p3):
      p = project_point_to_plane (gc2cc(pt), a,b,c)
      ! solve linear system for the barycentric coordinates
      !
      ! define matrix and right hand side
      DO i = 1, 3
        a%x(i) = a%x(i) - c%x(i)
        b%x(i) = b%x(i) - c%x(i)
        p%x(i) = p%x(i) - c%x(i)
      END DO
      B1(1:2) = (/ cc_dot_product(a,a) , cc_dot_product(a,b) /)   ! (a-c)^2     ,  (a-c)*(b-c)
      B2(1:2) = (/ B1(2)               , cc_dot_product(b,b) /)   ! (a-c)*(b-c) ,  (b-c)^2
      r(1:2)  = (/ cc_dot_product(p,a) , cc_dot_product(p,b) /)   ! [ (p-c)*(a-c), (p-c)*(b-c) ]
      ! solve using Cramer's rule:
      det     = B2(2)*B1(1) - B1(2)*B2(1)
      u(1)    = ( r(1)*B2(2) - r(2)*B1(2) )/det
      u(2)    = ( r(2)*B1(1) - r(1)*B2(1) )/det
      u(3)    = 1._wp - (u(1) + u(2))
    END SUBROUTINE compute_barycentric_coords
    FUNCTION inside_triangle(v, v1,v2,v3)
      INTEGER :: inside_triangle
      REAL(wp),       INTENT(IN)     :: v(3)          !< query point
      TYPE(t_point),  INTENT(IN)     :: v1,v2,v3      !< vertex longitudes/latitudes
      ! local variables
      INTEGER       :: c1,c2,c3
      TYPE(t_point) :: p

      p%x = v(1)
      p%y = v(2)
      p%z = v(3)

      c1  = ccw_spherical(v1,v2,p)
      c2  = ccw_spherical(v2,v3,p)
      c3  = ccw_spherical(v3,v1,p)
      IF ((c1 + c2 + c3) == 3 .OR. (c1 + c2 + c3) == -3) THEN
        inside_triangle = 1 ! numerically clearly inside
      ELSEIF (c1 == 0 .OR. c2 == 0 .OR. c3 == 0) THEN
        inside_triangle = 0 ! numerically indefinite
      ELSE
        inside_triangle = -1 ! numerically clearly outside
      ENDIF
    END FUNCTION inside_triangle
    FUNCTION inside_scaled_triangle(v, v1, v2, v3, rscale)
      ! Test if point v is inside a triangle around (v1, v2, v3), that is scaled in size by given rscale.
      ! The area of (v1, v2, v3) is approximately increased to (1+rscale)**2
      INTEGER :: inside_scaled_triangle
      REAL(wp),       INTENT(IN)     :: v(3)          !< query point
      REAL(wp),       INTENT(IN)     :: rscale        !< factor by which to increase the triangle size (>0 for growth)
      TYPE(t_point),  INTENT(IN)     :: v1,v2,v3      !< vertex longitudes/latitudes
      ! local variables
      TYPE(t_point) :: center, g1, g2, g3
      REAL(wp) :: rscale1

      center = (v1 + v2 + v3)/3._wp

      ! grow the v-triangle from the center by rscale
      ! \vec{g_i} = \vec{v_i} + scale * (\vec{v_i} - \vec{center}) = (scale+1) * \vec{v_i} - scale * \vec{center}
      rscale1 = 1._wp + rscale
      g1 = rscale1 * v1 + (-rscale) * center
      g2 = rscale1 * v2 + (-rscale) * center
      g3 = rscale1 * v3 + (-rscale) * center

      inside_scaled_triangle = inside_triangle(v, g1, g2, g3)
    END FUNCTION inside_scaled_triangle
end module
