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

! Adapted from ICON 8597da45ef4b86323f3fb844caedc4ae5e1ffc01.
! Serial numerical subset; MPI, I/O and unused list types removed.
! Type and object declarations, operators and utility routines for
! the construction of a Delaunay
! triangulation on the sphere.

MODULE icon_delaunay_types_mod

  USE types_mod, ONLY: wp => r8
  USE utilities_mod, ONLY: error_handler, E_ERR



  IMPLICIT NONE
  INTEGER, PARAMETER :: qp=selected_real_kind(33,4931), SUCCESS=0

  PRIVATE
  PUBLIC :: t_edge, t_point, t_triangle
  PUBLIC :: t_spherical_cap
  PUBLIC :: t_triangulation
  PUBLIC :: t_point_list
  PUBLIC :: point, triangle, point_list, spherical_cap
  PUBLIC :: OPERATOR(<), OPERATOR(/), OPERATOR(*), OPERATOR(==), OPERATOR(+)
  PUBLIC :: circum_circle_spherical, ccw_spherical, sagitta_on_unit_sphere

  CHARACTER(LEN=*), PARAMETER :: modname = 'icon_delaunay_types_mod'


  ! --------------------------------------------------------------------
  ! OBJECT DECLARATIONS
  ! --------------------------------------------------------------------

  !> type declaration: single triangle edge
  TYPE t_edge
    INTEGER        :: p1,p2                       !< edge vertex indices
  END TYPE t_edge

  !> type declaration: single point in R^3
  TYPE t_point
    REAL(wp)       :: x,y,z, ps
    INTEGER        :: gindex                      !< global index
  CONTAINS
    PROCEDURE :: norm2          => point_norm2
    PROCEDURE :: spherical_dist => point_spherical_dist
  END TYPE t_point

  !> type declaration: single triangle
  TYPE t_triangle
    INTEGER        :: p(0:2)                      !< vertex indices
    INTEGER        :: oedge                       !< local index of edge opposite to ghost point (if present)
    INTEGER        :: complete                    !< 0: flag not set; 1: triangle must not be visited again; 2: triangle invalid
    TYPE(t_point)  :: cc                          !< circumcircle center
    REAL(wp)       :: r                           !< circumcircle radius
    REAL(wp)       :: cap_distance                !< min distance from subset center to circumcircle
    REAL(wp)       :: rdiscard                    !< max distance from subset center, projected onto sort direction
  CONTAINS
    PROCEDURE, PUBLIC :: edge                 => triangle_edge
    PROCEDURE, PUBLIC :: compute_circumcenter => triangle_compute_circumcenter
    PROCEDURE, PUBLIC :: normalize_indices    => triangle_normalize_indices
  END TYPE t_triangle


  !> type declaration: spherical cap
  TYPE t_spherical_cap
    TYPE (t_point) :: point, sorting_direction
    REAL(wp)       :: radius                      !< cap radius (radius<0 means "global")
  END TYPE t_spherical_cap

  !> type declaration: list of abstract objects
  TYPE, ABSTRACT :: t_object_list
    INTEGER :: nentries
  CONTAINS
    PRIVATE
    PROCEDURE, PUBLIC :: initialize => initialize_object_list
    PROCEDURE(routine1),  DEFERRED, PASS :: destructor !< destructor... not all compilers support "FINAL"
  END TYPE t_object_list

  !> type declaration: list of abstract objects
  TYPE, ABSTRACT, EXTENDS(t_object_list) :: t_sortable_list
  CONTAINS
    PRIVATE
    PROCEDURE(arrayfct),    DEFERRED, PASS :: swap       !< swap operator for generic sort algorithm
    PROCEDURE(binaryop),    DEFERRED, PASS :: cmp        !< "<" operator for generic sort algorithm
    PROCEDURE(routine3),    DEFERRED, PASS :: reserve    !< re-allocate data structure
    PROCEDURE, PUBLIC :: quicksort  => quicksort_sortable_list
  END TYPE t_sortable_list

  !> type declaration: list of points
  TYPE, EXTENDS(t_sortable_list) :: t_point_list
    TYPE(t_point), ALLOCATABLE :: a(:)
  CONTAINS
    PROCEDURE :: swap       => swap_point_list
    PROCEDURE :: cmp        => cmp_point_list
    PROCEDURE :: destructor => destructor_point_list
    PROCEDURE :: reserve    => reserve_point_list
    PROCEDURE :: resize     => resize_point_list
    PROCEDURE :: push_back  => push_back_point_list
    PROCEDURE, PUBLIC :: quicksort  => quicksort_point_list
  END TYPE t_point_list

  !> type declaration: list of triangles
  TYPE, EXTENDS(t_sortable_list) :: t_triangulation
    TYPE(t_triangle), ALLOCATABLE :: a(:)
  CONTAINS
    PROCEDURE :: swap         => swap_triangulation
    PROCEDURE :: cmp          => cmp_triangulation
    PROCEDURE :: destructor   => destructor_triangulation
    PROCEDURE :: reserve      => reserve_triangulation
    PROCEDURE :: resize       => resize_triangulation
    PROCEDURE :: push_back    => push_back_triangulation

  END TYPE t_triangulation

  !> type declaration: list of triangles

  ! --------------------------------------------------------------------
  ! DECLARATION OF OPERATORS AND TYPE-BOUND PROCEDURES
  ! --------------------------------------------------------------------

  ABSTRACT INTERFACE
    SUBROUTINE routine1(this)
      IMPORT :: t_object_list
      CLASS(t_object_list), INTENT(INOUT) :: this
    END SUBROUTINE routine1
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE routine3(this, new_size)
      IMPORT :: t_sortable_list
      CLASS(t_sortable_list), INTENT(INOUT) :: this
      INTEGER,                INTENT(IN)    :: new_size
    END SUBROUTINE routine3
  END INTERFACE

  ABSTRACT INTERFACE
    LOGICAL PURE FUNCTION binaryop(this, i,j)
      IMPORT :: t_sortable_list
      CLASS(t_sortable_list), INTENT(IN) :: this
      INTEGER, INTENT(IN) :: i,j
    END FUNCTION binaryop
  END INTERFACE

  ABSTRACT INTERFACE
    PURE SUBROUTINE arrayfct(this, i,j)
      IMPORT :: t_sortable_list
      CLASS(t_sortable_list), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: i,j
    END SUBROUTINE arrayfct
  END INTERFACE

  !> operator: comparison of two edges
  INTERFACE OPERATOR ( == )
    MODULE PROCEDURE t_edge_equal
    MODULE PROCEDURE t_point_equal
    MODULE PROCEDURE t_triangle_equal
  END INTERFACE

  !> operator: subtraction of two points
  INTERFACE OPERATOR ( - )
    MODULE PROCEDURE t_point_minus
  END INTERFACE

  !> operator: subtraction of two points
  INTERFACE OPERATOR ( + )
    MODULE PROCEDURE t_point_plus
  END INTERFACE

  !> operator: dot product and product with scalar
  INTERFACE OPERATOR ( * )
    MODULE PROCEDURE t_point_mul
    MODULE PROCEDURE t_point_mul_scalar
  END INTERFACE

  !> operator: division of point by scalar
  INTERFACE OPERATOR ( / )
    MODULE PROCEDURE t_point_div_scalar
  END INTERFACE

  !> operator: comparison of points and triangles (wrt. sorting direction)
  INTERFACE OPERATOR ( < )
    MODULE PROCEDURE t_point_cmp
    MODULE PROCEDURE t_triangle_cmp
  END INTERFACE

CONTAINS

  ! --------------------------------------------------------------------
  ! GENERAL GEOMETRIC SUBROUTINES
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> Return TRUE if a point (xp,yp) is inside the circumcircle made up
  !  of the points ip[1,2,3].
  !
  !  See, e.g.,
  !  Renka, R. J. Interpolation of Data on the Surface of a Sphere
  !               ACM Trans. Math. Softw., ACM, 1984, 10, 417-436
  !  Renka's STRIPACK algorithm (http://www.netlib.org/toms/772)
  PURE &
  FUNCTION circum_circle_spherical(p, pxyz, ip)
    LOGICAL :: circum_circle_spherical
    TYPE (t_point),      INTENT(IN)   :: p
    TYPE (t_point_list), INTENT(IN)   :: pxyz
    INTEGER,             INTENT(IN)   :: ip(0:2)
    ! local variables
    REAL(wp) :: d1_x, d1_y, d1_z,d2_x, d2_y, d2_z,d3_x, d3_y, d3_z, d1, d2

    ! p lies above the plane of (p1,p3,p2) iff p2 lies above the plane
    ! of (p3,p1,p) iff Det(p2-p,p3-p,p1-p) = (p2-p,p3-p X p1-p) > 0.
    !
    ! det =      (pxyz%a(ip(1))%x - p%x)*(d3%y*d1%z - d1%y*d3%z) &
    !   &     -  (pxyz%a(ip(1))%y - p%y)*(d3%x*d1%z - d1%x*d3%z) &
    !   &     +  (pxyz%a(ip(1))%z - p%z)*(d3%x*d1%y - d1%x*d3%y)
    d1_x = pxyz%a(ip(0))%x - p%x
    d1_y = pxyz%a(ip(0))%y - p%y
    d1_z = pxyz%a(ip(0))%z - p%z
    d2_x = pxyz%a(ip(1))%x - p%x
    d2_y = pxyz%a(ip(1))%y - p%y
    d2_z = pxyz%a(ip(1))%z - p%z
    d3_x = pxyz%a(ip(2))%x - p%x
    d3_y = pxyz%a(ip(2))%y - p%y
    d3_z = pxyz%a(ip(2))%z - p%z

    d1 = d2_x*d3_y*d1_z  +  d2_y*d1_x*d3_z + d2_z*d3_x*d1_y
    d2 = d2_y*d3_x*d1_z + d2_x*d1_y*d3_z + d2_z*d1_x*d3_y
    IF (ABS(d1-d2) < 1.e-12_wp) THEN
      circum_circle_spherical = circum_circle_spherical_q128(p, pxyz, ip)
    ELSE
      circum_circle_spherical = (d1 > d2)
    END IF
  END FUNCTION circum_circle_spherical


  ! --------------------------------------------------------------------
  !> Return TRUE if a point (xp,yp) is inside the circumcircle made up
  !  of the points ip[1,2,3].
  !
  !  See, e.g.,
  !  Renka, R. J. Interpolation of Data on the Surface of a Sphere
  !               ACM Trans. Math. Softw., ACM, 1984, 10, 417-436
  !  Renka's STRIPACK algorithm (http://www.netlib.org/toms/772)
  PURE FUNCTION circum_circle_spherical_q128(p, pxyz, ip)
    LOGICAL :: circum_circle_spherical_q128
    TYPE (t_point),      INTENT(IN)   :: p
    TYPE (t_point_list), INTENT(IN)   :: pxyz
    INTEGER,             INTENT(IN)   :: ip(0:2)
    ! local variables
    REAL(qp) :: d1_x, d1_y, d1_z,d2_x, d2_y, d2_z,d3_x, d3_y, d3_z

    ! p lies above the plane of (p1,p3,p2) iff p2 lies above the plane
    ! of (p3,p1,p) iff Det(p2-p,p3-p,p1-p) = (p2-p,p3-p X p1-p) > 0.
    !
    ! det =      (pxyz%a(ip(1))%x - p%x)*(d3%y*d1%z - d1%y*d3%z) &
    !   &     -  (pxyz%a(ip(1))%y - p%y)*(d3%x*d1%z - d1%x*d3%z) &
    !   &     +  (pxyz%a(ip(1))%z - p%z)*(d3%x*d1%y - d1%x*d3%y)
    d1_x = pxyz%a(ip(0))%x - p%x
    d1_y = pxyz%a(ip(0))%y - p%y
    d1_z = pxyz%a(ip(0))%z - p%z
    d2_x = pxyz%a(ip(1))%x - p%x
    d2_y = pxyz%a(ip(1))%y - p%y
    d2_z = pxyz%a(ip(1))%z - p%z
    d3_x = pxyz%a(ip(2))%x - p%x
    d3_y = pxyz%a(ip(2))%y - p%y
    d3_z = pxyz%a(ip(2))%z - p%z

    IF ( d2_x*d3_y*d1_z  +  d2_y*d1_x*d3_z + d2_z*d3_x*d1_y &
      > d2_y*d3_x*d1_z + d2_x*d1_y*d3_z + d2_z*d1_x*d3_y ) THEN
      circum_circle_spherical_q128 = .TRUE.
    ELSE
      circum_circle_spherical_q128 = .FALSE.
    END IF
  END FUNCTION circum_circle_spherical_q128


  ! --------------------------------------------------------------------
  !> Locates a point relative to a directed arc. Let v1, v2, and v3 be
  !  distinct points on the sphere, and denote by v1 -> v2 the
  !  geodesic connecting v1 and v2 and directed toward v2. This test
  !  determines which of the two hemispheres defined by v1 -> v2
  !  contains v3, or, in other words, if we "turn left" when going
  !  from v1 to v2 to v3.
  PURE FUNCTION ccw_spherical(v1,v2,v3)
    INTEGER :: ccw_spherical ! -1: clearly left, 0: numerically uncertain, 1: clearly right
    TYPE (t_point), INTENT(IN)  :: v1,v2,v3
    REAL(wp), PARAMETER :: FLOAT_ERROR_MARGIN = 1.1e-14_wp
    REAL(wp) :: ccw

    ! det(v1,v2,v3) = <v1 x v2, v3> = | v1 x v2 | cos(a)
    !
    ! where a is the angle between v3 and the normal to the plane
    ! defined by v1 and v2.

    ccw =           v3%x*(v1%y*v2%z - v2%y*v1%z) &
      &        -    v3%y*(v1%x*v2%z - v2%x*v1%z) &
      &        +    v3%z*(v1%x*v2%y - v2%x*v1%y)


    ! we apply a static error of
    !   | e - e'| <= 3*2^-48
    ! to decide if a floating-point evaluation e' of an expression e
    ! has the correct sign, see Section 2.2 of
    !
    ! Burnikel, C.; Funke, S. & Seel, M.
    ! "Exact geometric computation using Cascading"
    ! International Journal of Computational Geometry & Applications,
    ! World Scientific, 2001, 11, 245-266
    IF (ccw < -FLOAT_ERROR_MARGIN) THEN
      ccw_spherical = -1
    ELSE IF (ccw > FLOAT_ERROR_MARGIN) THEN
      ccw_spherical = 1
    ELSE
      ccw_spherical = 0
    END IF
  END FUNCTION ccw_spherical



  ! --------------------------------------------------------------------
  !> Computes sagitta of great circle arc, which connects points a1 and a2 on a unit sphere
  REAL(wp) PURE FUNCTION sagitta_on_unit_sphere(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1, a2
    ! The sagitta of an arc is distance from the center of the arc to the center of its base.
    ! The base is the chord that connects a1 and a2.
    !
    ! This function assumes 3D points on a unit sphere.
    !
    !       a1          Let 2*q be the chord length between a1 and a2.
    !     , '\          Let s be sagitta that we are looking for and h the height of the
    !   ,   | \         triangle (center, a1, a2).
    !  ,   q|  \ 1      Then s = 1 - h
    ! ,     |   \         len(a1) = len(a2) = 1
    ! ,_____|___ \cen     h**2 + q**2 = 1**2 ! Pythagoras
    ! , s   | h  /ter     (2*q)**2 = len(a1)**2 + len(a2)**2 - 2 * len(a1) * len(a1) &
    !  ,   q|   /                  & * cos(angle(a1,a2)) ! Law of cosines
    !   ,   |  /1         cos(angle(a1,a2)) = abs(dot(a1, a2)) / (len(a1) * len(a2)) &
    !     , | /                    & ! definition of scalar product
    !       '/          Combine and simplify: q**2 = 0.5 - 0.5 * abs(dot(a1, a2))
    !       a2          Thus sagitta s = 1 - SQRT(0.5 * (1 + ABS(dot(a1, a2))))

    sagitta_on_unit_sphere = 1._wp - SQRT( 0.5_wp * (1._wp + ABS(a1 * a2)) )
  END FUNCTION sagitta_on_unit_sphere


  ! --------------------------------------------------------------------
  ! DEFINITION OF OVERLOADED OPERATORS
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> operator: comparison of two edges
  LOGICAL PURE FUNCTION t_edge_equal(e1,e2)
    TYPE(t_edge), INTENT(IN) :: e1,e2
    ! note: we test opposite directions first, since this is the
    ! situation for neighboring triangles
    t_edge_equal =                                            &
      &    ( (e1%p1 == e2%p2) .AND. (e1%p2 == e2%p1) )  .OR.  &
      &    ( (e1%p1 == e2%p1) .AND. (e1%p2 == e2%p2) )
  END FUNCTION t_edge_equal


  ! --------------------------------------------------------------------
  !> operator: comparison of two triangles
  LOGICAL PURE FUNCTION t_triangle_equal(t1,t2)
    TYPE(t_triangle), INTENT(IN) :: t1,t2
    t_triangle_equal =                                        &
      &    ( (t1%p(0) == t2%p(0)) .AND. (t1%p(1) == t2%p(1)) .AND. (t1%p(2) == t2%p(2)) )
  END FUNCTION t_triangle_equal


  ! --------------------------------------------------------------------
  !> operator: subtraction of two points
  TYPE(t_point) PURE FUNCTION t_point_minus(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1,a2
    t_point_minus = t_point(x=a1%x-a2%x, y=a1%y-a2%y, z=a1%z-a2%z, &
      &                     ps=-2._wp, gindex=0)
  END FUNCTION t_point_minus


  ! --------------------------------------------------------------------
  !> operator: addition of two points
  TYPE(t_point) PURE FUNCTION t_point_plus(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1,a2
    t_point_plus = t_point(x=a1%x+a2%x, y=a1%y+a2%y, z=a1%z+a2%z, &
      &                    ps=-2._wp, gindex=0)
  END FUNCTION t_point_plus


  ! --------------------------------------------------------------------
  !> operator: dot product
  REAL(wp) PURE FUNCTION t_point_mul(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1,a2
    t_point_mul = a1%x*a2%x + a1%y*a2%y + a1%z*a2%z
  END FUNCTION t_point_mul


  ! --------------------------------------------------------------------
  !> operator: product of scalar and point
  TYPE(t_point) PURE FUNCTION t_point_mul_scalar(v,t)
    TYPE(t_point), INTENT(IN) :: t
    REAL(wp),      INTENT(IN) :: v
    t_point_mul_scalar = t_point(x=v*t%x, y=v*t%y, z=v*t%z, &
      &                   ps=t%ps, gindex=t%gindex)
  END FUNCTION t_point_mul_scalar


  ! --------------------------------------------------------------------
  !> operator: division of point by scalar
  TYPE(t_point) PURE FUNCTION t_point_div_scalar(t,v)
    TYPE(t_point), INTENT(IN) :: t
    REAL(wp),      INTENT(IN) :: v
    t_point_div_scalar = t_point(x=t%x/v, y=t%y/v, z=t%z/v, &
      &                   ps=t%ps, gindex=t%gindex)
  END FUNCTION t_point_div_scalar


  ! --------------------------------------------------------------------
  ! DEFINITION OF "t_spherical_cap" SUBROUTINES
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> constructor for t_spherical_cap type
  PURE FUNCTION spherical_cap(pp, rr)
    TYPE(t_spherical_cap) :: spherical_cap
    TYPE(t_point), INTENT(IN)  :: pp
    REAL(wp),      INTENT(IN)  :: rr
    TYPE(t_point) :: point
    point = pp/pp%norm2()
    spherical_cap = t_spherical_cap(point=point, sorting_direction=(point/(-1._wp)), radius=rr)
  END FUNCTION spherical_cap


  ! --------------------------------------------------------------------
  ! DEFINITION OF "t_point" SUBROUTINES
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> constructor for t_point type
  PURE FUNCTION point(ix,iy,iz,iindex,ips)
    TYPE(t_point) :: point
    REAL(wp), INTENT(IN) :: ix, iy, iz
    INTEGER,  INTENT(IN), OPTIONAL :: iindex
    REAL(wp), INTENT(IN), OPTIONAL :: ips
    INTEGER  :: index
    REAL(wp) :: ps
    index =  0
    ps    = -2._wp
    IF (PRESENT(iindex))  index=iindex
    IF (PRESENT(ips))     ps=ips
    point = t_point(x=ix,y=iy,z=iz, ps=ps,gindex=index)
  END FUNCTION point

  ! --------------------------------------------------------------------
  !> operator: comparison of points (wrt. sorting direction)
  LOGICAL PURE FUNCTION t_point_cmp(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1,a2
    t_point_cmp = (a2%ps > a1%ps)
  END FUNCTION t_point_cmp


  ! --------------------------------------------------------------------
  !> operator: comparison of points (wrt. sorting direction)
  LOGICAL PURE FUNCTION t_point_equal(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1,a2
    t_point_equal = (a2%ps == a1%ps)
  END FUNCTION t_point_equal


  ! --------------------------------------------------------------------
  !> compute Euclidean norm in R^3
  REAL(wp) PURE FUNCTION point_norm2(this)
    CLASS(t_point), INTENT(IN) :: this
    point_norm2 = SQRT(this%x*this%x + this%y*this%y + this%z*this%z)
  END FUNCTION point_norm2


  ! --------------------------------------------------------------------
  !> compute spherical distance
  REAL(wp) PURE FUNCTION point_spherical_dist(a1,a2)
    CLASS(t_point), INTENT(IN) :: a1
    TYPE(t_point),  INTENT(IN) :: a2
    point_spherical_dist = ACOS(a1%x*a2%x + a1%y*a2%y + a1%z*a2%z)
  END FUNCTION point_spherical_dist


  ! --------------------------------------------------------------------
  ! DEFINITION OF "t_triangle" SUBROUTINES
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> constructor for t_triangle type
  PURE FUNCTION triangle(ip1, ip2, ip3, ioedge, icomplete)
    TYPE(t_triangle) :: triangle
    INTEGER, INTENT(IN)           :: ip1, ip2, ip3
    INTEGER, INTENT(IN), OPTIONAL :: ioedge
    INTEGER, INTENT(IN), OPTIONAL :: icomplete
    INTEGER :: oedge, complete

    oedge    = -1
    complete =  0
    IF (PRESENT(ioedge))     oedge=ioedge
    IF (PRESENT(icomplete))  complete=icomplete
    triangle = t_triangle(p=(/ip1,ip2,ip3/), oedge=oedge, complete=complete,  &
      &                   cc=point(-2._wp,-2._wp,-2._wp), r=0._wp, cap_distance=2._wp, rdiscard=-1._wp)
  END FUNCTION triangle


  ! --------------------------------------------------------------------
  !> operator: comparison of triangles (wrt. point indices)
  LOGICAL PURE FUNCTION t_triangle_cmp(a1,a2)
    TYPE(t_triangle), INTENT(IN) :: a1,a2
    t_triangle_cmp =  (a1%p(0) < a2%p(0))                        .OR.                       &
      &        ((a1%p(0) == a2%p(0)) .AND. (a1%p(1) < a2%p(1)))  .OR.                       &
      &        ((a1%p(0) == a2%p(0)) .AND. (a1%p(1) == a2%p(1)) .AND. (a1%p(2) < a2%p(2)))
  END FUNCTION t_triangle_cmp


  ! --------------------------------------------------------------------
  PURE FUNCTION triangle_edge(this, i)
    TYPE(t_edge) :: triangle_edge
    CLASS(t_triangle),  INTENT(IN) :: this
    INTEGER,            INTENT(IN) :: i
    triangle_edge = t_edge( p1=this%p(i), p2=this%p(MOD(i+1,3)) )
  END FUNCTION triangle_edge



  ! --------------------------------------------------------------------
  PURE SUBROUTINE triangle_compute_circumcenter(this, pxyz, subset)
    CLASS(t_triangle),     INTENT(INOUT) :: this
    TYPE(t_point_list),    INTENT(IN)    :: pxyz
    TYPE(t_spherical_cap), INTENT(IN)    :: subset
    ! local variables
    TYPE (t_point) :: a,b,cc0

    ASSOCIATE (p1 => pxyz%a(this%p(0)), p2 => pxyz%a(this%p(1)), p3 => pxyz%a(this%p(2)))
      ! compute circumcenter as (p2-p3) x (p1-p3) normalized:
      a = p2 - p3
      b = p1 - p3
      cc0 = t_point( x=(a%y*b%z - b%y*a%z), y=(a%z*b%x - b%z*a%x), z=(a%x*b%y - a%y*b%x), &
        &            ps=0._wp, gindex=0)
      this%cc = cc0/cc0%norm2()
      ! compute square radius
      this%r     = this%cc*p1                   ! dot product
      this%r     = SQRT(1._wp - this%r*this%r)    ! Pythagoras
      this%cc%ps = this%cc * subset%sorting_direction
      this%cap_distance = subset%point%spherical_dist(this%cc) - this%r
    END ASSOCIATE
  END SUBROUTINE triangle_compute_circumcenter


  ! --------------------------------------------------------------------
  !> cycle triangle indices st. smallest point index is first.
  PURE SUBROUTINE triangle_normalize_indices(this, pts)
    CLASS(t_triangle),     INTENT(INOUT) :: this
    TYPE(t_point_list),    INTENT(IN)    :: pts
    ! local variables
    INTEGER :: tp(0:2), ishift, j

    tp(0:2) = pts%a(this%p(0:2))%gindex
    ishift = 0
    IF      ((tp(1) <= tp(0)) .AND. (tp(1) <= tp(2))) THEN
      ishift = 1
    ELSE IF ((tp(2) <= tp(0)) .AND. (tp(2) <= tp(1))) THEN
      ishift = 2
    END IF
    DO j=0,2
      this%p(j) = tp(MOD(j + ishift, 3))
    END DO
  END SUBROUTINE triangle_normalize_indices


  ! --------------------------------------------------------------------
  ! SUBROUTINES FOR POINT/TRIANGLE LISTS (SORTABLE)
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> exchanges elements a(i), a(j)
  PURE SUBROUTINE initialize_object_list(this)
    CLASS(t_object_list), INTENT(INOUT) :: this
    this%nentries = 0
  END SUBROUTINE initialize_object_list

  ! --------------------------------------------------------------------
  !> constructor for t_point_list type
  FUNCTION point_list(list)
    TYPE(t_point_list) :: point_list
    TYPE(t_point_list), INTENT(IN)  :: list
    CALL point_list%resize(list%nentries)
    point_list%a(0:(list%nentries-1)) = list%a(0:(list%nentries-1))
  END FUNCTION point_list

  ! --------------------------------------------------------------------
  !> constructor for t_triangulation type

  ! --------------------------------------------------------------------
  !> constructor for t_triangulation type


  ! --------------------------------------------------------------------
  !> exchanges elements a(i), a(j)
  PURE SUBROUTINE swap_point_list(this, i,j)
    CLASS(t_point_list), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: i,j
    TYPE(t_point) :: t
    t         = this%a(i)
    this%a(i) = this%a(j)
    this%a(j) = t
  END SUBROUTINE swap_point_list


  ! --------------------------------------------------------------------
  !> @return .TRUE. if a(i) < a(j)
  LOGICAL PURE FUNCTION cmp_point_list(this, i,j)
    CLASS(t_point_list), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: i,j
    cmp_point_list = (this%a(i) < this%a(j))
  END FUNCTION cmp_point_list


  ! --------------------------------------------------------------------
  !> finalization subroutine
  SUBROUTINE destructor_point_list(this)
    CLASS(t_point_list), INTENT(INOUT) :: this
    CHARACTER(*), PARAMETER :: routine = modname//":destructor_point_list"
    INTEGER :: ierrstat
    IF (ALLOCATED(this%a)) THEN
      DEALLOCATE(this%a, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
    this%nentries = 0
  END SUBROUTINE destructor_point_list


  ! --------------------------------------------------------------------
  !> Reserve subroutine.
  !
  !  Requests that the vector capacity be at least enough to contain
  !  @p new_size elements. If new_size is greater than the current
  !  vector capacity, the function causes the array to reallocate its
  !  storage increasing its capacity to new_size.  In all other cases,
  !  the function call does not cause a reallocation and the capacity
  !  is not affected.
  SUBROUTINE reserve_point_list(this, new_size)
    CLASS(t_point_list), INTENT(INOUT) :: this
    INTEGER,             INTENT(IN)    :: new_size
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":reserve_point_list"
    INTEGER :: ierrstat
    TYPE(t_point), ALLOCATABLE :: tmp(:)

    ! allocate temporary storage
    IF (ALLOCATED(this%a)) THEN
      IF (SIZE(this%a) < new_size) THEN
        ALLOCATE(tmp(0:(this%nentries-1)), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        ! copy existing data and resize
        tmp(0:(this%nentries-1)) = this%a(0:(this%nentries-1))
        DEALLOCATE(this%a, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      END IF
    END IF
    IF (.NOT. ALLOCATED(this%a)) THEN
      ALLOCATE(this%a(0:(new_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ELSE IF (SIZE(this%a) < new_size) THEN
      ALLOCATE(this%a(0:(new_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF
    IF (ALLOCATED(tmp)) this%a(0:(this%nentries-1)) = tmp(0:(this%nentries-1))
    ! clean up
    IF (ALLOCATED(tmp)) THEN
      DEALLOCATE(tmp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE reserve_point_list


  ! --------------------------------------------------------------------
  !> Resize subroutine.
  !
  !  Resizes the container so that it contains n elements.  If n is
  !  smaller than the current container size, the content is reduced
  !  to its first n elements.  If n is greater than the current
  !  container size, the content is expanded by inserting at the end
  !  as many elements as needed to reach a size of n.  If n is also
  !  greater than the current container capacity, an automatic
  !  reallocation of the allocated storage space takes place.
  !
  !  This subroutine does not initialize new elements!
  SUBROUTINE resize_point_list(this, new_size)
    CLASS(t_point_list) :: this
    INTEGER,   INTENT(IN) :: new_size
    CALL this%reserve(new_size)
    this%nentries = new_size
  END SUBROUTINE resize_point_list


  ! --------------------------------------------------------------------
  !> "push_back" subroutine, appends element to end of list.
  SUBROUTINE push_back_point_list(this, element)
    CLASS(t_point_list) :: this
    TYPE(t_point), INTENT(IN) :: element
    IF (.NOT. ALLOCATED(this%a)) THEN
      CALL this%reserve(this%nentries + 1)
    ELSE IF (this%nentries == SIZE(this%a)) THEN
      CALL this%reserve(this%nentries + 1)
    END IF
    this%a(this%nentries) = element
    this%nentries         = this%nentries + 1
  END SUBROUTINE push_back_point_list


  ! --------------------------------------------------------------------
  !> exchanges elements a(i), a(j)
  PURE SUBROUTINE swap_triangulation(this, i,j)
    CLASS(t_triangulation), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: i,j
    TYPE(t_triangle) :: t
    t         = this%a(i)
    this%a(i) = this%a(j)
    this%a(j) = t
  END SUBROUTINE swap_triangulation


  ! --------------------------------------------------------------------
  !> @return .TRUE. if a(i) < a(j)
  LOGICAL PURE FUNCTION cmp_triangulation(this, i,j)
    CLASS(t_triangulation), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: i,j
    cmp_triangulation = (this%a(i) < this%a(j))
  END FUNCTION cmp_triangulation


  ! --------------------------------------------------------------------
  !> finalization subroutine
  SUBROUTINE destructor_triangulation(this)
    CLASS(t_triangulation), INTENT(INOUT) :: this
    CHARACTER(*), PARAMETER :: routine = modname//":destructor_triangulation"
    INTEGER :: ierrstat
    IF (ALLOCATED(this%a)) THEN
      DEALLOCATE(this%a, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
    this%nentries = 0
  END SUBROUTINE destructor_triangulation


  ! --------------------------------------------------------------------
  !> Reserve subroutine.
  !
  !  Requests that the vector capacity be at least enough to contain
  !  @p new_size elements. If new_size is greater than the current
  !  vector capacity, the function causes the array to reallocate its
  !  storage increasing its capacity to new_size.  In all other cases,
  !  the function call does not cause a reallocation and the capacity
  !  is not affected.
  SUBROUTINE reserve_triangulation(this, new_size)
    CLASS(t_triangulation), INTENT(INOUT) :: this
    INTEGER,                INTENT(IN)    :: new_size
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":reserve_triangulation"
    INTEGER :: ierrstat
    TYPE(t_triangle), ALLOCATABLE :: tmp(:)

    ! allocate temporary storage
    IF (ALLOCATED(this%a)) THEN
      IF (SIZE(this%a) < new_size) THEN
        ALLOCATE(tmp(0:(this%nentries-1)), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        ! copy existing data and resize
        tmp(0:(this%nentries-1)) = this%a(0:(this%nentries-1))
        DEALLOCATE(this%a, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      END IF
    END IF
    IF (.NOT. ALLOCATED(this%a)) THEN
      ALLOCATE(this%a(0:(new_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ELSE IF (SIZE(this%a) < new_size) THEN
      ALLOCATE(this%a(0:(new_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF
    IF (ALLOCATED(tmp)) this%a(0:(this%nentries-1)) = tmp(0:(this%nentries-1))
    ! clean up
    IF (ALLOCATED(tmp)) THEN
      DEALLOCATE(tmp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE reserve_triangulation


  ! --------------------------------------------------------------------
  !> Resize subroutine.
  !
  !  Resizes the container so that it contains n elements.  If n is
  !  smaller than the current container size, the content is reduced
  !  to its first n elements.  If n is greater than the current
  !  container size, the content is expanded by inserting at the end
  !  as many elements as needed to reach a size of n.  If n is also
  !  greater than the current container capacity, an automatic
  !  reallocation of the allocated storage space takes place.
  !
  !  This subroutine does not initialize new elements!
  SUBROUTINE resize_triangulation(this, new_size)
    CLASS(t_triangulation) :: this
    INTEGER,   INTENT(IN) :: new_size
    CALL this%reserve(new_size)
    this%nentries = new_size
  END SUBROUTINE resize_triangulation


  ! --------------------------------------------------------------------
  !> "push_back" subroutine, appends element to end of list.
  SUBROUTINE push_back_triangulation(this, element)
    CLASS(t_triangulation) :: this
    TYPE(t_triangle), INTENT(IN) :: element
    IF (.NOT. ALLOCATED(this%a)) THEN
      CALL this%reserve(this%nentries + 1)
    ELSE IF (this%nentries == SIZE(this%a)) THEN
      CALL this%reserve(this%nentries + 1)
    END IF

    this%a(this%nentries) = element
    this%nentries         = this%nentries + 1
  END SUBROUTINE push_back_triangulation


  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  !
  !  Ordering after the sorting process: smallest...largest.
  !
  !  The user provides two callback functions: "cmp"("<") and "swap"
  !  which emulate a "Fortran template".
  !
  RECURSIVE SUBROUTINE quicksort_sortable_list(this, in_l, in_r)
    CLASS(t_sortable_list) :: this
    INTEGER, INTENT(IN), OPTIONAL :: in_l,in_r     !< left, right partition indices
    ! local variables
    INTEGER :: i,j,m,l,r

    l = 0
    r = this%nentries-1
    IF (PRESENT(in_l))  l=in_l
    IF (PRESENT(in_r))  r=in_r

    IF (r>l) THEN
      ! median-of-three selection of partitioning element
      IF ((r-l) > 3) THEN
        m = (l+r)/2
        IF (this%cmp(m,l)) CALL this%swap(l,m)
        IF (this%cmp(r,l)) THEN
          CALL this%swap(l,r)
        ELSE IF (this%cmp(m,r)) THEN
          CALL this%swap(r,m)
        END IF
      END IF
      i = l-1
      j = r
      LOOP : DO
        CNTLOOP1 : DO
          i = i+1
          IF (.NOT. this%cmp(i,r)) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2 : DO
          j = j-1
          IF (.NOT. this%cmp(r,j) .OR. (j==0)) EXIT CNTLOOP2
        END DO CNTLOOP2
        CALL this%swap(i,j)
        IF (j <= i) EXIT LOOP
      END DO LOOP
      CALL this%swap(i,j) ! undo extra exchange
      CALL this%swap(i,r) ! put partitioning element into position
      CALL this%quicksort(l,i-1)
      CALL this%quicksort(i+1,r)
    END IF
  END SUBROUTINE quicksort_sortable_list


  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  !
  !  Ordering after the sorting process: smallest...largest.
  !
  RECURSIVE SUBROUTINE quicksort_point_list(this, in_l, in_r)
    CLASS(t_point_list) :: this
    INTEGER, INTENT(IN), OPTIONAL :: in_l,in_r     !< left, right partition indices
    ! local variables
    INTEGER :: i,j,m,l,r

    l = 0
    r = this%nentries-1
    IF (PRESENT(in_l))  l=in_l
    IF (PRESENT(in_r))  r=in_r

    IF (r>l) THEN
      ! median-of-three selection of partitioning element
      IF ((r-l) > 3) THEN
        m = (l+r)/2

        IF (this%a(m) < this%a(l)) CALL swap_point_list(this,l,m)
        IF (this%a(r) < this%a(l)) THEN
          CALL swap_point_list(this,l,r)
        ELSE IF (this%a(m) < this%a(r)) THEN
          CALL swap_point_list(this,r,m)
        END IF
      END IF
      i = l-1
      j = r
      LOOP : DO
        CNTLOOP1 : DO
          i = i+1
          IF (.NOT. (this%a(i) < this%a(r))) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2 : DO
          j = j-1
          IF (.NOT. (this%a(r) < this%a(j)) .OR. (j==0)) EXIT CNTLOOP2
        END DO CNTLOOP2
        CALL swap_point_list(this,i,j)
        IF (j <= i) EXIT LOOP
      END DO LOOP
      CALL swap_point_list(this,i,j) ! undo extra exchange
      CALL swap_point_list(this,i,r) ! put partitioning element into position
      CALL this%quicksort(l,i-1)
      CALL this%quicksort(i+1,r)
    END IF
  END SUBROUTINE quicksort_point_list

  ! DART infrastructure adapter; no MPI or file I/O in this module.
  SUBROUTINE finish(routine, message)
    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CALL error_handler(E_ERR, routine, message)
  END SUBROUTINE finish
END MODULE icon_delaunay_types_mod
