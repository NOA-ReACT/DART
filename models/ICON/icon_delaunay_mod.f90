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

! Parallel Bowyer/Watson algorithm for the construction of a Delaunay
! triangulation on the sphere.
!
! This algorithm borrows from
!
! (1) Paul Bourke: "Efficient Triangulation Algorithm Suitable for
!      Terrain Modelling or An Algorithm for Interpolating
!      Irregularly-Spaced Data with Applications in Terrain
!      Modelling." Pan Pacific Computer Conference, Beijing, China,
!      1989.
!
! (2) Jacobsen, D. W.; Gunzburger, M.; Ringler, T.; Burkardt, J.;
!      Peterson, J.: "Parallel algorithms for planar and spherical
!      Delaunay construction with an application to centroidal
!      Voronoi tessellations."  Geoscientific Model Development;2013,
!      Vol. 6 Issue 4, p.1353 (global Delaunay condition for
!      parallelization)
!
! (3) Renka, R. J. "Interpolation of Data on the Surface of a
!      Sphere." ACM Trans. Math. Softw., ACM, 1984, 10, 417-436
!
! (4) J. Shewchuck: "Lecture Notes on Delaunay Triangulation."
!      Univ. of California, Berkeley. 2012. (ghost point instead of
!      "super-triangle", as described in Section 3.4 of)
!
! with modifications:
! - sorting of points wrt. spherical cap (point subset for parallelization)
! - multi-threading and indexing shortcuts
! - global Delaunay condition is used as abort criterion for triangulation

! Adapted from ICON 8597da45ef4b86323f3fb844caedc4ae5e1ffc01.
! Serial execution of the selected repartition kernel, NOT ICON's different
! triangulate() routine. OpenMP/MPI and NEC branches are not included.
MODULE icon_delaunay_mod
  USE types_mod, ONLY: wp => r8
  USE utilities_mod, ONLY: error_handler, E_ERR
  USE icon_delaunay_types_mod, ONLY: t_edge, t_point, t_point_list, t_triangle, &
       t_triangulation, t_spherical_cap, point_list, triangle, point, &
       OPERATOR(/), OPERATOR(*), OPERATOR(==), ccw_spherical, circum_circle_spherical
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: triangulate_serial
CONTAINS
  SUBROUTINE triangulate_serial(points, tri, subset, ignore_completeness)

    TYPE(t_point_list),    INTENT(IN)            :: points
    TYPE(t_triangulation), INTENT(INOUT), TARGET :: tri
    TYPE(t_spherical_cap), INTENT(IN)            :: subset
    LOGICAL,               INTENT(IN)            :: ignore_completeness
    ! local variables

    !> triangle count that triggers clean-up of "tri" data structure
    INTEGER, PARAMETER :: cleanup_limit  = 50000
    INTEGER, PARAMETER :: MAX_EDGES      = 10000   !< size of edge buffer

    ! For the triangles that are created from the edges 0,1,2: compute
    ! the edge opposite to "ghost point" (if "ghost point" part of this
    ! triangle tri[j]). This lookup table is best figured out on a
    ! sheet of paper...
    INTEGER, PARAMETER :: next_oedge(0:2,-1:2) = &
      &  RESHAPE( (/ -1, -1, -1, &
      &              -1,  2,  1, &
      &               1, -1,  2, &
      &               2,  1, -1 /), (/3,4/) )
    TYPE(t_point_list)        :: pxyz
    TYPE(t_point)             :: ipoint
    TYPE(t_edge)              :: bedge
    TYPE(t_triangle), POINTER :: jtri
    INTEGER                   :: i,j, npts, gp, ninvalid, j0, ndiscard, ipt, ne, ntri, &
      &                          new_ninvalid, new_ndiscard, jmin,            &
      &                          oedge, jmin0, this_edge,k, endtri
    LOGICAL                   :: inside
    REAL(wp)                  :: cos_radius, rlimit

    ! edges buffer for re-Delaunay step (together with
    ! "edge-index-opposite-to-ghost-point")
    TYPE(t_edge)              :: edges(0:(MAX_EDGES-1))
    INTEGER                   :: edges_oedge(0:(MAX_EDGES-1))
    LOGICAL                   :: edges_valid(0:(MAX_EDGES-1))

    INTEGER                   :: jmin_t


    ! copy all points since they will be re-ordered:
    npts = points%nentries
    CALL pxyz%initialize()
    pxyz = point_list(points)
    DO i=0,(npts-1)
      pxyz%a(i)        = pxyz%a(i)/pxyz%a(i)%norm2()
      pxyz%a(i)%ps     = pxyz%a(i) * subset%sorting_direction
    END DO
    CALL pxyz%quicksort()
    ! assume triangulation covers entire sphere (data structure grows
    ! dynamically if not)
    CALL tri%reserve(2*npts-4 + cleanup_limit)
    tri%nentries = 0

    ! set up the initial triangle: insert the first three points
    IF (ccw_spherical(pxyz%a(0), pxyz%a(1), pxyz%a(2))<0) THEN
      CALL tri%push_back(triangle(0, 1, 2)) ! clockwise starting triangle
    ELSE
      CALL tri%push_back(triangle(0, 2, 1)) ! counter-clockwise starting triangle
    END IF

    ! connect triangle edges to a "ghost point"
    gp = npts
    CALL pxyz%push_back(point(-1._wp,-1._wp,-1._wp))
    ! add three initial "ghost triangles"
    CALL tri%push_back(triangle(tri%a(0)%p(1),tri%a(0)%p(0),gp, 0))
    CALL tri%push_back(triangle(tri%a(0)%p(2),tri%a(0)%p(1),gp, 0))
    CALL tri%push_back(triangle(tri%a(0)%p(0),tri%a(0)%p(2),gp, 0))
    CALL tri%a(0)%compute_circumcenter(pxyz, subset)

    ninvalid = 0       ! no. of invalid triangles in "tri" DATA structure
    j0       = 0       ! first triangle index that is not "complete"
    ndiscard = 0       ! no. of triangles that violate global Delaunay condition

    ipt      = 3           ! current point index
    ipoint   = pxyz%a(ipt) ! current point
    IF (subset%radius < 0._wp) THEN
      cos_radius = 99._wp
    ELSE
      cos_radius = -1._wp*COS(subset%radius)
    END IF

    ! Include each point one at a time into the existing mesh. Abort if
    ! all points have been processed or the requested radius has been
    ! reached without violating the global Delaunay condition.
    LOOP : DO ipt=3,(npts-1)
      !      IF (dbg_level > 10) THEN
      !        IF (MOD(ipt,1000) == 0) THEN
      !          WRITE (0,*) "ipt = ", ipt
      !        END IF
      !      END IF
      IF (ipoint%ps > cos_radius) THEN
        IF ((ndiscard==0) .OR. ignore_completeness)  EXIT LOOP
      END IF

      ipoint = pxyz%a(ipt)
      ! Set up the edge buffer.
      !
      ! If the point pxyz[i] lies inside the circumcircle then the
      ! three edges of that triangle are added to the edge buffer and
      ! that triangle is removed.

      ne           = 0 ! no. of entries in edge buffer
      ntri         = tri%nentries
      jmin_t       = ntri
      new_ninvalid = 0
      new_ndiscard = 0
      jmin    = ntri

      DO j=j0,(ntri-1)
        jtri => tri%a(j)

        IF (jtri%complete == 0) THEN
          ! test if triangle remains on "discard" list
          IF (-1._wp*ipoint%ps < jtri%rdiscard) THEN
            new_ndiscard  = new_ndiscard-1
            jtri%rdiscard = -1._wp
          END IF

          oedge = jtri%oedge
          jmin0 = jmin
          IF (jmin == ntri)  jmin = j
          IF (oedge /= -1) THEN
            ! boundary edge
            bedge  = jtri%edge(oedge)
            inside = (ccw_spherical(pxyz%a(bedge%p1), pxyz%a(bedge%p2), ipoint)<0)
          ELSE
            ! interior edge
            inside = circum_circle_spherical(ipoint, pxyz, jtri%p)
            IF (.NOT. inside .AND. ((ipoint%ps - jtri%cc%ps) > jtri%r)) THEN
              jtri%complete =     1
              jmin          = jmin0
            END IF
          END IF

          IF (inside) THEN
            ! push triangle's edges onto edge list
            IF (ne+3 > MAX_EDGES) CALL finish("triangulate_serial", "ICON edge buffer exceeded")
            this_edge = ne
            ne = ne + 3
            DO k=0,2
              edges(this_edge)       = jtri%edge(k)
              edges_oedge(this_edge) = next_oedge(k,oedge)
              edges_valid(this_edge) = .TRUE.
              this_edge = this_edge + 1
            END DO

            jtri%complete =       2 ! mark triangle for removal
            new_ninvalid = new_ninvalid + 1
            IF (jtri%rdiscard > -1._wp)  new_ndiscard=new_ndiscard-1
          END IF
        END IF
      END DO
      IF (jmin >= j0) THEN
        jmin_t = jmin
      ELSE
        jmin_t = ntri
      END IF

      ! remove multiple edges
      DO j=0,(ne-1)
        IF (edges_valid(j)) THEN
          DO k=j+1,(ne-1)
            IF (edges(j) == edges(k)) THEN
              edges_valid(j) = .FALSE.
              edges_valid(k) = .FALSE.
            END IF
          END DO
        END IF
      END DO

      ninvalid = ninvalid + new_ninvalid
      IF (jmin_t /= ntri)  j0 = jmin_t

      ! form new triangles for the current point. All edges are
      ! arranged in clockwise order.
      DO j=0,(ne-1)
        IF (edges_valid(j))  CALL tri%push_back( triangle(edges(j)%p1, edges(j)%p2, ipt, edges_oedge(j)) )
      END DO
      ndiscard     = ndiscard + new_ndiscard
      new_ndiscard = 0

      ! compute circumcircle for triangle without ghost point
      DO j=ntri,(tri%nentries-1)
        IF (tri%a(j)%oedge == -1)  CALL tri%a(j)%compute_circumcenter(pxyz, subset)
      END DO
      ! mark triangle if it violates the global Delaunay condition:
      IF (cos_radius > ipoint%ps) THEN
        DO j=ntri,(tri%nentries-1)
          IF ((tri%a(j)%oedge == -1) .AND. (tri%a(j)%cap_distance < subset%radius)) THEN
            new_ndiscard  = new_ndiscard + 1
            tri%a(j)%rdiscard = COS(tri%a(j)%cap_distance+2._wp*tri%a(j)%r)
          END IF
        END DO
      END IF
      ndiscard = ndiscard + new_ndiscard

      ! remove "invalid" triangles at regular intervals:
      IF (ninvalid > cleanup_limit) THEN
        ! first, replace invalid triangles by "complete" ones
        ntri   = tri%nentries
        endtri = ntri - 1
        LOOP3 : DO j=0,(ntri-1)
          IF (tri%a(j)%complete == 2) THEN
	    ! find the next "complete" triangle from the back
            LOOP4 : DO
              IF ((tri%a(endtri)%complete==1) .OR. (endtri==j))  EXIT LOOP4
              endtri = endtri - 1
            END DO LOOP4
            IF (endtri /= j) THEN
              tri%a(j)               = tri%a(endtri)
              tri%a(endtri)%complete = 2
            ELSE
              IF (j<j0)  j0=j
              EXIT LOOP3 ! only incomplete triangles follow
            END IF
          END IF
        END DO LOOP3
        ! if there are no "complete" triangles left at the list end,
        ! simply remove the invalid triangles:
        endtri = tri%nentries - 1
        LOOP5 : DO
          IF (j==endtri)  EXIT LOOP5
          IF (tri%a(j)%complete == 2) THEN
            tri%a(j) = tri%a(endtri)
            j      = j      - 1
            endtri = endtri - 1
            ntri   = ntri   - 1
          END IF
          j = j + 1
        END DO LOOP5
        CALL tri%resize(ntri)
        ninvalid = 0
      END IF
    END DO LOOP
    ! compute limit radius (last point inserted into triangulation)
    IF (subset%radius < 0._wp) THEN
      rlimit = 99._wp
    ELSE
      rlimit = subset%point%spherical_dist(ipoint)
    END IF
    ! discard triangles:
    ntri = tri%nentries

    j = 0
    LOOP6 : DO
      IF (j>=ntri)  EXIT LOOP6
      ! remove all ghost point triangles or triangles that will not
      ! hold the global Delaunay condition [Jacobson2013]:
      IF ((tri%a(j)%complete==2) .OR. &
        & (tri%a(j)%oedge /= -1) .OR. &
        & ((tri%a(j)%cap_distance+2._wp*tri%a(j)%r) > rlimit)) THEN
        ntri     = ntri - 1
        tri%a(j) = tri%a(ntri)
      ELSE
        CALL tri%a(j)%normalize_indices(pxyz)
        j = j + 1
      END IF
    END DO LOOP6
    CALL tri%resize(ntri)
    CALL pxyz%destructor()
  END SUBROUTINE triangulate_serial
  SUBROUTINE finish(routine, message)
    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CALL error_handler(E_ERR, routine, message)
  END SUBROUTINE finish
END MODULE icon_delaunay_mod
