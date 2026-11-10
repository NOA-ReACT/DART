! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

!> DART interface to the ICON model (icosahedral triangular grid, height based
!> terrain following vertical coordinate).
!>
!> The state is read from ICON output files on the native grid, i.e. variables with
!> the netCDF dimensions (time, height, ncells), (time, height_2, ncells) or
!> (time, ncells). All state variables live on cell centres. The horizontal winds must
!> be the diagnostic "u" and "v" on cells, not the prognostic "vn" on edges: ICON
!> itself turns the analysed u and v into vn when it starts from an analysis
!> (initicon_nml: init_mode = 1 or 5).
!>
!> Three files are needed:
!>   template_file  an ICON output file with the state variables
!>   grid_file      the ICON grid file (cell and vertex positions, connectivity)
!>   height_file    an ICON output file with the half level heights "z_ifc"
!>
!> Model levels are counted from the top in ICON (level 1 is the highest one).
!> Locations with a vertical coordinate of VERTISLEVEL count from the bottom, as the
!> forward operators that use them expect.

module model_mod

use        types_mod, only : r4, r8, i8, MISSING_R8, vtablenamelength, DEG2RAD, RAD2DEG

use time_manager_mod, only : time_type, set_time, set_date, set_calendar_type, &
                             GREGORIAN, operator(+)

use     location_mod, only : location_type, get_close_type,                    &
                             loc_get_close_obs => get_close_obs,               &
                             loc_get_close_state => get_close_state,           &
                             set_location, get_location, query_location,       &
                             VERTISUNDEF, VERTISSURFACE, VERTISLEVEL,          &
                             VERTISPRESSURE, VERTISHEIGHT,                     &
                             set_vertical_localization_coord

use xyz_location_mod, only : xyz_location_type, xyz_set_location,              &
                             xyz_get_close_type, xyz_get_close_init,           &
                             xyz_get_close_destroy, xyz_find_nearest

use    utilities_mod, only : error_handler, E_ERR, E_MSG,                      &
                             nmlfileunit, do_nml_file, do_nml_term, do_output, &
                             find_namelist_in_file, check_namelist_read

use     obs_kind_mod, only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT,       &
                             QTY_TEMPERATURE, QTY_PRESSURE,                    &
                             QTY_SURFACE_PRESSURE, QTY_SPECIFIC_HUMIDITY,      &
                             QTY_VAPOR_MIXING_RATIO, QTY_DENSITY,              &
                             QTY_POTENTIAL_TEMPERATURE,                        &
                             QTY_GEOPOTENTIAL_HEIGHT, QTY_GEOMETRIC_HEIGHT,    &
                             QTY_SURFACE_ELEVATION, QTY_2M_TEMPERATURE,        &
                             QTY_2M_SPECIFIC_HUMIDITY,                         &
                             QTY_10M_U_WIND_COMPONENT, QTY_10M_V_WIND_COMPONENT

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time,                  &
                                 nc_begin_define_mode, nc_end_define_mode,     &
                                 nc_open_file_readonly, nc_close_file, nc_check

use distributed_state_mod, only : get_state

use state_structure_mod, only : add_domain, get_dart_vector_index,             &
                                get_domain_size, get_model_variable_indices,   &
                                get_varid_from_kind, get_num_dims,             &
                                get_dim_length

use ensemble_manager_mod, only : ensemble_type

use default_model_mod, only : pert_model_copies, write_model_time,             &
                              init_time => fail_init_time,                     &
                              init_conditions => fail_init_conditions,         &
                              convert_vertical_state, adv_1step,               &
                              parse_variables_clamp,                           &
                              MAX_STATE_VARIABLE_FIELDS_CLAMP

use netcdf

implicit none
private

! routines required by DART code - will be called from filter and other DART executables.
public :: get_model_size,                      &
          get_state_meta_data,                 &
          model_interpolate,                   &
          end_model,                           &
          static_init_model,                   &
          nc_write_model_atts,                 &
          get_close_obs,                       &
          get_close_state,                     &
          pert_model_copies,                   &
          convert_vertical_obs,                &
          convert_vertical_state,              &
          read_model_time,                     &
          adv_1step,                           &
          init_time,                           &
          init_conditions,                     &
          shortest_time_between_assimilations, &
          write_model_time

character(len=*), parameter :: source = 'icon/model_mod.f90'

! Return codes of model_interpolate
integer, parameter :: ABOVE_TOP_LEVEL        = 2   ! VERTISLEVEL at or above the top; the
                                                    ! column operators rely on this value
integer, parameter :: OUTSIDE_DOMAIN         = 11
integer, parameter :: IN_BOUNDARY_ZONE       = 12
integer, parameter :: BELOW_SURFACE          = 15
integer, parameter :: ABOVE_MODEL_TOP        = 16
integer, parameter :: SURFACE_TOO_DIFFERENT  = 17
integer, parameter :: BAD_VERTICAL_FOR_FIELD = 18
integer, parameter :: NO_PRESSURE_IN_STATE   = 19
integer, parameter :: UNKNOWN_QUANTITY       = 20
integer, parameter :: UNKNOWN_VERTICAL       = 21
integer, parameter :: BAD_VALUE              = 22

! Physical constants, as in ICON (mo_physical_constants)
real(r8), parameter :: rd     = 287.04_r8       ! gas constant of dry air [J/K/kg]
real(r8), parameter :: cpd    = 1004.64_r8      ! heat capacity of dry air [J/K/kg]
real(r8), parameter :: vtmpc1 = 461.51_r8 / rd - 1.0_r8
real(r8), parameter :: p0ref  = 100000.0_r8     ! reference pressure [Pa]
real(r8), parameter :: earth_radius = 6371229.0_r8

character(len=512) :: string1, string2

logical :: module_initialized = .false.
integer :: dom_id
type(time_type) :: assimilation_time_step

! Grid, read once and kept by every task
integer :: ncells = 0
integer :: nverts = 0
integer :: nlev   = 0     ! number of full levels
integer :: nlevp1 = 0     ! number of half levels
real(r8), allocatable :: clon(:), clat(:)       ! cell centres [degrees], lon in [0, 360)
real(r8), allocatable :: cxyz(:,:)              ! (3, ncells) unit vectors of cell centres
real(r8), allocatable :: vxyz(:,:)              ! (3, nverts) unit vectors of vertices
integer,  allocatable :: vertex_of_cell(:,:)    ! (3, ncells)
integer,  allocatable :: cells_of_vertex(:,:)   ! (6, nverts) in order around the vertex
integer,  allocatable :: ncells_of_vertex(:)    ! how many of the 6 exist
integer,  allocatable :: refin_c_ctrl(:)        ! distance to the lateral boundary in cell rows
real(r4), allocatable :: z_mc(:,:)              ! (nlev, ncells) height of full levels [m]
real(r4), allocatable :: z_ifc(:,:)             ! (nlevp1, ncells) height of half levels [m]

! To find the cell centre that is closest to a location
type(xyz_get_close_type) :: cell_search
type(xyz_location_type), allocatable :: cell_locs(:)

! Namelist
character(len=256) :: template_file = 'icon_template.nc'
character(len=256) :: grid_file     = 'icon_grid.nc'
character(len=256) :: height_file   = 'icon_heights.nc'
integer  :: assimilation_period_days    = 0
integer  :: assimilation_period_seconds = 3600
!> Observations whose closest cell is within this many cell rows of the lateral
!> boundary of a limited area grid are rejected. 0 keeps all of them.
integer  :: boundary_width = 4
!> Surface observations are rejected if the model surface is more than this many
!> metres away from the height of the observation. Negative: no check.
real(r8) :: sfc_elev_max_diff = -1.0_r8
integer  :: debug = 0
character(len=vtablenamelength) :: variables(MAX_STATE_VARIABLE_FIELDS_CLAMP) = ' '

namelist /model_nml/ template_file, grid_file, height_file,         &
                     assimilation_period_days,                      &
                     assimilation_period_seconds,                   &
                     boundary_width, sfc_elev_max_diff, debug, variables

contains

!------------------------------------------------------------------
!> One time initialisation: namelist, state layout, grid.

subroutine static_init_model()

integer :: iunit, io

if (module_initialized) return
module_initialized = .true.

call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type(GREGORIAN)

assimilation_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

dom_id = add_domain(template_file, parse_variables_clamp(variables))

call read_grid()
call read_heights()
call init_cell_search()

! Every state element has a height, so that is what the vertical localisation uses
call set_vertical_localization_coord(VERTISHEIGHT)

end subroutine static_init_model

!------------------------------------------------------------------

function get_model_size()

integer(i8) :: get_model_size

if (.not. module_initialized) call static_init_model

get_model_size = get_domain_size(dom_id)

end function get_model_size

!------------------------------------------------------------------

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if (.not. module_initialized) call static_init_model

shortest_time_between_assimilations = assimilation_time_step

end function shortest_time_between_assimilations

!------------------------------------------------------------------
!> Read the horizontal grid. Connectivity arrays are (3, cell) or (6, vertex) in the
!> netCDF file, which the Fortran interface sees the other way round.

subroutine read_grid()

character(len=*), parameter :: routine = 'read_grid'

integer :: ncid, i, k
real(r8), allocatable :: lon(:), lat(:)
integer,  allocatable :: itmp(:,:)

ncid = nc_open_file_readonly(grid_file, routine)

ncells = dimension_size(ncid, 'cell', grid_file)
nverts = dimension_size(ncid, 'vertex', grid_file)

allocate(clon(ncells), clat(ncells), cxyz(3, ncells))
call get_real_1d(ncid, 'clon', clon, grid_file)
call get_real_1d(ncid, 'clat', clat, grid_file)
do i = 1, ncells
   cxyz(:, i) = unit_vector(clon(i), clat(i))
enddo
clon = clon * RAD2DEG
clat = clat * RAD2DEG
where (clon < 0.0_r8)    clon = clon + 360.0_r8
where (clon >= 360.0_r8) clon = clon - 360.0_r8

allocate(lon(nverts), lat(nverts), vxyz(3, nverts))
call get_real_1d(ncid, 'vlon', lon, grid_file)
call get_real_1d(ncid, 'vlat', lat, grid_file)
do i = 1, nverts
   vxyz(:, i) = unit_vector(lon(i), lat(i))
enddo
deallocate(lon, lat)

allocate(itmp(ncells, 3), vertex_of_cell(3, ncells))
call get_int_2d(ncid, 'vertex_of_cell', itmp, grid_file)
do k = 1, 3
   vertex_of_cell(k, :) = itmp(:, k)
enddo
deallocate(itmp)

allocate(itmp(nverts, 6), cells_of_vertex(6, nverts), ncells_of_vertex(nverts))
call get_int_2d(ncid, 'cells_of_vertex', itmp, grid_file)
do i = 1, nverts
   call order_cells_around_vertex(i, itmp(i, :))
enddo
deallocate(itmp)

allocate(refin_c_ctrl(ncells))
if (nf90_inq_varid(ncid, 'refin_c_ctrl', k) == NF90_NOERR) then
   call nc_check(nf90_get_var(ncid, k, refin_c_ctrl), routine, 'get refin_c_ctrl')
else
   refin_c_ctrl = 0
endif

call nc_close_file(ncid, routine)

if (do_output()) then
   write(string1, *) 'ICON grid with ', ncells, ' cells read from ', trim(grid_file)
   call error_handler(E_MSG, routine, string1, source)
endif

end subroutine read_grid

!------------------------------------------------------------------
!> Store the cells that share vertex iv in the order they are met when walking
!> around it. The centres of these cells are the corners of a polygon (a hexagon,
!> in most places) around the vertex. These polygons cover the grid without gaps
!> or overlaps, and each is cut into triangles that fan out from its first cell.

subroutine order_cells_around_vertex(iv, cells)

integer, intent(in) :: iv
integer, intent(in) :: cells(6)

integer  :: i, j, n, itmp
real(r8) :: angle(6), east(3), north(3), d(3), rtmp
integer  :: list(6)

n = 0
do i = 1, 6
   if (cells(i) >= 1 .and. cells(i) <= ncells) then
      n = n + 1
      list(n) = cells(i)
   endif
enddo

! Two directions at right angles to each other and to the vertex position
if (abs(vxyz(3, iv)) < 0.9_r8) then
   east = cross(vxyz(:, iv), (/ 0.0_r8, 0.0_r8, 1.0_r8 /))
else
   east = cross(vxyz(:, iv), (/ 1.0_r8, 0.0_r8, 0.0_r8 /))
endif
east  = east / sqrt(sum(east**2))
north = cross(east, vxyz(:, iv))

do i = 1, n
   d = cxyz(:, list(i)) - vxyz(:, iv)
   angle(i) = atan2(sum(d * north), sum(d * east))
enddo

! Insertion sort, there are at most 6
do i = 2, n
   j = i
   do while (j > 1)
      if (angle(j - 1) <= angle(j)) exit
      rtmp = angle(j);  angle(j) = angle(j - 1);  angle(j - 1) = rtmp
      itmp = list(j);   list(j)  = list(j - 1);   list(j - 1)  = itmp
      j = j - 1
   enddo
enddo

! Start with the cell that has the smallest index. The triangles are fanned out
! from the first cell, and this way the fan is the same one no matter which of the
! cells an observation is closest to.
if (n > 0) list(1:n) = cshift(list(1:n), minloc(list(1:n), 1) - 1)

cells_of_vertex(:, iv) = 0
cells_of_vertex(1:n, iv) = list(1:n)
ncells_of_vertex(iv) = n

end subroutine order_cells_around_vertex

!------------------------------------------------------------------
!> Read the half level heights, which do not change during a run.

subroutine read_heights()

character(len=*), parameter :: routine = 'read_heights'

integer :: ncid, varid, ndims, dimids(NF90_MAX_VAR_DIMS), dimlen(3), i, k
real(r4), allocatable :: tmp(:,:)

ncid = nc_open_file_readonly(height_file, routine)

call nc_check(nf90_inq_varid(ncid, 'z_ifc', varid), routine, &
              'no z_ifc in '//trim(height_file))
call nc_check(nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), &
              routine, 'inquire z_ifc')
if (ndims < 2 .or. ndims > 3) then
   call error_handler(E_ERR, routine, 'z_ifc must be (height, ncells) with or without time', &
                      source, text2=trim(height_file))
endif
dimlen = 1
do i = 1, ndims
   call nc_check(nf90_inquire_dimension(ncid, dimids(i), len=dimlen(i)), routine, 'dimension')
enddo
if (dimlen(1) /= ncells) then
   write(string1, *) 'z_ifc has ', dimlen(1), ' cells, the grid file has ', ncells
   call error_handler(E_ERR, routine, string1, source, text2=trim(height_file))
endif

nlevp1 = dimlen(2)
nlev   = nlevp1 - 1

allocate(tmp(ncells, nlevp1), z_ifc(nlevp1, ncells), z_mc(nlev, ncells))
! The first record, if there is a time dimension
call nc_check(nf90_get_var(ncid, varid, tmp, start=(/ 1, 1, 1 /), &
              count=(/ ncells, nlevp1, 1 /)), routine, 'get z_ifc')
do k = 1, nlevp1
   z_ifc(k, :) = tmp(:, k)
enddo
deallocate(tmp)
do k = 1, nlev
   z_mc(k, :) = 0.5_r4 * (z_ifc(k, :) + z_ifc(k + 1, :))
enddo

call nc_close_file(ncid, routine)

if (z_ifc(1, 1) < z_ifc(nlevp1, 1)) then
   call error_handler(E_ERR, routine, 'z_ifc must run from the model top downwards', &
                      source, text2=trim(height_file))
endif

end subroutine read_heights

!------------------------------------------------------------------

subroutine init_cell_search()

integer  :: i, k, iv, other
real(r8) :: maxdist

allocate(cell_locs(ncells))
do i = 1, ncells
   cell_locs(i) = xyz_set_location(clon(i), clat(i), 0.0_r8, earth_radius)
enddo

! Size of the search boxes: the largest distance between cells that share a vertex
maxdist = 0.0_r8
do i = 1, ncells, max(1, ncells / 20000)
   do k = 1, 3
      iv = vertex_of_cell(k, i)
      do other = 1, ncells_of_vertex(iv)
         maxdist = max(maxdist, sqrt(sum((cxyz(:, i) - cxyz(:, cells_of_vertex(other, iv)))**2)))
      enddo
   enddo
enddo

call xyz_get_close_init(cell_search, maxdist * earth_radius, ncells, cell_locs)

end subroutine init_cell_search

!------------------------------------------------------------------
!> Horizontal part of the interpolation: up to three cells around the location and
!> their weights (barycentric, in a triangle of cell centres). The weights are
!> never negative and add up to one, so the result stays between the values of the
!> three cells: no overshoots, and positive fields stay positive.

subroutine find_triangle(lon, lat, nc, cells, wts, ier)

real(r8), intent(in)  :: lon, lat     ! degrees
integer,  intent(out) :: nc
integer,  intent(out) :: cells(3)
real(r8), intent(out) :: wts(3)
integer,  intent(out) :: ier

type(xyz_location_type) :: point
integer  :: c, rc, k, iv, n, m, a, b, first
real(r8) :: p(3), w(3), total
real(r8), parameter :: tolerance = -1.0e-9_r8

nc    = 0
cells = 0
wts   = 0.0_r8
ier   = 0

point = xyz_set_location(lon, lat, 0.0_r8, earth_radius)
call xyz_find_nearest(cell_search, point, cell_locs, c, rc)
if (rc /= 0 .or. c < 1) then
   ier = OUTSIDE_DOMAIN
   return
endif

if (boundary_width > 0) then
   if (refin_c_ctrl(c) >= 1 .and. refin_c_ctrl(c) <= boundary_width) then
      ier = IN_BOUNDARY_ZONE
      return
   endif
endif

p = unit_vector(lon * DEG2RAD, lat * DEG2RAD)

! The location is in the polygon of cell centres around one of the vertices of
! cell c. Each polygon has one fixed set of triangles, so two locations next to each
! other get (nearly) the same value even if their nearest cells differ: the
! interpolated field is continuous.
do k = 1, 3
   iv = vertex_of_cell(k, c)
   n  = ncells_of_vertex(iv)
   if (n < 3) cycle
   first = cells_of_vertex(1, iv)

   do m = 2, n - 1
      a = cells_of_vertex(m, iv)
      b = cells_of_vertex(m + 1, iv)
      w(1) = triple(p, cxyz(:, a), cxyz(:, b))
      w(2) = triple(cxyz(:, first), p, cxyz(:, b))
      w(3) = triple(cxyz(:, first), cxyz(:, a), p)
      total = sum(w)
      if (abs(total) < tiny(1.0_r8)) cycle
      w = w / total
      if (all(w >= tolerance)) then
         nc    = 3
         cells = (/ first, a, b /)
         wts   = max(w, 0.0_r8)
         wts   = wts / sum(wts)
         return
      endif
   enddo
enddo

! On the rim of a limited area grid a location can be inside cell c without being
! between three cell centres
w(1) = triple(p, vxyz(:, vertex_of_cell(2, c)), vxyz(:, vertex_of_cell(3, c)))
w(2) = triple(vxyz(:, vertex_of_cell(1, c)), p, vxyz(:, vertex_of_cell(3, c)))
w(3) = triple(vxyz(:, vertex_of_cell(1, c)), vxyz(:, vertex_of_cell(2, c)), p)
total = sum(w)
if (abs(total) > tiny(1.0_r8)) then
   w = w / total
   if (all(w >= tolerance)) then
      nc       = 1
      cells(1) = c
      wts(1)   = 1.0_r8
      return
   endif
endif

ier = OUTSIDE_DOMAIN

end subroutine find_triangle

!------------------------------------------------------------------
!> Vertical part of the interpolation. For every cell of the triangle and every
!> ensemble member: the two full levels around the location (kup is the higher one,
!> so the smaller index) and how far the location is from kup towards kdn.

subroutine find_levels(state_handle, ens_size, location, nc, cells, wts, &
                       kup, kdn, frac, on_surface, ier)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: nc
integer,             intent(in)  :: cells(3)
real(r8),            intent(in)  :: wts(3)
integer,             intent(out) :: kup(3, ens_size), kdn(3, ens_size)
real(r8),            intent(out) :: frac(3, ens_size)
logical,             intent(out) :: on_surface
integer,             intent(out) :: ier(ens_size)

integer  :: which_vert, ic, k, e, c, level, pres_id, psfc_id
real(r8) :: llv(3), vert, surface, ztop, zbot
real(r8) :: pcol(nlev, ens_size), psfc(ens_size), logp
logical  :: have_psfc

kup = nlev
kdn = nlev
frac = 0.0_r8
ier = 0
on_surface = .false.

llv        = get_location(location)
vert       = llv(3)
which_vert = nint(query_location(location))

surface = 0.0_r8
do ic = 1, nc
   surface = surface + wts(ic) * z_ifc(nlevp1, cells(ic))
enddo

select case (which_vert)

case (VERTISUNDEF)
   on_surface = .true.

case (VERTISSURFACE)
   on_surface = .true.
   if (sfc_elev_max_diff >= 0.0_r8) then
      if (abs(vert - surface) > sfc_elev_max_diff) ier = SURFACE_TOO_DIFFERENT
   endif

case (VERTISLEVEL)
   ! Counted from the bottom. The lowest level is both 0 and 1, and the top level
   ! itself already counts as above the model, as in the WRF interface.
   if (vert >= real(nlev, r8)) then
      ier = ABOVE_TOP_LEVEL
      return
   endif
   level = max(1, int(vert))
   kup = nlev - level + 1
   kdn = kup

case (VERTISHEIGHT)
   if (vert < surface) then
      ier = BELOW_SURFACE
      return
   endif
   do ic = 1, nc
      c = cells(ic)
      if (vert > z_mc(1, c)) then
         ier = ABOVE_MODEL_TOP
         return
      endif
      if (vert <= z_mc(nlev, c)) cycle     ! in the lowest layer: its value
      do k = 1, nlev - 1
         ztop = z_mc(k, c)
         zbot = z_mc(k + 1, c)
         if (vert >= zbot) then
            kup(ic, :)  = k
            kdn(ic, :)  = k + 1
            frac(ic, :) = (ztop - vert) / (ztop - zbot)
            exit
         endif
      enddo
   enddo

case (VERTISPRESSURE)
   pres_id = get_varid_from_kind(dom_id, QTY_PRESSURE)
   if (pres_id <= 0) then
      ier = NO_PRESSURE_IN_STATE
      return
   endif
   if (get_num_levels(pres_id) /= nlev) then
      ier = NO_PRESSURE_IN_STATE
      return
   endif
   psfc_id   = get_varid_from_kind(dom_id, QTY_SURFACE_PRESSURE)
   have_psfc = psfc_id > 0
   do ic = 1, nc
      c = cells(ic)
      do k = 1, nlev
         pcol(k, :) = get_state(get_dart_vector_index(c, k, 1, dom_id, pres_id), state_handle)
      enddo
      if (have_psfc) then
         psfc = get_state(get_dart_vector_index(c, 1, 1, dom_id, psfc_id), state_handle)
      endif
      do e = 1, ens_size
         if (ier(e) /= 0) cycle
         if (vert < pcol(1, e)) then
            ier(e) = ABOVE_MODEL_TOP
            cycle
         endif
         if (vert >= pcol(nlev, e)) then
            ! In the lowest layer, if we can tell that it is above the ground
            if (have_psfc) then
               if (vert > psfc(e)) ier(e) = BELOW_SURFACE
            else if (vert > pcol(nlev, e) * 1.002_r8) then
               ier(e) = BELOW_SURFACE
            endif
            cycle
         endif
         logp = log(vert)
         do k = 1, nlev - 1
            if (vert <= pcol(k + 1, e)) then
               kup(ic, e)  = k
               kdn(ic, e)  = k + 1
               frac(ic, e) = (logp - log(pcol(k, e))) / (log(pcol(k + 1, e)) - log(pcol(k, e)))
               exit
            endif
         enddo
      enddo
   enddo

case default
   ier = UNKNOWN_VERTICAL

end select

end subroutine find_levels

!------------------------------------------------------------------
!> Number of levels of a state variable: nlev, nlevp1 or 1.

function get_num_levels(var_id)

integer, intent(in) :: var_id
integer :: get_num_levels

if (get_num_dims(dom_id, var_id) >= 2) then
   get_num_levels = get_dim_length(dom_id, var_id, 2)
else
   get_num_levels = 1
endif

end function get_num_levels

!------------------------------------------------------------------
!> Interpolate one state variable to the location described by the triangle and
!> the levels.

subroutine interpolate_variable(state_handle, ens_size, var_id, nc, cells, wts, &
                                kup, kdn, frac, values, ier)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: var_id
integer,             intent(in)  :: nc
integer,             intent(in)  :: cells(3)
real(r8),            intent(in)  :: wts(3)
integer,             intent(in)  :: kup(3, ens_size), kdn(3, ens_size)
real(r8),            intent(in)  :: frac(3, ens_size)
real(r8),            intent(out) :: values(ens_size)
integer,             intent(out) :: ier

integer  :: ic, e, c, nlevels, kup_now, kdn_now
real(r8) :: upper(ens_size), lower(ens_size)

values = 0.0_r8
ier = 0
nlevels = get_num_levels(var_id)

if (nlevels == 1) then
   do ic = 1, nc
      upper  = get_state(get_dart_vector_index(cells(ic), 1, 1, dom_id, var_id), state_handle)
      values = values + wts(ic) * upper
   enddo
   return
endif

if (nlevels /= nlev .and. nlevels /= nlevp1) then
   ier = BAD_VERTICAL_FOR_FIELD
   return
endif

do ic = 1, nc
   c = cells(ic)
   kup_now = -1
   kdn_now = -1
   do e = 1, ens_size
      ! The levels are the same for all members unless the location is in pressure,
      ! so this rarely fetches more than once
      if (kup(ic, e) /= kup_now) then
         kup_now = kup(ic, e)
         upper = level_values(c, kup_now)
      endif
      if (kdn(ic, e) /= kdn_now) then
         kdn_now = kdn(ic, e)
         if (kdn_now == kup_now) then
            lower = upper
         else
            lower = level_values(c, kdn_now)
         endif
      endif
      values(e) = values(e) + wts(ic) * &
                  (upper(e) * (1.0_r8 - frac(ic, e)) + lower(e) * frac(ic, e))
   enddo
enddo

contains

   !> All members' values at full level k. Variables on half levels (w) are averaged
   !> to the full level.
   function level_values(cell, k)
      integer, intent(in) :: cell, k
      real(r8) :: level_values(ens_size)

      level_values = get_state(get_dart_vector_index(cell, k, 1, dom_id, var_id), state_handle)
      if (nlevels == nlevp1) then
         level_values = 0.5_r8 * (level_values + &
            get_state(get_dart_vector_index(cell, k + 1, 1, dom_id, var_id), state_handle))
      endif
   end function level_values

end subroutine interpolate_variable

!------------------------------------------------------------------
!> Value of a quantity at the location: straight from a state variable if there is
!> one, otherwise computed from other quantities.

recursive subroutine interpolate_quantity(state_handle, ens_size, qty, nc, cells, wts, &
                                          kup, kdn, frac, on_surface, values, ier)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: qty
integer,             intent(in)  :: nc
integer,             intent(in)  :: cells(3)
real(r8),            intent(in)  :: wts(3)
integer,             intent(in)  :: kup(3, ens_size), kdn(3, ens_size)
real(r8),            intent(in)  :: frac(3, ens_size)
logical,             intent(in)  :: on_surface
real(r8),            intent(out) :: values(ens_size)
integer,             intent(out) :: ier

integer  :: var_id, ic, e, surface_qty, upper_air_qty
real(r8) :: pres(ens_size), temp(ens_size), hum(ens_size)

values = MISSING_R8
ier = 0

! Heights are not part of the state, they never change
if (qty == QTY_SURFACE_ELEVATION) then
   values = 0.0_r8
   do ic = 1, nc
      values = values + wts(ic) * z_ifc(nlevp1, cells(ic))
   enddo
   return
endif

if (qty == QTY_GEOPOTENTIAL_HEIGHT .or. qty == QTY_GEOMETRIC_HEIGHT) then
   values = 0.0_r8
   do ic = 1, nc
      do e = 1, ens_size
         if (on_surface) then
            values(e) = values(e) + wts(ic) * z_ifc(nlevp1, cells(ic))
         else
            values(e) = values(e) + wts(ic) * &
                        (z_mc(kup(ic, e), cells(ic)) * (1.0_r8 - frac(ic, e)) + &
                         z_mc(kdn(ic, e), cells(ic)) * frac(ic, e))
         endif
      enddo
   enddo
   return
endif

! At the surface, prefer the screen level diagnostics if they are in the state
if (on_surface) then
   surface_qty = surface_equivalent(qty)
   if (surface_qty /= qty) then
      var_id = get_varid_from_kind(dom_id, surface_qty)
      if (var_id > 0) then
         call interpolate_variable(state_handle, ens_size, var_id, nc, cells, wts, &
                                   kup, kdn, frac, values, ier)
         return
      endif
   endif
endif

var_id = get_varid_from_kind(dom_id, qty)
if (var_id > 0) then
   if (on_surface .and. get_num_levels(var_id) == nlevp1) then
      ier = BAD_VERTICAL_FOR_FIELD
      return
   endif
   ! With on_surface, kup and kdn are the lowest level
   call interpolate_variable(state_handle, ens_size, var_id, nc, cells, wts, &
                             kup, kdn, frac, values, ier)
   return
endif

! Not in the state: the screen level quantities fall back to the lowest model level
upper_air_qty = upper_air_equivalent(qty)
if (upper_air_qty /= qty) then
   call interpolate_quantity(state_handle, ens_size, upper_air_qty, nc, cells, wts, &
                             kup, kdn, frac, on_surface, values, ier)
   return
endif

select case (qty)

case (QTY_SPECIFIC_HUMIDITY)
   call interpolate_quantity(state_handle, ens_size, QTY_VAPOR_MIXING_RATIO, nc, cells, &
                             wts, kup, kdn, frac, on_surface, hum, ier)
   if (ier /= 0) return
   values = hum / (1.0_r8 + hum)

case (QTY_VAPOR_MIXING_RATIO)
   var_id = get_varid_from_kind(dom_id, QTY_SPECIFIC_HUMIDITY)
   if (var_id <= 0) then
      ier = UNKNOWN_QUANTITY
      return
   endif
   call interpolate_quantity(state_handle, ens_size, QTY_SPECIFIC_HUMIDITY, nc, cells, &
                             wts, kup, kdn, frac, on_surface, hum, ier)
   if (ier /= 0) return
   values = hum / (1.0_r8 - hum)

case (QTY_DENSITY, QTY_POTENTIAL_TEMPERATURE)
   call interpolate_quantity(state_handle, ens_size, QTY_PRESSURE, nc, cells, wts, &
                             kup, kdn, frac, on_surface, pres, ier)
   if (ier /= 0) return
   call interpolate_quantity(state_handle, ens_size, QTY_TEMPERATURE, nc, cells, wts, &
                             kup, kdn, frac, on_surface, temp, ier)
   if (ier /= 0) return
   if (qty == QTY_POTENTIAL_TEMPERATURE) then
      values = temp * (p0ref / pres) ** (rd / cpd)
   else
      call interpolate_quantity(state_handle, ens_size, QTY_SPECIFIC_HUMIDITY, nc, cells, &
                                wts, kup, kdn, frac, on_surface, hum, ier)
      if (ier /= 0) return
      values = pres / (rd * temp * (1.0_r8 + vtmpc1 * hum))
   endif

case default
   ier = UNKNOWN_QUANTITY

end select

end subroutine interpolate_quantity

!------------------------------------------------------------------

function surface_equivalent(qty)

integer, intent(in) :: qty
integer :: surface_equivalent

select case (qty)
case (QTY_TEMPERATURE);       surface_equivalent = QTY_2M_TEMPERATURE
case (QTY_SPECIFIC_HUMIDITY); surface_equivalent = QTY_2M_SPECIFIC_HUMIDITY
case (QTY_U_WIND_COMPONENT);  surface_equivalent = QTY_10M_U_WIND_COMPONENT
case (QTY_V_WIND_COMPONENT);  surface_equivalent = QTY_10M_V_WIND_COMPONENT
case (QTY_PRESSURE);          surface_equivalent = QTY_SURFACE_PRESSURE
case default;                 surface_equivalent = qty
end select

end function surface_equivalent

!------------------------------------------------------------------

function upper_air_equivalent(qty)

integer, intent(in) :: qty
integer :: upper_air_equivalent

select case (qty)
case (QTY_2M_TEMPERATURE);       upper_air_equivalent = QTY_TEMPERATURE
case (QTY_2M_SPECIFIC_HUMIDITY); upper_air_equivalent = QTY_SPECIFIC_HUMIDITY
case (QTY_10M_U_WIND_COMPONENT); upper_air_equivalent = QTY_U_WIND_COMPONENT
case (QTY_10M_V_WIND_COMPONENT); upper_air_equivalent = QTY_V_WIND_COMPONENT
case default;                    upper_air_equivalent = qty
end select

end function upper_air_equivalent

!------------------------------------------------------------------
!> Interpolate the state to a location. istatus is 0 where it worked.

subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: qty
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer  :: nc, cells(3), ier
integer  :: kup(3, ens_size), kdn(3, ens_size)
real(r8) :: wts(3), frac(3, ens_size), llv(3), values(ens_size)
logical  :: on_surface

if (.not. module_initialized) call static_init_model

expected_obs = MISSING_R8
istatus      = 0

llv = get_location(location)

call find_triangle(llv(1), llv(2), nc, cells, wts, ier)
if (ier /= 0) then
   istatus = ier
   return
endif

call find_levels(state_handle, ens_size, location, nc, cells, wts, &
                 kup, kdn, frac, on_surface, istatus)
if (all(istatus /= 0)) return

call interpolate_quantity(state_handle, ens_size, qty, nc, cells, wts, &
                          kup, kdn, frac, on_surface, values, ier)
if (ier /= 0) then
   istatus = ier
   return
endif

where (istatus == 0) expected_obs = values
! Anything that is not a number (missing values in the model files)
where (istatus == 0 .and. .not. (abs(expected_obs) < huge(1.0_r8)))
   istatus = BAD_VALUE
   expected_obs = MISSING_R8
end where

if (debug > 5 .and. do_output()) then
   write(string1, *) 'qty ', qty, ' at ', llv, ' cells ', cells(1:nc)
   write(string2, *) 'weights ', wts(1:nc), ' value ', expected_obs(1), ' status ', istatus(1)
   call error_handler(E_MSG, 'model_interpolate', string1, source, text2=string2)
endif

end subroutine model_interpolate

!------------------------------------------------------------------
!> Location and quantity of a state element. Everything gets a height, so the
!> state never needs a vertical conversion.

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: qty

integer  :: cell, level, unused, var_id, my_qty, nlevels
real(r8) :: height

if (.not. module_initialized) call static_init_model

call get_model_variable_indices(index_in, cell, level, unused, &
                                var_id=var_id, kind_index=my_qty)

nlevels = get_num_levels(var_id)
if (nlevels == nlev) then
   height = z_mc(level, cell)
else if (nlevels == nlevp1) then
   height = z_ifc(level, cell)
else
   height = z_ifc(nlevp1, cell)
endif

location = set_location(clon(cell), clat(cell), height, VERTISHEIGHT)

if (present(qty)) qty = my_qty

end subroutine get_state_meta_data

!------------------------------------------------------------------
!> Give the locations of observations a height, which is the vertical coordinate
!> of the localisation.

subroutine convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, status)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer,             intent(in)    :: loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status(:)

integer  :: i, current, nc, cells(3), ier
integer  :: kup(3, 1), kdn(3, 1), istatus(1)
real(r8) :: wts(3), frac(3, 1), llv(3), height(1)
logical  :: on_surface

if (.not. module_initialized) call static_init_model

if (which_vert /= VERTISHEIGHT) then
   call error_handler(E_ERR, 'convert_vertical_obs', &
                      'the ICON interface only localises in height', source)
endif

do i = 1, num
   status(i) = 0
   current = nint(query_location(locs(i)))
   if (current == VERTISHEIGHT .or. current == VERTISUNDEF) cycle

   llv = get_location(locs(i))

   if (current == VERTISSURFACE) then
      ! The vertical value of a surface location is its elevation
      locs(i) = set_location(llv(1), llv(2), llv(3), VERTISHEIGHT)
      cycle
   endif

   call find_triangle(llv(1), llv(2), nc, cells, wts, ier)
   if (ier /= 0) then
      status(i) = ier
      cycle
   endif
   call find_levels(state_handle, 1, locs(i), nc, cells, wts, kup, kdn, frac, &
                    on_surface, istatus)
   if (istatus(1) /= 0) then
      status(i) = istatus(1)
      cycle
   endif
   call interpolate_quantity(state_handle, 1, QTY_GEOMETRIC_HEIGHT, nc, cells, wts, &
                             kup, kdn, frac, on_surface, height, ier)
   if (ier /= 0) then
      status(i) = ier
      cycle
   endif
   locs(i) = set_location(llv(1), llv(2), height(1), VERTISHEIGHT)
enddo

end subroutine convert_vertical_obs

!------------------------------------------------------------------

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc
type(location_type),           intent(inout) :: base_loc
integer,                       intent(in)    :: base_type
type(location_type),           intent(inout) :: locs(:)
integer,                       intent(in)    :: loc_qtys(:)
integer,                       intent(in)    :: loc_types(:)
integer,                       intent(out)   :: num_close
integer,                       intent(out)   :: close_ind(:)
real(r8),            optional, intent(out)   :: dist(:)
type(ensemble_type), optional, intent(in)    :: ens_handle

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                       num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs

!------------------------------------------------------------------

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc
type(location_type),           intent(inout) :: base_loc
integer,                       intent(in)    :: base_type
type(location_type),           intent(inout) :: locs(:)
integer,                       intent(in)    :: loc_qtys(:)
integer(i8),                   intent(in)    :: loc_indx(:)
integer,                       intent(out)   :: num_close
integer,                       intent(out)   :: close_ind(:)
real(r8),            optional, intent(out)   :: dist(:)
type(ensemble_type), optional, intent(in)    :: ens_handle

call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                         num_close, close_ind, dist, ens_handle)

end subroutine get_close_state

!------------------------------------------------------------------
!> Time of the (last record of the) model file. ICON writes either a relative axis
!> ("minutes since 2026-9-13 00:00:00") or an absolute one ("day as %Y%m%d.%f").

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type) :: read_model_time

character(len=*), parameter :: routine = 'read_model_time'

integer  :: ncid, varid, dimids(1), ntimes, io, pos, i
integer  :: year, month, day, hour, minute, second, days, seconds
real(r8) :: value, factor, rest
character(len=256) :: units, datestring
integer(i8) :: total

if (.not. module_initialized) call static_init_model

ncid = nc_open_file_readonly(filename, routine)

call nc_check(nf90_inq_varid(ncid, 'time', varid), routine, 'no time in '//trim(filename))
call nc_check(nf90_inquire_variable(ncid, varid, dimids=dimids), routine, 'inquire time')
call nc_check(nf90_inquire_dimension(ncid, dimids(1), len=ntimes), routine, 'time dimension')
call nc_check(nf90_get_var(ncid, varid, value, start=(/ ntimes /)), routine, 'get time')
units = ' '
call nc_check(nf90_get_att(ncid, varid, 'units', units), routine, 'get time units')
call nc_close_file(ncid, routine)

pos = index(units, ' since ')

if (pos > 0) then
   select case (trim(adjustl(units(1:pos))))
   case ('seconds', 'second'); factor = 1.0_r8
   case ('minutes', 'minute'); factor = 60.0_r8
   case ('hours', 'hour');     factor = 3600.0_r8
   case ('days', 'day');       factor = 86400.0_r8
   case default
      call error_handler(E_ERR, routine, 'unknown time units "'//trim(units)//'"', &
                         source, text2=trim(filename))
   end select

   datestring = units(pos + 7:)
   do i = 1, len_trim(datestring)
      if (index('-:TZ', datestring(i:i)) > 0) datestring(i:i) = ' '
   enddo
   hour = 0;  minute = 0;  rest = 0.0_r8
   read(datestring, *, iostat=io) year, month, day, hour, minute, rest
   if (io /= 0) read(datestring, *, iostat=io) year, month, day, hour, minute
   if (io /= 0) read(datestring, *, iostat=io) year, month, day
   if (io /= 0) then
      call error_handler(E_ERR, routine, 'cannot read the date in "'//trim(units)//'"', &
                         source, text2=trim(filename))
   endif

   total   = nint(value * factor + rest, i8)
   days    = int(total / 86400_i8)
   seconds = int(total - int(days, i8) * 86400_i8)
   read_model_time = set_date(year, month, day, hour, minute, 0) + set_time(seconds, days)

else if (index(units, '%Y%m%d') > 0) then
   year    = int(value / 10000.0_r8)
   month   = int((value - year * 10000.0_r8) / 100.0_r8)
   day     = int(value - year * 10000.0_r8 - month * 100.0_r8)
   seconds = nint((value - aint(value)) * 86400.0_r8)
   read_model_time = set_date(year, month, day, 0, 0, 0) + set_time(seconds, 0)

else
   call error_handler(E_ERR, routine, 'unknown time units "'//trim(units)//'"', &
                      source, text2=trim(filename))
endif

end function read_model_time

!------------------------------------------------------------------

subroutine end_model()

if (allocated(cell_locs)) then
   call xyz_get_close_destroy(cell_search)
   deallocate(cell_locs)
endif

end subroutine end_model

!------------------------------------------------------------------
!> Attributes for the files that DART creates itself (mean, sd, inflation).

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid
integer, intent(in) :: domain_id

if (.not. module_initialized) call static_init_model

call nc_begin_define_mode(ncid)
call nc_add_global_creation_time(ncid)
call nc_add_global_attribute(ncid, "model_source", source)
call nc_add_global_attribute(ncid, "model", "ICON")
call nc_add_global_attribute(ncid, "icon_grid_file", trim(grid_file))
call nc_end_define_mode(ncid)
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!------------------------------------------------------------------
! Small helpers
!------------------------------------------------------------------

pure function unit_vector(lon, lat)
real(r8), intent(in) :: lon, lat     ! radians
real(r8) :: unit_vector(3)
unit_vector = (/ cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat) /)
end function unit_vector

pure function cross(a, b)
real(r8), intent(in) :: a(3), b(3)
real(r8) :: cross(3)
cross = (/ a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1) /)
end function cross

!> a . (b x c), positive when a, b, c are seen anticlockwise from outside the sphere
pure function triple(a, b, c)
real(r8), intent(in) :: a(3), b(3), c(3)
real(r8) :: triple
triple = sum(a * cross(b, c))
end function triple

function dimension_size(ncid, name, filename)
integer,          intent(in) :: ncid
character(len=*), intent(in) :: name, filename
integer :: dimension_size
integer :: dimid
call nc_check(nf90_inq_dimid(ncid, name, dimid), 'dimension_size', &
              'no dimension '//trim(name)//' in '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=dimension_size), 'dimension_size', name)
end function dimension_size

subroutine get_real_1d(ncid, name, values, filename)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: name, filename
real(r8),         intent(out) :: values(:)
integer :: varid
call nc_check(nf90_inq_varid(ncid, name, varid), 'get_real_1d', &
              'no variable '//trim(name)//' in '//trim(filename))
call nc_check(nf90_get_var(ncid, varid, values), 'get_real_1d', name)
end subroutine get_real_1d

subroutine get_int_2d(ncid, name, values, filename)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: name, filename
integer,          intent(out) :: values(:,:)
integer :: varid
call nc_check(nf90_inq_varid(ncid, name, varid), 'get_int_2d', &
              'no variable '//trim(name)//' in '//trim(filename))
call nc_check(nf90_get_var(ncid, varid, values), 'get_int_2d', name)
end subroutine get_int_2d

end module model_mod
