! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module icon_geometry_mod

! Static, native-grid geometry only. No history/state I/O or interpolation.
! Arrays use NetCDF Fortran order: (cell, level), (cell, corner).
use types_mod, only : r8, PI
use utilities_mod, only : error_handler, E_ERR, E_MSG, find_namelist_in_file, &
                          check_namelist_read, do_nml_file, do_nml_term, nmlfileunit
use netcdf_utilities_mod, only : nc_check
use netcdf
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
implicit none
private

public :: init_icon_geometry, end_icon_geometry, geometry

! ICON 2026.04 mo_impl_constants_grf: grf_bdywidth_c=4,
! grf_nudge_start_c=5. Cell rows, NOT edge rows.
integer, parameter :: specified_zone_width = 4
character(len=*), parameter :: routine = 'icon_geometry_mod'

type :: geometry_type
   integer :: ncells = 0, nvertices = 0, nlevels = 0
   logical :: global_grid = .false.
   character(len=128) :: horizontal_uuid = '', vertical_uuid = ''
   real(r8) :: sphere_radius = 0.0_r8
   real(r8), allocatable :: clon(:), clat(:), vlon(:), vlat(:), xyz(:,:)
   real(r8), allocatable :: z_mc(:,:), z_ifc(:,:)
   integer, allocatable :: vertices(:,:), neighbors(:,:), boundary_row(:)
   logical, allocatable :: excluded(:)
end type

! Public read-only view; ownership and lifecycle stay in this module.
type(geometry_type), protected :: geometry
logical :: initialized = .false.

contains

subroutine init_icon_geometry()
integer :: unit, ios, grid, vertical, cell_dim, vertex_dim, corner_dim
integer :: mc_dims(2), ifc_dims(2), varid, global_flag, ninterfaces
character(len=128) :: uuid
character(len=512) :: message
character(len=512) :: horizontal_grid_file, vertical_geometry_file
integer :: nudge_zone_width
namelist /icon_geometry_nml/ horizontal_grid_file, vertical_geometry_file, nudge_zone_width

if (initialized) return
horizontal_grid_file = 'icon_grid.nc'
vertical_geometry_file = 'icon_geometry.nc'
nudge_zone_width = -1 ! Required explicitly for LAM; never infer from metadata.
call find_namelist_in_file('input.nml', 'icon_geometry_nml', unit)
read(unit, nml=icon_geometry_nml, iostat=ios)
call check_namelist_read(unit, ios, 'icon_geometry_nml')
if (do_nml_file()) write(nmlfileunit, nml=icon_geometry_nml)
if (do_nml_term()) write(*, nml=icon_geometry_nml)

call nc_check(nf90_open(trim(horizontal_grid_file), nf90_nowrite, grid), routine, 'open horizontal grid')
call nc_check(nf90_get_att(grid, nf90_global, 'uuidOfHGrid', geometry%horizontal_uuid), routine)
call require(len_trim(geometry%horizontal_uuid)>0, 'Empty horizontal grid UUID')
call nc_check(nf90_get_att(grid, nf90_global, 'global_grid', global_flag), routine)
call require(global_flag==0 .or. global_flag==1, 'global_grid must be 0 or 1')
geometry%global_grid = global_flag==1
call nc_check(nf90_get_att(grid, nf90_global, 'semi_major_axis', geometry%sphere_radius), routine)
call require(ieee_is_finite(geometry%sphere_radius), 'Nonfinite sphere radius')
call require(geometry%sphere_radius>0.0_r8, 'Sphere radius must be positive')
call dimension(grid, 'cell', cell_dim, geometry%ncells)
call dimension(grid, 'vertex', vertex_dim, geometry%nvertices)
call dimension(grid, 'nv', corner_dim, ios)
call require(ios==3, 'Native ICON cells must have three vertices')
allocate(geometry%clon(geometry%ncells), geometry%clat(geometry%ncells))
allocate(geometry%vlon(geometry%nvertices), geometry%vlat(geometry%nvertices))
call read_coordinate(grid, 'clon', cell_dim, geometry%clon, .false.)
call read_coordinate(grid, 'clat', cell_dim, geometry%clat, .true.)
call read_coordinate(grid, 'vlon', vertex_dim, geometry%vlon, .false.)
call read_coordinate(grid, 'vlat', vertex_dim, geometry%vlat, .true.)
allocate(geometry%vertices(geometry%ncells,3), geometry%neighbors(geometry%ncells,3))
call variable(grid, 'vertex_of_cell', [cell_dim,corner_dim], varid, .true.)
call nc_check(nf90_get_var(grid, varid, geometry%vertices), routine, 'vertex_of_cell')
call require(all(geometry%vertices>=1 .and. geometry%vertices<=geometry%nvertices), 'Invalid vertex_of_cell index')
call variable(grid, 'neighbor_cell_index', [cell_dim,corner_dim], varid, .true.)
call nc_check(nf90_get_var(grid, varid, geometry%neighbors), routine, 'neighbor_cell_index')
! Preserve -1; callers must never dereference this native missing-neighbor sentinel.
call require(all(geometry%neighbors==-1 .or. &
   (geometry%neighbors>=1 .and. geometry%neighbors<=geometry%ncells)), 'Invalid neighbor_cell_index')
allocate(geometry%boundary_row(geometry%ncells), geometry%excluded(geometry%ncells))
geometry%boundary_row = 0
geometry%excluded = .false.
if (.not. geometry%global_grid) then
   call require(nudge_zone_width>=0, 'LAM requires explicit nonnegative nudge_zone_width in icon_geometry_nml')
   call variable(grid, 'refin_c_ctrl', [cell_dim], varid, .true.)
   call nc_check(nf90_get_var(grid, varid, geometry%boundary_row), routine, 'refin_c_ctrl')
   ! This reader supports single-domain LAM, not nested-grid row classifications.
   call require(all(geometry%boundary_row>=0), 'Negative boundary rows: nested-grid mapping is not supported')
   call require(maxval(geometry%boundary_row)>=specified_zone_width, 'Insufficient specified boundary row coverage')
   call require(nudge_zone_width<=maxval(geometry%boundary_row)-specified_zone_width, &
                'Insufficient refin_c_ctrl coverage for specified plus nudging width')
   geometry%excluded = geometry%boundary_row>=1 .and. &
                      geometry%boundary_row<=specified_zone_width+nudge_zone_width
endif
call nc_check(nf90_close(grid), routine)

call nc_check(nf90_open(trim(vertical_geometry_file), nf90_nowrite, vertical), routine, 'open vertical geometry')
call nc_check(nf90_get_att(vertical, nf90_global, 'uuidOfHGrid', uuid), routine)
call require(uuid==geometry%horizontal_uuid, 'Horizontal UUID mismatch between geometry inputs')
call nc_check(nf90_get_att(vertical, nf90_global, 'uuidOfVGrid', geometry%vertical_uuid), routine)
call require(len_trim(geometry%vertical_uuid)>0, 'Empty vertical grid UUID')
! The native stream can name level axes height/height_2 in either order.
! Associate them with z_mc/z_ifc, not with a presumed suffix.
call height_dimensions(vertical, 'z_mc', mc_dims, geometry%nlevels)
call height_dimensions(vertical, 'z_ifc', ifc_dims, ninterfaces)
call require(mc_dims(1)==ifc_dims(1), 'Height fields must share the cell dimension')
call require(ninterfaces==geometry%nlevels+1, 'Expected K full levels and K+1 interfaces')
allocate(geometry%z_mc(geometry%ncells,geometry%nlevels))
allocate(geometry%z_ifc(geometry%ncells,ninterfaces))
call read_height(vertical, 'z_mc', mc_dims, geometry%z_mc)
call read_height(vertical, 'z_ifc', ifc_dims, geometry%z_ifc)
call nc_check(nf90_close(vertical), routine)
call require(all(geometry%z_ifc(:,:geometry%nlevels)>geometry%z_mc) .and. &
             all(geometry%z_mc>geometry%z_ifc(:,2:)), 'Full-level heights must lie strictly between their interfaces')
call require(all(geometry%z_mc(:,:geometry%nlevels-1)>geometry%z_mc(:,2:)), &
             'Full-level heights must decrease from top to bottom')

allocate(geometry%xyz(3,geometry%ncells))
geometry%xyz(1,:) = cos(geometry%clat)*cos(geometry%clon)
geometry%xyz(2,:) = cos(geometry%clat)*sin(geometry%clon)
geometry%xyz(3,:) = sin(geometry%clat)
initialized = .true.
write(message,'(A,I0,A,I0,A,I0)') 'Geometry: cells=',geometry%ncells, &
   ', full levels=',geometry%nlevels, ', excluded cells=',count(geometry%excluded)
call error_handler(E_MSG, routine, trim(message))
end subroutine

subroutine dimension(ncid, name, dimid, length)
integer, intent(in) :: ncid
character(len=*), intent(in) :: name
integer, intent(out) :: dimid, length
call nc_check(nf90_inq_dimid(ncid, name, dimid), routine, name)
call nc_check(nf90_inquire_dimension(ncid, dimid, len=length), routine, name)
call require(length>0, 'Empty dimension '//name)
end subroutine

subroutine variable(ncid, name, expected_dims, varid, integer_field)
integer, intent(in) :: ncid, expected_dims(:)
character(len=*), intent(in) :: name
integer, intent(out) :: varid
logical, intent(in) :: integer_field
integer :: rank, xtype, dims(nf90_max_var_dims), status
call nc_check(nf90_inq_varid(ncid, name, varid), routine, name)
call nc_check(nf90_inquire_variable(ncid,varid,xtype=xtype,ndims=rank,dimids=dims), routine, name)
call require(rank==size(expected_dims), 'Wrong rank for '//name)
call require(all(dims(:rank)==expected_dims), 'Wrong dimension order for '//name)
if (integer_field) then
   call require(xtype==nf90_int, 'Expected native integer field '//name)
else
   call require(xtype==nf90_float .or. xtype==nf90_double, 'Expected floating-point field '//name)
endif
status = nf90_inquire_attribute(ncid,varid,'scale_factor')
call require(status==nf90_enotatt, 'Packed fields are unsupported: '//name)
status = nf90_inquire_attribute(ncid,varid,'add_offset')
call require(status==nf90_enotatt, 'Packed fields are unsupported: '//name)
end subroutine

subroutine read_coordinate(ncid, name, dimid, values, latitude)
integer, intent(in) :: ncid, dimid
character(len=*), intent(in) :: name
real(r8), intent(out) :: values(:)
logical, intent(in) :: latitude
integer :: varid
character(len=64) :: units
call variable(ncid, name, [dimid], varid, .false.)
call nc_check(nf90_get_att(ncid,varid,'units',units), routine, name)
call require(trim(adjustl(units))=='radian' .or. trim(adjustl(units))=='radians', 'Expected radians for '//name)
call nc_check(nf90_get_var(ncid,varid,values), routine, name)
call valid_values(ncid,varid,name,values)
if (latitude) then
   call require(all(abs(values)<=PI/2.0_r8), 'Latitude outside valid range: '//name)
else
   call require(all(abs(values)<=2.0_r8*PI), 'Longitude outside valid range: '//name)
endif
end subroutine

subroutine height_dimensions(ncid, name, dims, levels)
integer, intent(in) :: ncid
character(len=*), intent(in) :: name
integer, intent(out) :: dims(2), levels
integer :: varid, rank, all_dims(nf90_max_var_dims), cells
call nc_check(nf90_inq_varid(ncid,name,varid), routine, name)
call nc_check(nf90_inquire_variable(ncid,varid,ndims=rank,dimids=all_dims), routine, name)
call require(rank==2, 'Expected static two-dimensional '//name)
dims = all_dims(:2)
call nc_check(nf90_inquire_dimension(ncid,dims(1),len=cells), routine, name)
call nc_check(nf90_inquire_dimension(ncid,dims(2),len=levels), routine, name)
call require(cells==geometry%ncells, 'Cell count/order mismatch for '//name)
call require(levels>0, 'Empty vertical axis for '//name)
end subroutine

subroutine read_height(ncid, name, dims, values)
integer, intent(in) :: ncid, dims(2)
character(len=*), intent(in) :: name
real(r8), intent(out) :: values(:,:)
integer :: varid, k
character(len=64) :: units
call variable(ncid,name,dims,varid,.false.)
call nc_check(nf90_get_att(ncid,varid,'units',units), routine, name)
call require(trim(adjustl(units))=='m', 'Expected metres for '//name)
call nc_check(nf90_get_var(ncid,varid,values), routine, name)
do k=1,size(values,2)
   call valid_values(ncid,varid,name,values(:,k))
enddo
end subroutine

subroutine valid_values(ncid, varid, name, values)
integer, intent(in) :: ncid, varid
character(len=*), intent(in) :: name
real(r8), intent(in) :: values(:)
integer :: status, i, length
real(r8) :: missing
character(len=13), parameter :: attributes(2) = [character(len=13) :: '_FillValue', 'missing_value']
call require(all(ieee_is_finite(values)), 'Nonfinite values in '//name)
call require(all(values/=real(nf90_fill_double,r8)), 'Default fill values in '//name)
call require(all(values/=real(nf90_fill_float,r8)), 'Default fill values in '//name)
do i=1,2
   status = nf90_inquire_attribute(ncid,varid,trim(attributes(i)),len=length)
   if (status==nf90_enotatt) cycle
   call nc_check(status, routine, name)
   call require(length==1, 'Expected scalar missing-value attribute for '//name)
   call nc_check(nf90_get_att(ncid,varid,trim(attributes(i)),missing), routine, name)
   call require(all(values/=missing), 'Missing values in '//name)
enddo
end subroutine

subroutine require(condition, message)
logical, intent(in) :: condition
character(len=*), intent(in) :: message
if (.not. condition) call error_handler(E_ERR, routine, message)
end subroutine

subroutine end_icon_geometry()
type(geometry_type) :: empty
! Intrinsic assignment deallocates all allocatable components.
geometry = empty
initialized = .false.
end subroutine

end module icon_geometry_mod
