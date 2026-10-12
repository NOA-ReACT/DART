! DART software - Copyright UCAR. See the DART license for terms of use.
module icon_state_mod

! ICON-specific field semantics and locations; DART owns flattening and I/O.
use types_mod, only : r8, i8, RAD2DEG, vtablenamelength
use utilities_mod, only : error_handler, E_ERR
use default_model_mod, only : parse_variables, state_var_type, MAX_STATE_VARIABLE_FIELDS
use state_structure_mod, only : add_domain, get_model_variable_indices, get_domain_size
use obs_kind_mod, only : get_index_for_quantity
use location_mod, only : location_type, set_location, VERTISHEIGHT
use netcdf_utilities_mod, only : nc_check
use netcdf
use icon_geometry_mod, only : geometry
implicit none
private
public :: register_icon_state, icon_state_location, end_icon_state, default_icon_state

integer, parameter :: nfields=17
character(len=8), parameter :: names(nfields) = [character(len=8) :: &
   'u','v','w','temp','qv','qc','qi','qr','qs','qg','qh','qnc','qni','qnr','qns','qng','qnh']
character(len=40), parameter :: quantities(nfields) = [character(len=40) :: &
   'QTY_U_WIND_COMPONENT', 'QTY_V_WIND_COMPONENT', 'QTY_VERTICAL_VELOCITY', &
   'QTY_TEMPERATURE', 'QTY_SPECIFIC_HUMIDITY', 'QTY_CLOUDWATER_MIXING_RATIO', &
   'QTY_ICE_MIXING_RATIO', 'QTY_RAINWATER_MIXING_RATIO', 'QTY_SNOW_MIXING_RATIO', &
   'QTY_GRAUPEL_MIXING_RATIO', 'QTY_HAIL_MIXING_RATIO', 'QTY_DROPLET_NUMBER_CONCENTR', &
   'QTY_ICE_NUMBER_CONCENTRATION', 'QTY_RAIN_NUMBER_CONCENTR', 'QTY_SNOW_NUMBER_CONCENTR', &
   'QTY_GRAUPEL_NUMBER_CONCENTR', 'QTY_HAIL_NUMBER_CONCENTR']
character(len=8), parameter :: units(nfields) = [character(len=8) :: &
   'm s-1','m s-1','m s-1','K','kg kg-1','kg kg-1','kg kg-1','kg kg-1','kg kg-1', &
   'kg kg-1','kg kg-1','kg-1','kg-1','kg-1','kg-1','kg-1','kg-1']
integer :: domain_id=-1
logical, allocatable :: on_interfaces(:)
character(len=*), parameter :: routine='icon_state_mod'
contains

subroutine default_icon_state(table)
character(len=vtablenamelength), intent(out) :: table(MAX_STATE_VARIABLE_FIELDS)
integer :: j
table=''
do j=1,7
   table(3*j-2:3*j)=[character(len=vtablenamelength) :: names(j), quantities(j), 'UPDATE']
enddo
end subroutine

function register_icon_state(filename, table) result(id)
character(len=*), intent(in) :: filename
character(len=vtablenamelength), intent(in) :: table(MAX_STATE_VARIABLE_FIELDS)
integer :: id, ncid, j, entry, k, last, rank, xtype, varid, dims(nf90_max_var_dims)
integer :: sizes(3), timedim, cell_dim, full_dim, ifc_dim, expected_levels, status
character(len=128) :: uuid, actual_units, dimname
type(state_var_type) :: vars

call require(domain_id<0, 'State already registered')
! Check table shape before the standard parser (which stops at the first blank row).
last=0
do j=1,size(table)
   if (len_trim(table(j))>0) last=j
enddo
call require(last>0 .and. mod(last,3)==0, 'state_variables requires complete name/quantity/update triples')
call require(last<=3*nfields, 'Too many ICON state entries')
call require(all(len_trim(table(:last))>0), 'Empty entry inside state_variables')
vars=parse_variables(table)
allocate(on_interfaces(vars%nvars))
call nc_check(nf90_open(trim(filename),nf90_nowrite,ncid),routine,filename)
call nc_check(nf90_get_att(ncid,nf90_global,'uuidOfHGrid',uuid),routine)
call require(uuid==geometry%horizontal_uuid,'State horizontal UUID mismatch')
call nc_check(nf90_get_att(ncid,nf90_global,'uuidOfVGrid',uuid),routine)
call require(uuid==geometry%vertical_uuid,'State vertical UUID mismatch')
call nc_check(nf90_inquire(ncid,unlimitedDimId=timedim),routine)
call require(timedim>=0,'State requires a singleton unlimited time dimension')
cell_dim=-1
full_dim=-1
ifc_dim=-1
do j=1,vars%nvars
   entry=0
   do k=1,nfields
      if (vars%netcdf_var_names(j)==names(k)) entry=k
   enddo
   call require(entry>0,'Unsupported state field (pressure deferred): '//trim(vars%netcdf_var_names(j)))
   call require(vars%qtys(j)==get_index_for_quantity(quantities(entry)), &
                'Incorrect quantity for '//trim(names(entry)))
   call require(.not. any(vars%netcdf_var_names(:j-1)==vars%netcdf_var_names(j)), 'Duplicate state variable')
   on_interfaces(j)=names(entry)=='w'
   call nc_check(nf90_inq_varid(ncid,trim(names(entry)),varid),routine,trim(names(entry)))
   call nc_check(nf90_inquire_variable(ncid,varid,xtype=xtype,ndims=rank,dimids=dims),routine)
   call require(rank==3,'Expected CDL field(time,level,cell): '//trim(names(entry)))
   call require(xtype==nf90_float .or. xtype==nf90_double,'Expected floating-point state field')
   do k=1,3
      call nc_check(nf90_inquire_dimension(ncid,dims(k),len=sizes(k)),routine)
   enddo
   call require(dims(3)==timedim .and. sizes(3)==1,'Expected one unlimited time record')
   call nc_check(nf90_inquire_dimension(ncid,timedim,name=dimname),routine)
   call require(trim(dimname)=='time','Expected time dimension named time')
   call require(sizes(1)==geometry%ncells,'State cell count/order mismatch')
   if (cell_dim<0) cell_dim=dims(1)
   call require(cell_dim==dims(1),'State variables must share the cell dimension')
   expected_levels=geometry%nlevels
   if (on_interfaces(j)) then
      expected_levels=expected_levels+1
      if (ifc_dim<0) ifc_dim=dims(2)
      call require(ifc_dim==dims(2),'Inconsistent interface-level axis')
   else
      if (full_dim<0) full_dim=dims(2)
      call require(full_dim==dims(2),'Inconsistent full-level axis')
   endif
   call require(sizes(2)==expected_levels,'Incorrect vertical staggering/level count: '//trim(names(entry)))
   call nc_check(nf90_get_att(ncid,varid,'units',actual_units),routine,trim(names(entry)))
   call require(trim(adjustl(actual_units))==units(entry),'Incorrect units for '//trim(names(entry)))
   status=nf90_inquire_attribute(ncid,varid,'scale_factor')
   call require(status==nf90_enotatt,'Packed state fields are unsupported')
   status=nf90_inquire_attribute(ncid,varid,'add_offset')
   call require(status==nf90_enotatt,'Packed state fields are unsupported')
enddo
call nc_check(nf90_close(ncid),routine)
! Preserve DART's variable order, indexing, quantity table and copy-back flags.
domain_id=add_domain(filename,vars)
id=domain_id
end function

subroutine icon_state_location(index_in, location, qty)
integer(i8), intent(in) :: index_in
type(location_type), intent(out) :: location
integer, intent(out), optional :: qty
integer :: cell, level, unused, varid, domid, kind
real(r8) :: height, longitude
call require(domain_id>0,'State has not been registered')
! DART permits negative indices for identity observations.
call require(index_in/=0_i8 .and. index_in<=get_domain_size(domain_id) .and. &
             index_in>=-get_domain_size(domain_id),'State index out of range')
call get_model_variable_indices(index_in,cell,level,unused,varid,domid,kind)
call require(domid==domain_id,'State index belongs to another domain')
if (on_interfaces(varid)) then
   height=geometry%z_ifc(cell,level)
else
   height=geometry%z_mc(cell,level)
endif
longitude=modulo(geometry%clon(cell)*RAD2DEG,360.0_r8)
location=set_location(longitude,geometry%clat(cell)*RAD2DEG,height,VERTISHEIGHT)
if (present(qty)) qty=kind
end subroutine

subroutine end_icon_state()
if (allocated(on_interfaces)) deallocate(on_interfaces)
! DART owns the domain registry; end_model is terminal for this model instance.
end subroutine

subroutine require(condition,message)
logical, intent(in) :: condition
character(len=*), intent(in) :: message
if (.not. condition) call error_handler(E_ERR,routine,message)
end subroutine
end module icon_state_mod
