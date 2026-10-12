! DART software - Copyright UCAR. See the DART license for terms of use.
program state_metadata_check
use types_mod, only : r8, i8, RAD2DEG
use utilities_mod, only : initialize_utilities, finalize_utilities, error_handler, E_ERR
use model_mod, only : static_init_model, get_model_size, get_state_meta_data, end_model
use state_structure_mod, only : get_num_domains, get_num_variables, get_dim_length, &
   get_index_start, get_index_end, get_dart_vector_index, get_model_variable_indices, &
   get_variable_name, get_kind_index
use location_mod, only : location_type, get_location, query_location, VERTISHEIGHT
use icon_geometry_mod, only : geometry
implicit none
integer :: v, c, k, levs, qty, ci, ki, unused, vi, di, stride
integer(i8) :: idx, expected_size, expected_index
real(r8) :: coords(3), expected_height
type(location_type) :: loc
character(len=128) :: name, argument

call initialize_utilities('state_metadata_check')
call static_init_model()
call static_init_model()
call check(get_num_domains()==1,'Repeated initialization registered extra domains')
! Optional index lets the external test harness exercise invalid bounds.
call get_command_argument(1,argument)
if (len_trim(argument)>0) then
   read(argument,*) idx
   call get_state_meta_data(idx,loc)
   call check(.false.,'Invalid index was accepted')
endif
expected_size=0_i8
stride=1
if (geometry%ncells>1000) stride=geometry%ncells-1
do v=1,get_num_variables(1)
   name=get_variable_name(1,v)
   levs=geometry%nlevels
   if (trim(name)=='w') levs=levs+1
   call check(get_dim_length(1,v,1)==geometry%ncells,'Wrong cell dimension')
   call check(get_dim_length(1,v,2)==levs,'Wrong level dimension')
   call check(get_index_start(1,v)==expected_size+1_i8,'Wrong field offset')
   do k=1,levs
      do c=1,geometry%ncells,stride
         idx=get_dart_vector_index(c,k,1,1,v)
         expected_index=expected_size+int(k-1,i8)*geometry%ncells+c
         call check(idx==expected_index,'Flattened index differs from native cell-fast order')
         call get_model_variable_indices(idx,ci,ki,unused,vi,di)
         call check(ci==c .and. ki==k .and. vi==v .and. di==1,'Inverse indexing mismatch')
         call get_state_meta_data(idx,loc,qty)
         coords=get_location(loc)
         expected_height=geometry%z_mc(c,min(k,geometry%nlevels))
         if (trim(name)=='w') expected_height=geometry%z_ifc(c,k)
         call check(abs(coords(1)-modulo(geometry%clon(c)*RAD2DEG,360.0_r8))<1.e-10_r8,'Longitude mismatch')
         call check(abs(coords(2)-geometry%clat(c)*RAD2DEG)<1.e-10_r8,'Latitude mismatch')
         call check(abs(coords(3)-expected_height)<1.e-10_r8,'Height mismatch')
         call check(nint(query_location(loc,'which_vert'))==VERTISHEIGHT,'Incorrect vertical type')
         call check(qty==get_kind_index(1,v),'Quantity mismatch')
         call get_state_meta_data(-idx,loc) ! DART identity-index convention; optional qty omitted.
         call check(all(abs(get_location(loc)-coords)<1.e-10_r8),'Identity-index location mismatch')
      enddo
   enddo
   expected_size=expected_size+int(geometry%ncells,i8)*levs
   call check(get_index_end(1,v)==expected_size,'Wrong field end index')
   write(*,'(A,I0,3A)') 'FIELD ',v,': ',trim(name),' OK'
enddo
call check(get_model_size()==expected_size,'Wrong state size')
write(*,'(A,I0)') 'STATE SIZE ',expected_size
call end_model()
write(*,'(A)') 'state_metadata_check PASS'
call finalize_utilities()
contains
subroutine check(condition,message)
logical, intent(in) :: condition
character(len=*), intent(in) :: message
if (.not. condition) call error_handler(E_ERR,'state_metadata_check',message)
end subroutine
end program
