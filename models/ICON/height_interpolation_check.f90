! DART software - Copyright UCAR. See the DART license for terms of use.
program height_interpolation_check
! Public-interface integration test. Only this driver accesses local ensemble
! storage, to populate known synthetic fields before creating DART's window.
use types_mod, only : r8,i8,PI,rad2deg,MISSING_R8
use mpi_utilities_mod, only : initialize_mpi_utilities,finalize_mpi_utilities,my_task_id,task_count
use model_mod, only : static_init_model,get_model_size,get_state_meta_data,model_interpolate,end_model
use location_mod, only : location_type,set_location,get_location,VERTISHEIGHT,VERTISPRESSURE
use obs_kind_mod, only : QTY_TEMPERATURE,QTY_U_WIND_COMPONENT,QTY_V_WIND_COMPONENT,QTY_VERTICAL_VELOCITY,QTY_PRESSURE
use ensemble_manager_mod, only : ensemble_type,init_ensemble_manager,end_ensemble_manager,set_num_extra_copies, &
   get_var_owner_index
use distributed_state_mod, only : create_state_window,free_state_window
use state_structure_mod, only : get_num_variables,get_kind_index,get_varid_from_kind,get_dart_vector_index
use icon_geometry_mod, only : geometry
use icon_domain_mod, only : native_contains,nearest_icon_cell
use icon_vertical_mod, only : height_bracket
use icon_interpolate_mod
use icon_search_init_mod, only : horizontal_search
use icon_horizontal_search_mod, only : query_horizontal_search,SEARCH_OK
use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan,ieee_next_after
implicit none
type(ensemble_type) :: ensemble
type(location_type) :: location
type(observation_stencil_type) :: stencil
integer :: i,j,k,m,q,phase,rc,counts(0:8),unit,lower,upper,ref_lower,ref_upper
integer :: qty,stat(3),offset
integer :: coverage(9),wind_kinds(2)
integer(i8) :: bad_index
real(r8) :: xyz(3),lon,lat,h,values(3),fraction,ref_fraction,z(5),xyzv(3,3)
real(r8) :: heights(8)
logical :: valid,ref_valid,found
character(len=64) :: filename
call initialize_mpi_utilities('height_interpolation_check')
call static_init_model()
counts=0
coverage=0
wind_kinds=[QTY_U_WIND_COMPONENT,QTY_V_WIND_COMPONENT]
write(filename,'(a,i0)') 'operator_results.',my_task_id()
open(newunit=unit,file=trim(filename),status='replace')

! Independent ascending MPAS default loop, with exact endpoints/interior and
! adjacent representable numbers. It is independent of the production scan.
z=[3000._r8,1800._r8,701._r8,500._r8,20._r8]
do i=1,5
   do j=-1,1
      h=z(i)
      if (j/=0) h=ieee_next_after(h,real(j,r8)*huge(h))
      call height_bracket(z,h,lower,upper,fraction,valid)
      call mpas_reference(z,h,ref_lower,ref_upper,ref_fraction,ref_valid)
      call require(valid.eqv.ref_valid,'bracket acceptance')
      if (.not. valid) cycle
      call require(lower==ref_lower.and.upper==ref_upper,'MPAS bracket indices')
      call require(fraction==ref_fraction,'MPAS bracket fraction')
   enddo
enddo
call height_bracket(z(1:1),z(1),lower,upper,fraction,valid)
call require(.not.valid,'one level rejected')

! Membership and nearest-cell index versus exhaustive independent distances
! and edge determinants. Vertices, native edge midpoints and cell centers.
do i=1,geometry%ncells,max(1,geometry%ncells/128)
   do j=0,6
      call native_vertices(i,xyzv)
      if (j==0) then
         xyz=geometry%xyz(:,i)
      else if (j<=3) then
         xyz=xyzv(:,j)
      else
         k=j-3
         xyz=xyzv(:,k)+xyzv(:,mod(k,3)+1)
         xyz=xyz/norm2(xyz)
      endif
      call check_domain(xyz)
   enddo
enddo
! Regular lon/lat probes across and beyond a regional bounding box also
! exercise holes and concavities. Global fixtures sample the whole sphere.
do i=0,20
   do j=0,20
      lon=minval(geometry%clon)+(maxval(geometry%clon)-minval(geometry%clon))*(real(i,r8)-2)/16
      lat=minval(geometry%clat)+(maxval(geometry%clat)-minval(geometry%clat))*(real(j,r8)-2)/16
      lat=max(-PI/2,min(PI/2,lat))
      xyz=[cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)]
      call check_domain(xyz)
   enddo
enddo

call init_ensemble_manager(ensemble,3,get_model_size())
call set_num_extra_copies(ensemble,0)
bad_index=0_i8
! Every phase creates one window collectively. Queries inside the window are
! deliberately different on each task (including early rejection paths).
do phase=1,4
   call populate(phase,0_i8)
   call create_state_window(ensemble)
   do q=1,get_num_variables(1)
      qty=get_kind_index(1,q)
      do i=my_task_id()+1,geometry%ncells,max(1,geometry%ncells/96)*task_count()
         xyz=geometry%xyz(:,i)
         call check_point(xyz,1000._r8,qty,phase)
         ! A zero-weight neighbor can have its lowest level above this query.
         call check_point(xyz,500._r8,qty,phase)
         call check_point(xyz,3000._r8,qty,phase)
      enddo
      do i=my_task_id()+1,size(horizontal_search%triangles%vertex,2), &
         max(1,size(horizontal_search%triangles%vertex,2)/96)*task_count()
         xyz=sum(horizontal_search%triangles%xyz(:,horizontal_search%triangles%vertex(:,i)),dim=2)
         xyz=xyz/norm2(xyz)
         heights=[1000._r8,3000._r8,1800._r8,0._r8,3500._r8, &
            ieee_next_after(3000._r8,huge(h)),ieee_next_after(3000._r8,-huge(h)),500._r8]
         do k=1,size(heights)
            call check_point(xyz,heights(k),qty,phase)
         enddo
      enddo
   enddo
   location=set_location(0._r8,0._r8,1000._r8,VERTISPRESSURE)
   call model_interpolate(ensemble,3,location,QTY_TEMPERATURE,values,stat)
   call require(all(stat==INTERP_UNSUPPORTED).and.all(values==MISSING_R8),'unsupported coordinate')
   location=set_location(0._r8,0._r8,1000._r8,VERTISHEIGHT)
   call model_interpolate(ensemble,3,location,QTY_PRESSURE,values,stat)
   call require(all(stat==INTERP_UNAVAILABLE).and.all(values==MISSING_R8),'unavailable pressure')
   location=set_location(0._r8,0._r8,ieee_value(0._r8,ieee_quiet_nan),VERTISHEIGHT)
   call model_interpolate(ensemble,3,location,QTY_TEMPERATURE,values,stat)
   call require(all(stat==INTERP_LOCATION).and.all(values==MISSING_R8),'nonfinite location')
   do j=1,2
      qty=wind_kinds(j)
      if (get_varid_from_kind(1,qty)>0) cycle
      call model_interpolate(ensemble,3,location,qty,values,stat)
      call require(all(stat==INTERP_UNAVAILABLE),'unconfigured quantity')
   enddo
   call free_state_window(ensemble)
enddo

! Find a valid full-level query. Corrupt each of its six dependencies in turn
! for members 2 (MISSING_R8) and 3 (NaN); at the common top, lower coefficients
! are exactly zero. Both public calls and direct stencil evaluation are used.
q=get_varid_from_kind(1,QTY_TEMPERATURE)
if (q>0) then
 found=.false.
 do i=1,geometry%ncells
   call lon_lat(geometry%xyz(:,i),lon,lat)
   location=set_location(lon,lat,3000._r8,VERTISHEIGHT)
   call build_observation_stencil(location,stencil,rc)
   if (rc/=0) cycle
   found=.true.
   exit
 enddo
 if (found) then
   do j=1,3
     do k=1,2
       offset=stencil%lower(j)
       if (k==2) offset=stencil%upper(j)
       bad_index=get_dart_vector_index(stencil%cells(j),offset,1,1,q)
       call populate(2,bad_index)
       call create_state_window(ensemble)
       call model_interpolate(ensemble,3,location,QTY_TEMPERATURE,values,stat)
       fraction=1._r8-stencil%fraction(j)
       if (k==2) fraction=stencil%fraction(j)
       call require(stat(1)==0,'uncorrupted member')
       call require(abs(values(1)-field(geometry%xyz(:,i),3000._r8,QTY_TEMPERATURE,1,2))<1.e-10_r8, &
          'uncorrupted member value')
       if (stencil%weights(j)*fraction==0._r8) then
          call require(all(stat==0),'zero coefficient ignores missing/NaN')
          coverage(5)=coverage(5)+1
          do m=1,3
            call require(abs(values(m)-field(geometry%xyz(:,i),3000._r8,QTY_TEMPERATURE,m,2))<1.e-10_r8, &
               'zero coefficient leaves value unchanged')
          enddo
       else
          call require(all(stat(2:3)==INTERP_INVALID_DATA),'member-local invalid data')
          call require(all(values(2:3)==MISSING_R8),'invalid members return missing')
          coverage(6)=coverage(6)+1
       endif
       call free_state_window(ensemble)
     enddo
   enddo
 endif
endif
call end_ensemble_manager(ensemble)
if (task_count()>1) call require(coverage(9)>0,'remote state samples exercised')
write(unit,'(a,9(1x,i0))') 'COVERAGE',coverage
write(unit,'(a,9(1x,i0))') 'STATUSES',counts
close(unit)
call end_model()
write(*,'(a,i0)') 'height_interpolation_check PASS task ',my_task_id()
call finalize_mpi_utilities()
contains
subroutine require(ok,message)
logical,intent(in) :: ok
character(len=*),intent(in) :: message
if (.not.ok) then
   write(*,*) 'FAIL task',my_task_id(),': ',message
   error stop 1
endif
end subroutine
subroutine mpas_reference(bounds,height,lo,hi,f,ok)
real(r8),intent(in) :: bounds(:),height
integer,intent(out) :: lo,hi
real(r8),intent(out) :: f
logical,intent(out) :: ok
real(r8) :: ascending(size(bounds))
integer :: a,n
n=size(bounds)
ascending=bounds(n:1:-1)
lo=0; hi=0; f=0
ok=.false.
if (n<2) return
if (height<ascending(1).or.height>ascending(n)) return
do a=2,n
 if (height<=ascending(a)) then
   lo=n-(a-1)+1; hi=n-a+1
   f=(height-ascending(a-1))/(ascending(a)-ascending(a-1))
   ok=.true.
   return
 endif
enddo
end subroutine
subroutine native_vertices(cell,v)
integer,intent(in) :: cell
real(r8),intent(out) :: v(3,3)
integer :: a,b
real(r8) :: x,y
do a=1,3
 b=geometry%vertices(cell,a)
 x=geometry%vlon(b); y=geometry%vlat(b)
 v(:,a)=[cos(y)*cos(x),cos(y)*sin(x),sin(y)]
enddo
end subroutine
function cross(a,b) result(c)
real(r8),intent(in) :: a(3),b(3)
real(r8) :: c(3)
c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
end function
logical function exhaustive(p) result(inside)
real(r8),intent(in) :: p(3)
real(r8) :: v(3,3),orient,det
integer :: a,b
logical :: ok
inside=.true.
if (geometry%global_grid) return
inside=.false.
do a=1,geometry%ncells
 call native_vertices(a,v)
 orient=dot_product(cross(v(:,1),v(:,2)),v(:,3))
 if (abs(orient)<=1.1e-14_r8) cycle
 ok=.true.
 do b=1,3
   det=dot_product(cross(v(:,b),v(:,mod(b,3)+1)),p)
   if (abs(det)>1.1e-14_r8.and.det*orient<0._r8) ok=.false.
 enddo
 if (ok) then
   inside=.true.
   return
 endif
enddo
end function
subroutine check_domain(p)
real(r8),intent(in) :: p(3)
real(r8) :: dmin,d
integer :: a,cell,code
call require(native_contains(p).eqv.exhaustive(p),'indexed versus exhaustive native membership')
call nearest_icon_cell(p,cell,code)
call require(code==0.and.cell>0,'nearest service')
dmin=huge(dmin)
do a=1,geometry%ncells
 d=norm2(p-geometry%xyz(:,a))
 dmin=min(dmin,d)
enddo
if (abs(norm2(p-geometry%xyz(:,cell))-dmin)>=2.e-14_r8) then
 write(*,*) 'nearest mismatch query/cell/dist/min: ',p,cell,norm2(p-geometry%xyz(:,cell)),dmin
endif
call require(abs(norm2(p-geometry%xyz(:,cell))-dmin)<2.e-14_r8,'nearest versus brute force chord distance')
end subroutine
subroutine lon_lat(p,x,y)
real(r8),intent(in) :: p(3)
real(r8),intent(out) :: x,y
x=modulo(atan2(p(2),p(1))*rad2deg,360._r8)
y=atan2(p(3),norm2(p(1:2)))*rad2deg
end subroutine
real(r8) function field(p,height,kind,member,mode) result(value)
real(r8),intent(in) :: p(3),height
integer,intent(in) :: kind,member,mode
value=real(10*member+kind,r8)
if (mode>=2) value=value+height*.01_r8
! A different horizontal profile tests weighting after per-column brackets.
if (mode>=3) value=value+2*p(1)-3*p(2)+4*p(3)
if (mode==4) value=value+height*height*1.e-5_r8
end function
subroutine populate(mode,bad)
integer,intent(in) :: mode
integer(i8),intent(in) :: bad
type(location_type) :: loc
real(r8) :: coordinates(3),p(3),x,y
integer :: a,b,kind
do a=1,ensemble%my_num_vars
 call get_state_meta_data(ensemble%my_vars(a),loc,kind)
 coordinates=get_location(loc)
 x=coordinates(1)/rad2deg; y=coordinates(2)/rad2deg
 p=[cos(y)*cos(x),cos(y)*sin(x),sin(y)]
 do b=1,3
   ensemble%copies(b,a)=field(p,coordinates(3),kind,b,mode)
 enddo
 if (ensemble%my_vars(a)==bad) then
   ensemble%copies(2,a)=MISSING_R8
   ensemble%copies(3,a)=ieee_value(0._r8,ieee_quiet_nan)
 endif
enddo
end subroutine
subroutine check_point(p,height,kind,mode)
real(r8),intent(in) :: p(3),height
integer,intent(in) :: kind,mode
type(location_type) :: loc
integer :: status(3),expected,a,b,lo,hi,cell,search_status,tri,corner(3),lower_levels(3),nvalid,owner,local_index,varid
real(r8) :: x,y,answer(3),reference(3),w(3),f,low,high,coord(3),minimum(3),maximum(3)
logical :: ok
call lon_lat(p,x,y)
loc=set_location(x,y,height,VERTISHEIGHT)
! Match the actual public-location round trip when obtaining reference weights.
coord=get_location(loc)
x=coord(1)/rad2deg; y=coord(2)/rad2deg
expected=INTERP_OK
reference=0._r8
minimum=huge(1._r8)
maximum=-huge(1._r8)
lower_levels=0
if (kind==QTY_VERTICAL_VELOCITY) then
 expected=INTERP_UNSUPPORTED
else if (.not.exhaustive(p)) then
 expected=INTERP_DOMAIN
else
 cell=minloc(sum((geometry%xyz-spread(p,2,geometry%ncells))**2,dim=1),dim=1)
 if (geometry%excluded(cell)) then
   expected=INTERP_BOUNDARY
   coverage(1)=coverage(1)+1
 else
   call query_horizontal_search(horizontal_search,x,y,tri,corner,w,search_status)
   if (search_status/=SEARCH_OK) then
     expected=INTERP_SEARCH
   else if (any(geometry%excluded(corner))) then
     expected=INTERP_BOUNDARY
     coverage(2)=coverage(2)+1
     if (any(geometry%excluded(corner).and.w==0._r8)) coverage(3)=coverage(3)+1
   else
     nvalid=0
     varid=get_varid_from_kind(1,kind)
     do a=1,3
       call mpas_reference(geometry%z_mc(corner(a),:),height,lo,hi,f,ok)
       if (.not.ok) then
         expected=INTERP_VERTICAL
         if (w(a)==0._r8) coverage(8)=coverage(8)+1
         cycle
       endif
       nvalid=nvalid+1
       lower_levels(a)=lo
       if (w(a)*(1-f)/=0._r8) then
          call get_var_owner_index(ensemble,get_dart_vector_index(corner(a),lo,1,1,varid),owner,local_index)
          if (owner/=my_task_id()) coverage(9)=coverage(9)+1
       endif
       if (w(a)*f/=0._r8) then
          call get_var_owner_index(ensemble,get_dart_vector_index(corner(a),hi,1,1,varid),owner,local_index)
          if (owner/=my_task_id()) coverage(9)=coverage(9)+1
       endif
       do b=1,3
         low=field(geometry%xyz(:,corner(a)),geometry%z_mc(corner(a),lo),kind,b,mode)
         high=field(geometry%xyz(:,corner(a)),geometry%z_mc(corner(a),hi),kind,b,mode)
         reference(b)=reference(b)+w(a)*((1-f)*low+f*high)
         minimum(b)=min(minimum(b),low,high)
         maximum(b)=max(maximum(b),low,high)
       enddo
     enddo
     if (nvalid>0.and.nvalid<3) coverage(7)=coverage(7)+1
   endif
 endif
endif
call model_interpolate(ensemble,3,loc,kind,answer,status)
if (any(status/=expected)) write(*,*) 'query ',p,height,kind,'expected',expected,'got',status
call require(all(status==expected),'public status')
if (expected==INTERP_OK) then
 if (any(lower_levels/=lower_levels(1))) coverage(4)=coverage(4)+1
 call require(all(answer>=minimum-1.e-10_r8).and.all(answer<=maximum+1.e-10_r8),'bounded result')
 call require(maxval(abs(answer-reference))<1.e-10_r8,'synthetic ensemble values')
 if (mode<3) then
   do b=1,3
     call require(abs(answer(b)-field(p,height,kind,b,mode))<1.e-10_r8,'analytic constant/height linear')
   enddo
 endif
else
 call require(all(answer==MISSING_R8),'rejected values missing')
endif
counts(expected)=counts(expected)+1
write(unit,'(4i8,3es25.16)') mode,kind,expected,my_task_id(),answer
end subroutine
end program
