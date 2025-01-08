! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! LIDAR_EXTINCTION,  QTY_ATM_LIDAR_EXTINCTION
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_atm_lidar_mod, only : get_extinction
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(LIDAR_EXTINCTION)
!         call get_extinction(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(LIDAR_EXTINCTION)
!      continue
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(LIDAR_EXTINCTION)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(LIDAR_EXTINCTION)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_atm_lidar_mod

   use        types_mod, only : r8, PI, metadatalength, MISSING_R8, DEG2RAD
   use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
      logfileunit, get_unit, open_file, close_file, &
      file_exist, ascii_file_format
   use     location_mod, only : location_type, set_location, get_location, &
      VERTISHEIGHT, VERTISLEVEL, set_location_missing
   use     obs_kind_mod, only : QTY_ATM_LIDAR_EXTINCTION, QTY_DENSITY, &
      QTY_GC_DUST_BIN1, QTY_GC_DUST_BIN2, QTY_GC_DUST_BIN3, &
      QTY_GC_DUST_BIN4, QTY_GC_DUST_BIN5
   use  assim_model_mod, only : interpolate

   use obs_def_utilities_mod, only : track_status
   use ensemble_manager_mod,  only : ensemble_type

   implicit none
   private

   public :: get_extinction

! version controlled file description for error handling, do not edit
   character(len=256), parameter :: source   = "$URL$"
   character(len=32 ), parameter :: revision = "$Revision$"
   character(len=128), parameter :: revdate  = "$Date$"

   character(len=256) :: string1, string2
   logical, save      :: module_initialized = .false.

   ! For loops
   integer :: i

   type optical_properties
      private
      integer :: wavelength ! Wavelength of the laser [nm]
      integer :: size_bin ! Size bin number
      real(r8) :: qext ! Extinction efficiency [unitless]
      real(r8) :: rho ! Particle density [kg/m^3]
      real(r8) :: effective_diameter ! Effective diameter [m]
   end type optical_properties

   type(optical_properties), allocatable, dimension(:) :: optical_props

   logical :: debug = .FALSE.

contains

   ! Read the optical properties from a CSV file
   subroutine initialize_module
      integer :: f, fstat ! File handle for optical properties
      type(optical_properties) :: line

      if (module_initialized) return
      call register_module(source, revision, revdate)

      open(action='read', file='atm_lidar_optical_properties.csv', newunit=f, iostat=fstat)
      if (fstat /= 0) then
         write(string1, *) 'Could not open optical properties file atm_lidar_optical_properties.csv'
         call error_handler(E_ERR, 'initialize_module', string1, source, revision, revdate)
      endif

      ! Skip the first line, should be the header
      read(f, *)

      ! Read each line into an optical_properties object
      allocate(optical_props(0))
      do
         read(f, *, iostat=fstat) line
         if (fstat /= 0) exit

         ! Fix units
         line%rho = line%rho * 1e6
         line%effective_diameter = line%effective_diameter * 1e-6

         optical_props = [optical_props, line]
      end do

      ! Print out the optical properties for debugging
      if (debug) then
         do i = 1, size(optical_props)
            write(*, *) 'Wavelength: ', optical_props(i)%wavelength
            write(*, *) 'Size bin: ', optical_props(i)%size_bin
            write(*, *) 'Qext: ', optical_props(i)%qext
            write(*, *) 'Rho: ', optical_props(i)%rho
            write(*, *) 'Effective diameter: ', optical_props(i)%effective_diameter
         end do
      endif

      close(f)
      module_initialized = .true.
   end subroutine initialize_module

   ! Returns the optical properties for the given wavelength and size bin
   subroutine get_optical_props(wavelength, bin, out)
      integer, intent(in) :: wavelength
      integer, intent(in) :: bin
      type(optical_properties), intent(out) :: out

      do i = 1, size(optical_props)
         if (optical_props(i)%wavelength == wavelength .and. optical_props(i)%size_bin == bin) then
            out = optical_props(i)
            return
         endif
      end do

      write(string1, *) 'Could not find optical properties for wavelength ', wavelength, ' and size bin ', bin
      call error_handler(E_ERR, 'get_optical_props', string1, source, revision, revdate)
   end subroutine get_optical_props

   ! Computes the extinction from the dust mixing ratios
   subroutine get_extinction(state_handle, ens_size, location, key, extinction, istatus)
      type(ensemble_type), intent(in) :: state_handle
      integer, intent(in) :: ens_size
      type(location_type), intent(in) :: location
      integer, intent(in) :: key
      real(r8), intent(inout) :: extinction(ens_size)
      integer, intent(out) :: istatus(ens_size)

      type(optical_properties) :: props
      real(r8) :: rho(ens_size), dust_1(ens_size), dust_2(ens_size), dust_3(ens_size), dust_4(ens_size), dust_5(ens_size)
      integer :: this_istatus(ens_size)
      logical :: return_now = .false.
      real(r8) lon, lat, height, obsloc(3)

      ! obsloc   = get_location(location)
      ! lon      = obsloc(1)                       ! degree: 0 to 360
      ! lat      = obsloc(2)                       ! degree: -90 to 90
      ! height   = obsloc(3)                       ! (m)
      ! ! write(*,*) 'Location: (x, y, z)', lon, lat, height

      ! Initialize arrays
      rho = 0.0_r8
      dust_1 = 0.0_r8
      dust_2 = 0.0_r8
      dust_3 = 0.0_r8
      dust_4 = 0.0_r8
      dust_5 = 0.0_r8
      istatus = 0 ! All success
      this_istatus = 0

      call initialize_module

      ! Grab air density (rho, [kg/m^3]) and dust mixing ratios from the state
      call interpolate(state_handle, ens_size, location, QTY_DENSITY, rho, this_istatus)
      call track_status(ens_size, this_istatus, extinction, istatus, return_now)
      if (return_now) return

      call interpolate(state_handle, ens_size, location, QTY_GC_DUST_BIN1, dust_1, this_istatus)
      call track_status(ens_size, this_istatus, extinction, istatus, return_now)
      if (return_now) return

      call interpolate(state_handle, ens_size, location, QTY_GC_DUST_BIN2, dust_2, this_istatus)
      call track_status(ens_size, this_istatus, extinction, istatus, return_now)
      if (return_now) return

      call interpolate(state_handle, ens_size, location, QTY_GC_DUST_BIN3, dust_3, this_istatus)
      call track_status(ens_size, this_istatus, extinction, istatus, return_now)
      if (return_now) return

      call interpolate(state_handle, ens_size, location, QTY_GC_DUST_BIN4, dust_4, this_istatus)
      call track_status(ens_size, this_istatus, extinction, istatus, return_now)
      if (return_now) return

      call interpolate(state_handle, ens_size, location, QTY_GC_DUST_BIN5, dust_5, this_istatus)
      call track_status(ens_size, this_istatus, extinction, istatus, return_now)
      if (return_now) return

      ! Compute the extinction for each size bin
      ! Extinction = (3 * Qext * MixingRatio * rho) / (2 * EffectiveDiameter * ParticleDensity)

      call get_optical_props(355, 1, props)
      where (istatus == 0) &
         dust_1 = (3 * props%qext * dust_1 * rho * 1e-6) / (2 * props%effective_diameter * props%rho)

      call get_optical_props(355, 2, props)
      where (istatus == 0) &
         dust_2 = (3 * props%qext * dust_2 * rho * 1e-6) / (2 * props%effective_diameter * props%rho)

      call get_optical_props(355, 3, props)
      where (istatus == 0) &
         dust_3 = (3 * props%qext * dust_3 * rho * 1e-6) / (2 * props%effective_diameter * props%rho)

      call get_optical_props(355, 4, props)
      where (istatus == 0) &
         dust_4 = (3 * props%qext * dust_4 * rho * 1e-6) / (2 * props%effective_diameter * props%rho)

      call get_optical_props(355, 5, props)
      where (istatus == 0) &
         dust_5 = (3 * props%qext * dust_5 * rho * 1e-6) / (2 * props%effective_diameter * props%rho)

      ! Sum the extinction from each size bin
      where (istatus == 0) &
         extinction = dust_1 + dust_2 + dust_3 + dust_4 + dust_5
      where (istatus /= 0) &
         extinction = MISSING_R8
   end subroutine get_extinction


   ! Sanity check for iostat, if it is not 0, call the error handler
   subroutine check_iostat(istat, routine, varname, msgstring)

      integer,          intent(in) :: istat
      character(len=*), intent(in) :: routine
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: msgstring

      if ( istat /= 0 ) then
         write(string1,*) 'istat should be 0 but is ', istat, ' for ' // varname
         call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=msgstring)
      end if

   end subroutine check_iostat


end module obs_def_atm_lidar_mod
! END DART PREPROCESS MODULE CODE
