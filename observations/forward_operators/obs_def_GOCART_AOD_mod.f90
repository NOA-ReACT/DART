! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! AIRSENSE_AOD,  QTY_AOD
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_GOCART_AOD_mod, only : get_aod
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(AIRSENSE_AOD)
!         call get_aod(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(AIRSENSE_AOD)
!      continue
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(AIRSENSE_AOD)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(AIRSENSE_AOD)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_GOCART_AOD_mod

   use types_mod, only : r8, PI, metadatalength, MISSING_R8, DEG2RAD
   use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
      logfileunit, get_unit, open_file, close_file, &
      file_exist, ascii_file_format
   use     location_mod, only : location_type, set_location, get_location, &
      VERTISHEIGHT, VERTISLEVEL, VERTISUNDEF, set_location_missing
   use     obs_kind_mod, only : QTY_AOD, QTY_DENSITY, QTY_GEOPOTENTIAL_HEIGHT, &
      QTY_GC_DUST_BIN1, QTY_GC_DUST_BIN2, QTY_GC_DUST_BIN3, &
      QTY_GC_DUST_BIN4, QTY_GC_DUST_BIN5
   use  assim_model_mod, only : interpolate
   use obs_def_utilities_mod, only : track_status
   use ensemble_manager_mod,  only : ensemble_type

   public :: get_aod

! version controlled file description for error handling, do not edit
   character(len=256), parameter :: source   = "$URL$"
   character(len=32 ), parameter :: revision = "$Revision$"
   character(len=128), parameter :: revdate  = "$Date$"

   character(len=256) :: string1, string2
   logical, save      :: module_initialized = .false.

   ! For loops
   integer :: i, obs_kind

   ! Represents a row in the optical properties CSV file
   ! Represents a row in the optical properties CSV file
   type :: optical_properties_csv
      integer :: wavelength ! Wavelength of the laser [nm]
      character(len=6) :: species ! Species name
      real(r8) :: qext ! Extinction efficiency [unitless]
      real(r8) :: rho ! Particle density [kg/m^3]
      real(r8) :: effective_diameter ! Effective diameter [m]
   end type optical_properties_csv

   ! Represents the optical properties of an aerosol bin
   type, extends(optical_properties_csv) :: optical_properties
      integer :: qty ! DART quantity type
   end type optical_properties

   type(optical_properties), allocatable, dimension(:) :: optical_props

   logical :: debug = .TRUE.

contains

   ! Read the optical properties from a CSV file
   subroutine initialize_module
      integer :: f, fstat ! File handle for optical properties
      type(optical_properties_csv) :: line
      type(optical_properties) :: props

      if (module_initialized) return
      call register_module(source, revision, revdate)

      open(action='read', file='GOCART_AOD_optical_properties.csv', newunit=f, iostat=fstat)
      if (fstat /= 0) then
         write(string1, *) 'Could not open optical properties file GOCART_AOD_optical_properties.csv'
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
         line%rho = line%rho * 1e6 ! Convert from g/cm^3 to kg/m^3
         line%effective_diameter = line%effective_diameter * 1e-6 ! Convert from um to m

         ! Copy to props
         props%wavelength = line%wavelength
         props%species = line%species
         props%qext = line%qext
         props%rho = line%rho
         props%effective_diameter = line%effective_diameter

         ! Identify the quantity type
         if (line%species == 'DUST_1') then
            props%qty = QTY_GC_DUST_BIN1
         else if (line%species == 'DUST_2') then
            props%qty = QTY_GC_DUST_BIN2
         else if (line%species == 'DUST_3') then
            props%qty = QTY_GC_DUST_BIN3
         else if (line%species == 'DUST_4') then
            props%qty = QTY_GC_DUST_BIN4
         else if (line%species == 'DUST_5') then
            props%qty = QTY_GC_DUST_BIN5
         else
            write(*, *) 'Unknown species ignored: ', line%species
            cycle
         endif

         optical_props = [optical_props, props]
      end do

      print*, 'Read ', size(optical_props), ' optical properties'
      if (size(optical_props) == 0) then
         write(string1, *) 'No optical properties read from file'
         call error_handler(E_ERR, 'initialize_module', string1, source, revision, revdate)
      endif

      ! Print out the optical properties for debugging
      if (debug) then
         do i = 1, size(optical_props)
            write(*, *) 'Wavelength: ', optical_props(i)%wavelength
            write(*, *) 'Species: ', optical_props(i)%species
            write(*, *) 'DART quantity: ', optical_props(i)%qty
            write(*, *) 'Qext: ', optical_props(i)%qext
            write(*, *) 'Rho: ', optical_props(i)%rho
            write(*, *) 'Effective diameter: ', optical_props(i)%effective_diameter
         end do
      endif

      close(f)
      module_initialized = .true.
   end subroutine initialize_module

   ! Returns the optical properties for the given wavelength and DART quantity
   subroutine get_optical_props(wavelength, qty, out)
      integer, intent(in) :: wavelength
      integer, intent(in) :: qty
      type(optical_properties), intent(out) :: out

      do i = 1, size(optical_props)
         if (optical_props(i)%wavelength == wavelength .and. optical_props(i)%qty == qty) then
            out = optical_props(i)
            return
         endif
      end do

      write(string1, *) 'Could not find optical properties for wavelength ', wavelength, ' and size bin ', qty
      call error_handler(E_ERR, 'get_optical_props', string1, source, revision, revdate)
   end subroutine get_optical_props

   ! Forward model for Aerosol Optical Depth (AOD)
   subroutine get_aod(state_handle, ens_size, location, key, aod, istatus)
      type(ensemble_type), intent(in) :: state_handle
      integer, intent(in) :: ens_size
      type(location_type), intent(in) :: location
      integer, intent(in) :: key
      real(r8), intent(inout) :: aod(ens_size)
      integer, intent(out) :: istatus(ens_size)

      ! Optical properties for the current dust bin
      type(optical_properties) :: props

      ! For storing the `interpolator` output
      real(r8) :: interp_val(ens_size)

      ! For vertical level summing
      integer :: model_levels, current_level
      type(location_type) :: vert_location
      real(r8), allocatable :: level_heights(:, :) ! (ens_size, model_levels)

      ! Error tracker
      logical :: return_now = .false.
      integer :: this_istatus(ens_size)

      ! Which observation kinds to use for AOD
      integer :: obs_kinds(5) = [QTY_GC_DUST_BIN1, QTY_GC_DUST_BIN2, QTY_GC_DUST_BIN3, QTY_GC_DUST_BIN4, QTY_GC_DUST_BIN5]

      ! Storage for model variables
      real(r8), allocatable :: concentrations(:, :, :), extinctions(:, :, :) ! (ens_size, model_levels, size(obs_kinds))
      real(r8), allocatable :: rho(:, :) ! (ens_size, model_levels)

      ! Initialize the module
      call initialize_module()

      ! Determine how many levels the model has
      call determine_model_levels(location, state_handle, ens_size, 500, model_levels, istatus)
      print*, 'Model has ', model_levels, ' levels'

      allocate(level_heights(ens_size, model_levels))
      allocate(rho(ens_size, model_levels))
      allocate(concentrations(ens_size, model_levels, size(obs_kinds)))
      allocate(extinctions(ens_size, model_levels, size(obs_kinds)))
      level_heights(:, :) = MISSING_R8
      rho(:, :) = MISSING_R8
      concentrations(:, :, :) = MISSING_R8
      extinctions(:, :, :) = MISSING_R8

      ! Get the geopotential height of each model level
      call get_model_heights(state_handle, ens_size, location, model_levels, level_heights, istatus)
      if (debug) then
         ! Print geopotential heights of first member
         write(*, *) 'Geopotential heights of first member:'
         do current_level = 1, model_levels
            write(*, *) level_heights(1, current_level)
         end do
      end if

      ! Get model density for all levels
      if (debug) write(*, *) 'Interpolating density'
      do current_level = 1, model_levels
         call make_location_vertislevel(location, vert_location, real(current_level - 1, r8))
         call interpolate(state_handle, ens_size, vert_location, QTY_DENSITY, rho(:, current_level), istatus)
         call track_status(ens_size, istatus, aod, istatus, return_now)
         if (return_now) return
      end do

      ! For each obs_type, get the concentrations and compute the extinctions
      if (debug) write(*, *) 'Computing extinctions'
      do obs_kind = 1, size(obs_kinds)
         if (debug) write (*, *) 'Computing for obs_kind ', obs_kinds(obs_kind)

         call get_optical_props(532, obs_kinds(obs_kind), props)

         do current_level = 1, model_levels
            ! Get the mixing ratios for this obs_type
            call make_location_vertislevel(location, vert_location, real(current_level - 1, r8))
            call interpolate(state_handle, ens_size, vert_location, obs_kinds(obs_kind), concentrations(:, current_level, obs_kind), istatus)
            call track_status(ens_size, this_istatus, concentrations(:, current_level, obs_kind), istatus, return_now)
            if (return_now) return

            ! Multiply with density to get concentrations
            concentrations(:, current_level, obs_kind) = concentrations(:, current_level, obs_kind) * rho(:, current_level) * 1e-6 ! Convert from ug/m^3 to g/m^3

            ! Use optical properties to compute extinction
            extinctions(:, current_level, obs_kind) = (3 * concentrations(:, current_level, obs_kind) * props%qext) / (2 * props%rho * props%effective_diameter)
         end do
      end do

      ! To compute AOD, integrate the extinction over all levels
      aod = 0.0_r8
      do current_level = 1, model_levels
         aod = aod + sum(extinctions(:, current_level, :), dim=2) * (level_heights(:, current_level + 1) - level_heights(:, current_level))
      end do
      if (debug) write(*, *) 'AOD: ', aod
   end subroutine get_aod

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

   ! Finds out how many vertical levels the model has
   ! This is done by calling interpolate() with VERTISLEVEL locations until it errors out.
   ! If the level count exceeds `max_levels`, that is the returned value.
   ! Inspiration: https://docs.dart.ucar.edu/en/latest/guide/forward_operator.html#summing-or-integrating-a-column-value
   ! istatus values:
   ! 0: No error
   ! 999: Max levels reached without finding model top
   ! The rest of the values are passed from the interpolate() call
   subroutine determine_model_levels(location, state_handle, ens_size, max_levels, levels_found, istatus)
      type(location_type), intent(in) :: location
      type(ensemble_type), intent(in) :: state_handle
      integer, intent(in) :: ens_size
      integer, intent(in) :: max_levels
      integer, intent(out) :: levels_found
      integer, intent(out) :: istatus(ens_size)

      integer :: level
      type(location_type) :: vert_location
      real(r8), dimension(ens_size) :: interp_val
      logical :: has_error

      levels_found = 0
      has_error = .false.

      do level = 0, max_levels - 1
         call make_location_vertislevel(location, vert_location, real(level, r8))
         call interpolate(state_handle, ens_size, vert_location, QTY_GC_DUST_BIN1, interp_val, istatus)

         if (all(istatus == 0)) then
            ! Valid interpolation at this level
            levels_found = levels_found + 1
         else if (all(istatus == 2)) then
            ! Reached the end of available levels
            exit
         else
            ! Encountered an error
            has_error = .true.
            exit
         end if
      end do

      ! Write a warning if max_levels is reached
      if (levels_found == max_levels) then
         write(string1, *) 'Reached max_levels without finding an error', istatus
      end if

      if (has_error) then
         write(string1, *) 'Error in vertical level detection', istatus
      end if
   end subroutine determine_model_levels

   ! Get the geopotential height of each model level, for each ensemble member.
   ! This is done by interpolating the geopotential height at the given location and for all model levels.
   ! The return array (`level_heights`) should be of dimenions (ens_size, model_levels).
   subroutine get_model_heights(state_handle, ens_size, location, model_levels, level_heights, istatus)
      type(ensemble_type), intent(in) :: state_handle
      integer, intent(in) :: ens_size
      type(location_type), intent(in) :: location
      integer, intent(in) :: model_levels
      real(r8), intent(out) :: level_heights(ens_size, model_levels)
      integer, intent(out) :: istatus(ens_size)

      integer :: current_level
      type(location_type) :: vert_location
      real(r8), dimension(ens_size) :: interp_val

      do current_level = 1, model_levels
         call make_location_vertislevel(location, vert_location, real(current_level - 1, r8))
         call interpolate(state_handle, ens_size, vert_location, QTY_GEOPOTENTIAL_HEIGHT, interp_val, istatus)

         if (any(istatus /= 0)) then
            write(string1, *) 'Error in geopotential height interpolation'
            call error_handler(E_ERR, 'get_model_heights', string1, source, revision, revdate)
         end if

         print *, 'height for level ', current_level, ' is ', sum(interp_val) / ens_size
         level_heights(:, current_level) = interp_val
      end do
   end subroutine get_model_heights

   ! Sets `loc_out` to the same latitude and longitude as `loc_in`, but at the given vertical level.
   subroutine make_location_vertislevel(loc_in, loc_out, level)
      type(location_type), intent(in) :: loc_in
      type(location_type), intent(out) :: loc_out
      real(r8), intent(in) :: level

      real(r8) :: loc_contents(3)

      loc_contents = get_location(loc_in)
      loc_out = set_location(loc_contents(1), loc_contents(2), level, VERTISLEVEL)
   end subroutine make_location_vertislevel
end module
! END DART PREPROCESS MODULE CODE
