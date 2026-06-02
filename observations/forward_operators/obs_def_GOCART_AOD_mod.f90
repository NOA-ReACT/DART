! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! AIRSENSE_AOD,         QTY_AOD
! AIRSENSE_AOD_FINE,    QTY_AOD
! AIRSENSE_AOD_COARSE,  QTY_AOD
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_GOCART_AOD_mod, only : get_aod, &
!      AOD_MODE_TOTAL, AOD_MODE_FINE, AOD_MODE_COARSE
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(AIRSENSE_AOD)
!         call get_aod(AOD_MODE_TOTAL, state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!      case(AIRSENSE_AOD_FINE)
!         call get_aod(AOD_MODE_FINE, state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!      case(AIRSENSE_AOD_COARSE)
!         call get_aod(AOD_MODE_COARSE, state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(AIRSENSE_AOD, AIRSENSE_AOD_FINE, AIRSENSE_AOD_COARSE)
!      continue
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(AIRSENSE_AOD, AIRSENSE_AOD_FINE, AIRSENSE_AOD_COARSE)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(AIRSENSE_AOD, AIRSENSE_AOD_FINE, AIRSENSE_AOD_COARSE)
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
      QTY_GC_DUST_BIN4, QTY_GC_DUST_BIN5, &
      QTY_GC_SEAS_BIN1, QTY_GC_SEAS_BIN2, QTY_GC_SEAS_BIN3, &
      QTY_GC_SEAS_BIN4
   use  assim_model_mod, only : interpolate
   use obs_def_utilities_mod, only : track_status
   use ensemble_manager_mod,  only : ensemble_type

   public :: get_aod, AOD_MODE_TOTAL, AOD_MODE_FINE, AOD_MODE_COARSE

   ! AOD modes: which aerosol size bins to integrate for the AOD
   !   TOTAL  - all dust and sea salt bins
   !   FINE   - DUST_1 and SEAS_1, SEAS_2
   !   COARSE - DUST_2..DUST_5 and SEAS_3, SEAS_4
   integer, parameter :: AOD_MODE_TOTAL  = 0
   integer, parameter :: AOD_MODE_FINE   = 1
   integer, parameter :: AOD_MODE_COARSE = 2

! version controlled file description for error handling, do not edit
   character(len=256), parameter :: source   = "$URL$"
   character(len=32 ), parameter :: revision = "$Revision$"
   character(len=128), parameter :: revdate  = "$Date$"

   character(len=256) :: string1, string2
   logical, save      :: module_initialized = .false.

   ! For loops
   integer :: i, obs_kind

   ! Represents a row in the optical properties CSV file
   type :: optical_properties_csv
      integer :: wavelength ! Wavelength of the laser [nm]
      character(len=6) :: species ! Species name
      real(r8) :: ext_eff ! Extinction efficiency [meter ** 2 / kilogram]
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

         ! Copy to props
         props%wavelength = line%wavelength
         props%species = line%species
         props%ext_eff = line%ext_eff

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
         else if (line%species == 'SEAS_1') then
            props%qty = QTY_GC_SEAS_BIN1
         else if (line%species == 'SEAS_2') then
            props%qty = QTY_GC_SEAS_BIN2
         else if (line%species == 'SEAS_3') then
            props%qty = QTY_GC_SEAS_BIN3
         else if (line%species == 'SEAS_4') then
            props%qty = QTY_GC_SEAS_BIN4
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
            write(*, *) 'Ext_eff: ', optical_props(i)%ext_eff
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
   ! `aod_mode` selects which aerosol size bins contribute (see AOD_MODE_* parameters)
   subroutine get_aod(aod_mode, state_handle, ens_size, location, key, aod, istatus)
      integer, intent(in) :: aod_mode
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

      ! Which observation kinds to use for AOD (selected below based on aod_mode)
      integer, allocatable :: obs_kinds(:)

      ! Storage for model variables
      real(r8), allocatable :: concentrations(:, :, :), extinctions(:, :, :) ! (ens_size, model_levels, size(obs_kinds))
      real(r8), allocatable :: rho(:, :) ! (ens_size, model_levels)
      real(r8) :: current_level_height(ens_size)
      real(r8) :: current_level_extinction(ens_size)

      ! Initialize the module
      call initialize_module()

      ! Select which aerosol size bins contribute to the AOD
      select case (aod_mode)
      case (AOD_MODE_TOTAL)
         obs_kinds = [QTY_GC_DUST_BIN1, QTY_GC_DUST_BIN2, QTY_GC_DUST_BIN3, QTY_GC_DUST_BIN4, QTY_GC_DUST_BIN5, &
                      QTY_GC_SEAS_BIN1, QTY_GC_SEAS_BIN2, QTY_GC_SEAS_BIN3, QTY_GC_SEAS_BIN4]
      case (AOD_MODE_FINE)
         obs_kinds = [QTY_GC_DUST_BIN1, QTY_GC_SEAS_BIN1, QTY_GC_SEAS_BIN2]
      case (AOD_MODE_COARSE)
         obs_kinds = [QTY_GC_DUST_BIN2, QTY_GC_DUST_BIN3, QTY_GC_DUST_BIN4, QTY_GC_DUST_BIN5, &
                      QTY_GC_SEAS_BIN3, QTY_GC_SEAS_BIN4]
      case default
         write(string1, *) 'Unknown aod_mode ', aod_mode
         call error_handler(E_ERR, 'get_aod', string1, source, revision, revdate)
      end select

      ! Determine how many levels the model has
      call determine_model_levels(location, state_handle, ens_size, 500, model_levels, istatus)
      ! print*, 'Model has ', model_levels, ' levels'

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
      ! if (debug) then
      !    ! Print geopotential heights of first member
      !    ! write(*, *) 'Geopotential heights of first member:'
      !    do current_level = 1, model_levels
      !       write(*, *) level_heights(1, current_level)
      !    end do
      ! end if

      ! Get model density for all levels
      ! if (debug) write(*, *) 'Interpolating density'
      do current_level = 1, model_levels
         call make_location_vertislevel(location, vert_location, real(current_level - 1, r8))
         call interpolate(state_handle, ens_size, vert_location, QTY_DENSITY, rho(:, current_level), istatus)
         call track_status(ens_size, istatus, aod, istatus, return_now)
         if (return_now) return
      end do

      ! For each obs_type, get the concentrations and compute the extinctions
      ! if (debug) print *, 'Computing extinctions'
      do obs_kind = 1, size(obs_kinds)
         ! if (debug) print *, 'Computing for obs_kind ', obs_kinds(obs_kind)

         call get_optical_props(532, obs_kinds(obs_kind), props)

         do current_level = 1, model_levels
            ! Get the mixing ratios for this obs_type
            call make_location_vertislevel(location, vert_location, real(current_level - 1, r8))
            call interpolate(state_handle, ens_size, vert_location, obs_kinds(obs_kind), concentrations(:, current_level, obs_kind), this_istatus)
            call track_status(ens_size, this_istatus, concentrations(:, current_level, obs_kind), istatus, return_now)

            if (return_now) return

            ! Multiply with density to get concentrations
            concentrations(:, current_level, obs_kind) = concentrations(:, current_level, obs_kind) * rho(:, current_level) * 1e-9 ! Convert from ug/m^3 to kg/m^3

            ! Use optical properties to compute extinction
            extinctions(:, current_level, obs_kind) = concentrations(:, current_level, obs_kind) * props%ext_eff
         end do
      end do

      ! Print total extinction & concentration for each model level
      ! do current_level = 1, model_levels
      !    print*, 'Level', current_level, ", ", sum(concentrations(:, current_level, :), dim=2) / ens_size, "kg/m^3", ", ", sum(extinctions(:, current_level, :), dim=2) / ens_size, "1/m", "density ", sum(rho(:, current_level)) / ens_size
      ! end do

      ! print *, 'Total concentration for whole column', sum(sum(sum(concentrations(:, :, :), dim=2) / ens_size, dim=1), dim=1), "kg/m^3"

      ! To compute AOD, integrate the extinction over all levels
      aod = 0.0_r8
      do current_level = 1, model_levels - 1
         current_level_height = (level_heights(:, current_level + 1) - level_heights(:, current_level))
         current_level_extinction = (sum(extinctions(:, current_level + 1, :), dim=2) + sum(extinctions(:, current_level, :), dim=2)) / 2
         aod = aod + (current_level_extinction * current_level_height)
         ! print*, "Current level: ", current_level, " height: ", current_level_height, " total_extinction: ", current_level_extinction
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
         write(string1, *) 'Reached max_levels without finding an error'
      end if

      if (has_error) then
         print*, 'Error encountered during level detection'
         print*, 'Levels found: ', levels_found
         print*, 'istatus: ', istatus
         write(string1, *) 'Error in vertical level detection'
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

         ! print *, 'height for level ', current_level, ' is ', sum(interp_val) / ens_size
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
