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

   use types_mod, only : r8, MISSING_R8
   use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
      find_namelist_in_file, check_namelist_read, NAMELIST_NOT_PRESENT, &
      do_nml_file, do_nml_term, nmlfileunit
   use location_mod, only : location_type, set_location, get_location, VERTISLEVEL
   use obs_kind_mod, only : QTY_DENSITY, QTY_GEOPOTENTIAL_HEIGHT, QTY_GC_DUST_BIN1
   use assim_model_mod, only : interpolate
   use obs_def_utilities_mod, only : track_status
   use ensemble_manager_mod, only : ensemble_type
   use gocart_optics_mod, only : gocart_optics_init, N_GOCART_BINS, GOCART_BIN_QTYS, &
      bin_optics_type, get_bin_optics, bin_ext_eff, bin_fine_fraction, &
      optics_are_hygroscopic, optics_have_fine_fraction, get_expected_growth_factor

   implicit none
   private

   public :: get_aod, AOD_MODE_TOTAL, AOD_MODE_FINE, AOD_MODE_COARSE

   ! AOD modes: how much of each aerosol bin's extinction contributes to the AOD.
   !   TOTAL  - all of it
   !   FINE   - the share of it coming from particles below the ambient-size cutoff the optical
   !            properties table was built with, i.e. the table's fine_fraction[1]
   !   COARSE - the remainder
   ! Every bin contributes to every mode: the fine fraction of a bin is generally neither 0 nor
   ! 1, and for the hygroscopic species it follows the ambient growth of the particles.
   integer, parameter :: AOD_MODE_TOTAL  = 0
   integer, parameter :: AOD_MODE_FINE   = 1
   integer, parameter :: AOD_MODE_COARSE = 2

   ! Give up looking for the model top after this many levels
   integer, parameter :: MAX_MODEL_LEVELS = 500

! version controlled file description for error handling, do not edit
   character(len=256), parameter :: source   = "$URL$"
   character(len=32 ), parameter :: revision = "$Revision$"
   character(len=128), parameter :: revdate  = "$Date$"

   character(len=512) :: string1
   logical, save      :: module_initialized = .false.

   ! Namelist
   integer :: wavelength = 532   ! [nm], has to be one of the wavelengths of the optics table
   logical :: debug      = .false.

   namelist /obs_def_GOCART_AOD_nml/ wavelength, debug

contains

   ! Reads the namelist and the optical properties table
   subroutine initialize_module()
      integer :: iunit, io

      if (module_initialized) return
      module_initialized = .true.

      call register_module(source, revision, revdate)

      call find_namelist_in_file('input.nml', 'obs_def_GOCART_AOD_nml', iunit, optional_nml = .true.)
      if (iunit /= NAMELIST_NOT_PRESENT) then
         read(iunit, nml = obs_def_GOCART_AOD_nml, iostat = io)
         call check_namelist_read(iunit, io, 'obs_def_GOCART_AOD_nml')
      endif
      if (do_nml_file()) write(nmlfileunit, nml = obs_def_GOCART_AOD_nml)
      if (do_nml_term()) write(     *     , nml = obs_def_GOCART_AOD_nml)

      call gocart_optics_init()
   end subroutine initialize_module

   ! Forward model for Aerosol Optical Depth (AOD)
   ! `aod_mode` selects how much of each bin contributes (see the AOD_MODE_* parameters)
   subroutine get_aod(aod_mode, state_handle, ens_size, location, key, aod, istatus)
      integer, intent(in) :: aod_mode
      type(ensemble_type), intent(in) :: state_handle
      integer, intent(in) :: ens_size
      type(location_type), intent(in) :: location
      integer, intent(in) :: key
      real(r8), intent(inout) :: aod(ens_size)
      integer, intent(out) :: istatus(ens_size)

      ! Optical properties of the current aerosol bin
      type(bin_optics_type) :: props
      type(location_type) :: vert_location

      integer :: model_levels, current_level, bin
      logical :: return_now
      integer :: this_istatus(ens_size)

      ! Storage for model variables, all (ens_size, model_levels)
      real(r8), allocatable :: level_heights(:, :), rho(:, :), growth_factor(:, :), &
                               extinctions(:, :)

      real(r8), dimension(ens_size) :: concentration, ext_eff, ext, fine_fraction, &
                                       layer_thickness, layer_extinction

      call initialize_module()

      select case (aod_mode)
      case (AOD_MODE_TOTAL)
         continue
      case (AOD_MODE_FINE, AOD_MODE_COARSE)
         if (.not. optics_have_fine_fraction()) then
            write(string1, *) 'The fine and coarse mode AODs need the fine_fraction[1] column ' // &
                              'of the optical properties table, which this one does not have'
            call error_handler(E_ERR, 'get_aod', string1, source, revision, revdate)
         endif
      case default
         write(string1, *) 'Unknown aod_mode ', aod_mode
         call error_handler(E_ERR, 'get_aod', string1, source, revision, revdate)
      end select

      ! Determine how many levels the model has. This leaves the istatus of the interpolation
      ! that ran off the top of the model behind, so zero it before track_status(), which reads
      ! istatus before it writes it.
      call determine_model_levels(location, state_handle, ens_size, MAX_MODEL_LEVELS, &
                                  model_levels, istatus)
      istatus(:) = 0

      allocate(level_heights(ens_size, model_levels), rho(ens_size, model_levels), &
               growth_factor(ens_size, model_levels), extinctions(ens_size, model_levels))
      level_heights(:, :) = MISSING_R8
      rho(:, :) = MISSING_R8
      growth_factor(:, :) = 1.0_r8
      extinctions(:, :) = 0.0_r8

      ! Get the geopotential height of each model level
      call get_model_heights(state_handle, ens_size, location, model_levels, level_heights, &
                             this_istatus)

      ! Get model density for all levels
      do current_level = 1, model_levels
         call make_location_vertislevel(location, vert_location, real(current_level - 1, r8))
         call interpolate(state_handle, ens_size, vert_location, QTY_DENSITY, &
                          rho(:, current_level), this_istatus)
         call track_status(ens_size, this_istatus, aod, istatus, return_now)
         if (return_now) return
      end do

      ! Get the hygroscopic growth of the particles, but only if the optical properties are
      ! resolved over it. A table without a growth factor axis keeps the old, dry behaviour.
      if (optics_are_hygroscopic()) then
         do current_level = 1, model_levels
            call make_location_vertislevel(location, vert_location, real(current_level - 1, r8))
            call get_expected_growth_factor(state_handle, ens_size, vert_location, &
                                            growth_factor(:, current_level), this_istatus)
            call track_status(ens_size, this_istatus, aod, istatus, return_now)
            if (return_now) return
         end do
      endif

      ! Build the extinction profile, one aerosol bin at a time
      do bin = 1, N_GOCART_BINS
         call get_bin_optics(wavelength, GOCART_BIN_QTYS(bin), props)

         do current_level = 1, model_levels
            ! Get the mixing ratio of this bin
            call make_location_vertislevel(location, vert_location, real(current_level - 1, r8))
            call interpolate(state_handle, ens_size, vert_location, GOCART_BIN_QTYS(bin), &
                             concentration, this_istatus)
            call track_status(ens_size, this_istatus, aod, istatus, return_now)
            if (return_now) return

            ! Multiply with density to get a concentration, ug/kg -> kg/m^3, which is what the
            ! mass extinction efficiencies of the table expect
            concentration = concentration * rho(:, current_level) * 1.0e-9_r8

            call bin_ext_eff(props, growth_factor(:, current_level), ext_eff)
            ext = concentration * ext_eff ! m^2/kg * kg/m^3 -> 1/m

            ! Attribute this bin's extinction to the requested size mode. The fine and the
            ! coarse mode stay exactly additive by taking the coarse part as the remainder.
            if (aod_mode /= AOD_MODE_TOTAL) then
               call bin_fine_fraction(props, growth_factor(:, current_level), fine_fraction)
               if (aod_mode == AOD_MODE_FINE) then
                  ext = fine_fraction * ext
               else
                  ext = ext - fine_fraction * ext
               endif
            endif

            extinctions(:, current_level) = extinctions(:, current_level) + ext
         end do
      end do

      ! To compute AOD, integrate the extinction over all levels
      aod = 0.0_r8
      do current_level = 1, model_levels - 1
         layer_thickness  = level_heights(:, current_level + 1) - level_heights(:, current_level)
         layer_extinction = (extinctions(:, current_level + 1) + &
                             extinctions(:, current_level)) / 2.0_r8
         aod = aod + layer_extinction * layer_thickness
      end do
      where (istatus /= 0) aod = MISSING_R8

      if (debug) write(*, *) 'AOD: ', aod
   end subroutine get_aod

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
