! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! LIDAR_EXTINCTION_355nm,  QTY_ATM_LIDAR_EXTINCTION
! LIDAR_EXTINCTION_532nm,  QTY_ATM_LIDAR_EXTINCTION
! LIDAR_EXTINCTION_1064nm, QTY_ATM_LIDAR_EXTINCTION
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_atm_lidar_mod, only : get_extinction
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(LIDAR_EXTINCTION_355nm)
!         call get_extinction(355, state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!      case(LIDAR_EXTINCTION_532nm)
!         call get_extinction(532, state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!      case(LIDAR_EXTINCTION_1064nm)
!         call get_extinction(1064, state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(LIDAR_EXTINCTION_355nm)
!      continue
!   case(LIDAR_EXTINCTION_532nm)
!      continue
!   case(LIDAR_EXTINCTION_1064nm)
!      continue
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(LIDAR_EXTINCTION_355nm)
!      continue
!   case(LIDAR_EXTINCTION_532nm)
!      continue
!   case(LIDAR_EXTINCTION_1064nm)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(LIDAR_EXTINCTION_355nm)
!      continue
!   case(LIDAR_EXTINCTION_532nm)
!      continue
!   case(LIDAR_EXTINCTION_1064nm)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_atm_lidar_mod

  use types_mod, only : r8, MISSING_R8
  use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
  use     location_mod, only : location_type
  use     obs_kind_mod, only : QTY_DENSITY
  use  assim_model_mod, only : interpolate
  use obs_def_utilities_mod, only : track_status
  use ensemble_manager_mod,  only : ensemble_type
  use gocart_optics_mod, only : gocart_optics_init, N_GOCART_BINS, GOCART_BIN_QTYS, &
    bin_optics_type, get_bin_optics, bin_ext_eff, optics_are_hygroscopic, &
    get_expected_growth_factor

  implicit none
  private

  public :: get_extinction

! version controlled file description for error handling, do not edit
  character(len=256), parameter :: source   = "$URL$"
  character(len=32 ), parameter :: revision = "$Revision$"
  character(len=128), parameter :: revdate  = "$Date$"

  logical, save :: module_initialized = .false.

contains

  ! Reads the optical properties table shared with the GOCART AOD operator
  subroutine initialize_module()
    if (module_initialized) return
    module_initialized = .true.

    call register_module(source, revision, revdate)
    call gocart_optics_init()
  end subroutine initialize_module

  ! Computes the extinction from the aerosol mixing ratios
  subroutine get_extinction(wavelength, state_handle, ens_size, location, key, extinction, istatus)
    integer, intent(in) :: wavelength
    type(ensemble_type), intent(in) :: state_handle
    integer, intent(in) :: ens_size
    type(location_type), intent(in) :: location
    integer, intent(in) :: key
    real(r8), intent(inout) :: extinction(ens_size)
    integer, intent(out) :: istatus(ens_size)

    ! Optical properties of the current aerosol bin
    type(bin_optics_type) :: props

    integer :: bin

    ! Storage for model state variables
    real(r8), dimension(ens_size) :: rho, growth_factor, concentration, ext_eff, ext

    ! Error handling
    integer :: this_istatus(ens_size)
    logical :: return_now

    call initialize_module()

    ! istatus is intent(out) and is first read inside track_status (`where istatus == 0`),
    ! so it MUST be zeroed here I think.
    istatus(:) = 0
    extinction(:) = 0.0_r8
    growth_factor(:) = 1.0_r8

    ! Grab airdensity from state
    call interpolate(state_handle, ens_size, location, QTY_DENSITY, rho, this_istatus)
    call track_status(ens_size, this_istatus, extinction, istatus, return_now)
    if (return_now) return

    ! Get the hygroscopic growth of the particles, but only if the optical properties are
    ! resolved over it. A table without a growth factor axis keeps the old, dry behaviour.
    if (optics_are_hygroscopic()) then
      call get_expected_growth_factor(state_handle, ens_size, location, growth_factor, this_istatus)
      call track_status(ens_size, this_istatus, extinction, istatus, return_now)
      if (return_now) return
    endif

    ! For each size bin, get mixing ratio, compute concentration and then extinction
    do bin = 1, N_GOCART_BINS
      call get_bin_optics(wavelength, GOCART_BIN_QTYS(bin), props)

      call interpolate(state_handle, ens_size, location, GOCART_BIN_QTYS(bin), &
        concentration, this_istatus)
      call track_status(ens_size, this_istatus, extinction, istatus, return_now)
      if (return_now) return

      ! Multiply w/ density to get concentration
      concentration = concentration * rho * 1.0e-9_r8  ! Convert from ug/m^3 to kg/m^3

      call bin_ext_eff(props, growth_factor, ext_eff)
      ext = max(concentration * ext_eff, 0.0_r8)  ! 1/m, no negative extinctions

      extinction = extinction + ext
    end do

    where (istatus /= 0) extinction = MISSING_R8
  end subroutine get_extinction

end module obs_def_atm_lidar_mod
! END DART PREPROCESS MODULE CODE
