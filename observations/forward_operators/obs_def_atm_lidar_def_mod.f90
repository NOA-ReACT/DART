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
    QTY_GC_DUST_BIN4, QTY_GC_DUST_BIN5, &
    QTY_GC_SEAS_BIN1, QTY_GC_SEAS_BIN2, QTY_GC_SEAS_BIN3, &
    QTY_GC_SEAS_BIN4
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

  logical :: debug = .FALSE.

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

  ! Returns the optical properties for the given wavelength and size bin
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

  ! Computes the extinction from the dust mixing ratios
  subroutine get_extinction(state_handle, ens_size, location, key, extinction, istatus)
    type(ensemble_type), intent(in) :: state_handle
    integer, intent(in) :: ens_size
    type(location_type), intent(in) :: location
    integer, intent(in) :: key
    real(r8), intent(inout) :: extinction(ens_size)
    integer, intent(out) :: istatus(ens_size)

    ! Optical property temp. storage
    type(optical_properties) :: props

    ! For looping
    integer :: obs_kind = 0

    ! Which observation kinds to use for AOD
    integer :: obs_kinds(9) = [QTY_GC_DUST_BIN1, QTY_GC_DUST_BIN2, QTY_GC_DUST_BIN3, QTY_GC_DUST_BIN4, QTY_GC_DUST_BIN5, &
      QTY_GC_SEAS_BIN1, QTY_GC_SEAS_BIN2, QTY_GC_SEAS_BIN3, QTY_GC_SEAS_BIN4]

    ! Storage for model state variables
    real(r8) :: concentrations(ens_size, size(obs_kinds)), &
      extinctions(ens_size, size(obs_kinds)), &
      rho(ens_size)

    ! Error handling
    integer :: this_istatus(ens_size)
    logical :: return_now = .false.

    real(r8) lon, lat, height, obsloc(3)

    ! Initialise
    call initialize_module()
    concentrations(:, :) = MISSING_R8
    extinctions(:, :) = MISSING_R8
    rho(:) = MISSING_R8

    ! Grab airdensity from state
    call interpolate(state_handle, ens_size, location, QTY_DENSITY, rho, this_istatus)
    call track_status(ens_size, this_istatus, extinction, istatus, return_now)
    if (return_now) return

    ! For each size bin, get mixing ratio, compute concentration and then extinction
    do obs_kind = 1, size(obs_kinds)
      call get_optical_props(355, obs_kinds(obs_kind), props)

      call interpolate(state_handle, ens_size, location, obs_kinds(obs_kind), &
        concentrations(:, obs_kind), this_istatus)
      call track_status(ens_size, this_istatus, extinction, istatus, return_now)
      if (return_now) return

      ! Multiply w/ density to get concentration
      concentrations(:, obs_kind) = concentrations(:, obs_kind) * rho(:) * 1e-9  ! Convert from ug/m^3 to kg/m^3
      extinctions(:, obs_kind) = concentrations(:, obs_kind) * props%ext_eff ! 1/m
      extinctions(:, obs_kind) = max(extinctions(:, obs_kind), 0.0_r8) ! No negative extinctions
    end do

    ! Sum the extinction from each size bin
    where (istatus == 0) &
      extinction = sum(extinctions, dim=2)
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
