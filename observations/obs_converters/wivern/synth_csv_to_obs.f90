! An observational converter for creating empty obs_seq files based on a CSV file of locations.
! The CSV should include the following columns:
!   observation_type, lat, lon, year, month, day, hour, minute, z, obserr
! The `observation_type` field refers to which kind of observation:
!   RADIOSONDE_DEWPOINT -> RDEWP
!   RADIOSONDE_GEOPOTENTIAL_HGT (hPa) -> RGPHG
!   RADIOSONDE_RELATIVE_HUMIDITY -> RRHUM
!   RADIOSONDE_SURFACE_PRESSURE -> RSPRS
!   RADIOSONDE_TEMPERATURE -> RTEMP
!   RADIOSONDE_U_WIND_COMPONENT -> RUWND
!   RADIOSONDE_V_WIND_COMPONENT -> RVWND
!   SYNOP_DEWPOINT -> SDEWP
!   SYNOP_SPECIFIC_HUMIDITY -> SSPHM
!   SYNOP_SURFACE_PRESSURE -> SSPRS
!   SYNOP_TEMPERATURE -> STEMP
!   SYNOP_U_WIND_COMPONENT -> SUWND
!   SYNOP_V_WIND_COMPONENT -> SVWND
!   QKSWND_U_WIND_COMPONENT -> ASCTU
!   QKSWND_V_WIND_COMPONENT -> ASCTV
! You need to pass the path to the input CSV and where to store the obs_seq as command line arguments.

program synth_csv_to_obs
  use types_mod, only: r8
  use utilities_mod, only: initialize_utilities, finalize_utilities, &
    open_file, close_file, &
    error_handler, E_ERR, E_MSG
  use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
    operator(>=), increment_time, get_time, &
    operator(-), GREGORIAN, operator(+), print_date
  use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
    static_init_obs_sequence, init_obs, write_obs_seq, &
    init_obs_sequence, get_num_obs, &
    set_copy_meta_data, set_qc_meta_data, set_obs_def
  use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
    set_obs_def_location
  use obs_utilities_mod, only : add_obs_to_seq
  use obs_kind_mod, only: RADIOSONDE_DEWPOINT, RADIOSONDE_GEOPOTENTIAL_HGT, &
    RADIOSONDE_RELATIVE_HUMIDITY, RADIOSONDE_SURFACE_PRESSURE, &
    RADIOSONDE_TEMPERATURE, RADIOSONDE_U_WIND_COMPONENT, &
    RADIOSONDE_V_WIND_COMPONENT, SYNOP_DEWPOINT, &
    SYNOP_SPECIFIC_HUMIDITY, SYNOP_SURFACE_PRESSURE, &
    SYNOP_TEMPERATURE, SYNOP_U_WIND_COMPONENT, &
    SYNOP_V_WIND_COMPONENT, QKSWND_U_WIND_COMPONENT, &
    QKSWND_V_WIND_COMPONENT, AIRSENSE_AOD, AEOLUS_L2B_HLOS
  use obs_def_aeolus_hlos_mod, only: set_aeolus_metadata

  implicit none
  ! Type representing a single line of a CSV file.
  type :: csv_line
    character(len=5) :: obs_code
    real(r8) :: lat, lon
    integer :: year, month, day, hour, minute, z
    real(r8) :: observation_error, obs_meta
  end type csv_line
  type(csv_line) :: buf
  integer :: line_obs_type, line_vloc_kind

  ! -- Local variables
  ! I/O paths
  character(len=256) :: path_input_csv = "", path_output_obs_seq = ""

  ! Observation sequence
  type(obs_sequence_type) :: obs_seq
  type(obs_type)          :: obs, prev_obs
  type(time_type)         :: time_obs, prev_time
  integer :: num_copies = 0, num_qc = 0
  logical :: first_obs = .true.
  ! Both are 0 because they will be filled by perfect_model_obs

  ! For observations that have metadata (e.g. AEOLUS), we need a key
  integer :: obs_key = -1
  logical :: obs_key_set = .false.

  ! Time variables
  integer :: osec, oday

  ! File handlers
  integer :: iunit_csv
  integer :: iostat_csv

  ! ---
  ! Start of main

  call initialize_utilities('synth_csv_to_obs')
  call set_calendar_type(GREGORIAN)

  ! Check command line arguments
  if (command_argument_count() /= 2) then
    call error_handler(E_ERR, 'Usage: synth_csv_to_obs <input_csv> <output_obs_seq>', 'main')
  end if
  call get_command_argument(1, path_input_csv)
  call get_command_argument(2, path_output_obs_seq)

  ! Initialise the obs. sequence
  call static_init_obs_sequence()
  call init_obs(obs, num_copies, num_qc)
  call init_obs(prev_obs, num_copies, num_qc)
  call init_obs_sequence(obs_seq, num_copies, num_qc, 200000) ! Max 200000 obs

  ! Open CSV file
  open(action='read', file=path_input_csv, newunit=iunit_csv, iostat=iostat_csv)
  if (iostat_csv /= 0) then
    call error_handler(E_ERR, 'Failed to open input CSV file', 'main')
  end if
  read(iunit_csv, *) ! Skip header

  ! Loop over lines, creating observations
  obsloop: do
    ! Read the line into a buffer
    read(iunit_csv, *, iostat=iostat_csv) buf
    if (iostat_csv /= 0) exit obsloop

    ! Get observation type
    call get_obs_kind(buf%obs_code, line_obs_type, line_vloc_kind)
    if (line_obs_type == -1) then
      call error_handler(E_ERR, 'Unknown observation type', 'main')
      cycle obsloop
    end if

    ! Sanity check location
    if (buf%lon > 180.0_r8 .or. buf%lon < -180.0_r8) cycle obsloop
    if (buf%lon < 0.0_r8) buf%lon = buf%lon + 360.0_r8 ! changes into 0-360

    ! Parse time
    time_obs = set_date(buf%year, buf%month, buf%day, buf%hour, buf%minute, 0) ! 0 seconds
    call get_time(time_obs, osec, oday)
    if (first_obs) then
      ! Print time of first observation for debugging
      call print_date(time_obs)
    end if

    ! Handle metadata for obs. that require it (e.g. AEOLUS)
    if (buf%obs_code == 'ADMT0' .or. buf%obs_code == 'ADMT1') then
      ! The second argument should be the unique observation ID from L2B files
      ! but it's not actually used anywhere. So we just pass -1 because we don't
      ! have access to the IDs in the CSV files.
      call set_aeolus_metadata(obs_key, -1, buf%obs_meta)
      obs_key_set = .true.
    end if

    ! Add to obs. sequence (only pass obs_key if required!)
    if (obs_key_set) then
      call create_obs(buf%lat, buf%lon, line_vloc_kind, buf%z, line_obs_type, &
        oday, osec, buf%observation_error, obs, obs_key)
    else
      call create_obs(buf%lat, buf%lon, line_vloc_kind, buf%z, line_obs_type, &
        oday, osec, buf%observation_error, obs)
    end if
    call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

    first_obs = .false.
    obs_key_set = .false.
  end do obsloop

  ! If we added any observations to the sequence, write the output file now
  if (get_num_obs(obs_seq) == 0) then
    call error_handler(E_ERR, 'No observations added to sequence', 'main')
  end if
  call write_obs_seq(obs_seq, path_output_obs_seq)

  call finalize_utilities()

contains
  ! Get the DART observation kind from the string ID/code we used in the CSV file
  ! Check output for -1, which signifies an unknown observation kind.
  subroutine get_obs_kind(obs_code, observation_kind, vloc_kind)
    use location_mod, only : VERTISPRESSURE, VERTISHEIGHT, VERTISUNDEF
    implicit none

    character(len=*), intent(in) :: obs_code
    integer, intent(out) :: observation_kind
    integer, intent(out) :: vloc_kind

    select case(obs_code)
     case('RDEWP')
      observation_kind = RADIOSONDE_DEWPOINT
      vloc_kind = VERTISPRESSURE
     case('RGPHG')
      observation_kind = RADIOSONDE_GEOPOTENTIAL_HGT
      vloc_kind = VERTISPRESSURE
     case('RRHUM')
      observation_kind = RADIOSONDE_RELATIVE_HUMIDITY
      vloc_kind = VERTISPRESSURE
     case('RSPRS')
      observation_kind = RADIOSONDE_SURFACE_PRESSURE
      vloc_kind = VERTISPRESSURE
     case('RTEMP')
      observation_kind = RADIOSONDE_TEMPERATURE
      vloc_kind = VERTISPRESSURE
     case('RUWND')
      observation_kind = RADIOSONDE_U_WIND_COMPONENT
      vloc_kind = VERTISPRESSURE
     case('RVWND')
      observation_kind = RADIOSONDE_V_WIND_COMPONENT
      vloc_kind = VERTISPRESSURE
     case('SDEWP')
      observation_kind = SYNOP_DEWPOINT
      vloc_kind = VERTISHEIGHT
     case('SSPHM')
      observation_kind = SYNOP_SPECIFIC_HUMIDITY
      vloc_kind = VERTISHEIGHT
     case('SSPRS')
      observation_kind = SYNOP_SURFACE_PRESSURE
      vloc_kind = VERTISHEIGHT
     case('STEMP')
      observation_kind = SYNOP_TEMPERATURE
      vloc_kind = VERTISHEIGHT
     case('SUWND')
      observation_kind = SYNOP_U_WIND_COMPONENT
      vloc_kind = VERTISHEIGHT
     case('SVWND')
      observation_kind = SYNOP_V_WIND_COMPONENT
      vloc_kind = VERTISHEIGHT
     case('ASCTU')
      observation_kind = QKSWND_U_WIND_COMPONENT
      vloc_kind = VERTISHEIGHT
     case('ASCTV')
      observation_kind = QKSWND_V_WIND_COMPONENT
      vloc_kind = VERTISHEIGHT
     case('MODIS')
      observation_kind = AIRSENSE_AOD
      vloc_kind = VERTISUNDEF
     case('ADMT0')
      observation_kind = AEOLUS_L2B_HLOS
      vloc_kind = VERTISHEIGHT
     case('ADMT1')
      observation_kind = AEOLUS_L2B_HLOS
      vloc_kind = VERTISHEIGHT
     case default
      observation_kind = -1  ! Or some error code
      vloc_kind = -1
      print *, 'Unknown observation kind: ', obs_code
    end select
  end subroutine get_obs_kind

  ! Create an observation w/out a value
  subroutine create_obs(lat, lon, vloc_kind, vloc, obstype, day, sec, variance, new_obs, key)
    use types_mod, only : r8
    use obs_def_mod, only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
      set_obs_def_location, set_obs_def_error_variance, set_obs_def_key
    use obs_sequence_mod, only : obs_type, set_obs_def
    use time_manager_mod, only : set_time
    use location_mod, only : set_location

    real(r8), intent(in) :: lat, lon ! Hoz. location
    integer, intent(in) :: vloc_kind ! Vertical location kind
    integer, intent(in) :: vloc ! Vertical location
    integer, intent(in) :: obstype, day, sec ! Obs type and time
    real(r8), intent(in) :: variance ! Observation variance
    type(obs_type), intent(out) :: new_obs ! Observation to fill in
    integer, optional, intent(in) :: key ! Metadata key

    type(obs_def_type) :: obs_def

    call set_obs_def_location(obs_def, set_location(lon, lat, REAL(vloc, r8), vloc_kind))
    call set_obs_def_type_of_obs(obs_def, obstype)
    call set_obs_def_time(obs_def, set_time(sec, day))
    call set_obs_def_error_variance(obs_def, variance)
    if (present(key)) then
      call set_obs_def_key(obs_def, key)
    endif
    call set_obs_def(new_obs, obs_def)
  end subroutine create_obs

end program synth_csv_to_obs
