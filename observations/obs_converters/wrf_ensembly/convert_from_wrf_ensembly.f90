! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_universal_csv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_from_wrf_ensembly - program that reads CSV files in a specific format
!                               and writes a DART obs_seq file. Implemented so you
!                               can do most of your data wrangling in Python, where I
!                               think it's more convenient. Compatible with the
!                               observation management suite of wrf-ensembly.
!
! The CSV file should have the following columns:
! obs_type, longitude, latitude, vert, year, month, day, hour, minute, second,
! obs_value, obs_error, obs_meta
!
! Where obs_type is the OBS_TYPE string, longitude and latitude are in degrees,
! the date should be in UTC and each component an integer (no leading zeros),
! vert is the vertical coordinate and depends on the obs_type. The obs_meta column
! is only required for some obs_type, can be empty otherwise.
!
!     created May. 2025 Thanasis Georgiou, National Observatory of Athens
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use types_mod, only : r8
  use utilities_mod, only : initialize_utilities, register_module, finalize_utilities, &
    open_file, close_file, error_handler, E_ERR, E_MSG
  use time_manager_mod, only : time_type, set_calendar_type, set_date, &
    operator(>=), increment_time, get_time, &
    operator(-), GREGORIAN, operator(+), print_date
  use obs_sequence_mod, only : obs_sequence_type, obs_type
  use location_mod, only : location_type, VERTISSURFACE, VERTISLEVEL, &
    VERTISHEIGHT, VERTISPRESSURE, VERTISUNDEF
  use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq
  use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
    static_init_obs_sequence, init_obs, write_obs_seq, &
    init_obs_sequence, get_num_obs, &
    set_copy_meta_data, set_qc_meta_data
  use obs_kind_mod, only : AIRSENSE_AOD

  implicit none

  ! version controlled file description for error handling, do not edit
  character(len=256), parameter :: source   = "$URL$"
  character(len=32 ), parameter :: revision = "$Revision$"
  character(len=128), parameter :: revdate  = "$Date$"

  ! Describes one line of the CSV file
  type :: csv_line
    character(len=32) :: obs_type
    real(r8) :: longitude
    real(r8) :: latitude
    real(r8) :: vert
    integer :: year
    integer :: month
    integer :: day
    integer :: hour
    integer :: minute
    integer :: second
    real(r8) :: obs_value
    real(r8) :: obs_error
    character(len=200) :: obs_meta
  end type csv_line
  ! When parsing the CSV file
  integer, parameter :: max_line_length = 1024

  ! Describes a key-value pair in obs metadata
  type :: metadata_pair
    character(len=50) :: key
    character(len=100) :: value
  end type metadata_pair

  !! Local variables
  ! CSV parsing
  integer, parameter :: max_records = 10000
  integer :: num_records ! How many records were read
  type(csv_line) :: observations(max_records), r

  ! I/O Paths & unit numbers
  character(len=255) :: csv_input_path = ""
  character(len=255) :: obsseq_output_path = "obs_seq.csv"
  integer :: funit

  ! For converting to DART date-times
  integer :: oday, osec
  ! I/O status for reading the CSV file
  integer :: rcio

  ! Number of copies and QC variables, maximum size of obs_seq
  integer :: num_copies, num_qc, max_obs

  ! QC flag is always 0
  real(r8) :: qc = 0.0_r8

  ! obs_type and vertical location kind
  integer :: obs_type_id, vert_loc_kind

  ! For initialising the obs_sequence on the first observation
  logical :: first_obs

  ! Observation key for obs_types that need metadata
  ! We use -1 as a sentinel value for "no key", for those that don't need them.
  integer :: obs_key = -1

  ! The obs_sequence object from DART
  type(obs_sequence_type) :: obs_seq
  type(obs_type)          :: obs, prev_obs
  type(time_type)         :: time_obs, prev_time

  ! For looping
  integer :: i

  !! Executable bit
  ! Ceremony
  call initialize_utilities('convert_universal_csv')
  call register_module(source, revision, revdate)

  ! Handle arguments
  if (command_argument_count() .ge. 1) then
    call get_command_argument(1, value=csv_input_path)
    if (command_argument_count() .ge. 2) then
      call get_command_argument(2, value=obsseq_output_path)
    end if
  else
    ! No arguments provided - use stdin for input
    csv_input_path = ""
  end if

  ! Check for help flag
  if (command_argument_count() .ge. 1) then
    if (trim(csv_input_path) == '-h' .or. trim(csv_input_path) == '--help') then
      write(*,*) 'Usage: convert_universal_csv [input_file] [output_file]'
      write(*,*) 'input_file: CSV file containing the columns (obs_type, longitude, latitude, vert, year, month, day, hour, minute, second, obs_value, obs_error, obs_meta) IN THIS ORDER.'
      write(*,*) '           If omitted or "-", reads from stdin. First row is assumed to be the header.'
      write(*,*) 'output_file: Output obs_seq file. If omitted, defaults to obs_seq.csv'
      write(*,*) ''
      write(*,*) 'Examples:'
      write(*,*) '  convert_universal_csv data.csv obs_seq.out'
      write(*,*) '  cat data.csv | convert_universal_csv'
      write(*,*) '  convert_universal_csv - obs_seq.out < data.csv'
      stop
    end if
  end if

  ! Handle stdin input case
  if (trim(csv_input_path) == "-") then
    csv_input_path = ""
  end if

  ! Time setup
  call set_calendar_type(GREGORIAN)

  ! Initialise the observation sequence in memory
  num_copies = 1 ! One data value
  num_qc = 1 ! One QC value
  max_obs = 100000 ! Max 100000 obs

  call static_init_obs_sequence()
  call init_obs(obs, num_copies, num_qc)
  call init_obs(prev_obs, num_copies, num_qc)
  first_obs = .true.
  call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)
  call set_copy_meta_data(obs_seq, 1, 'observation')
  call set_qc_meta_data(obs_seq, 1, 'Data QC')

  ! Read the CSV file
  call read_csv_file(csv_input_path, observations, num_records, max_records)

  ! Loop over all CSV lines, creating observations
  obsloop: do i = 1, num_records
    r = observations(i)

    call map_observation_type(trim(r%obs_type), obs_type_id, vert_loc_kind)

    ! Some observations need metadata, now it's the time to handle this
    call handle_obs_metadata(obs_type_id, trim(r%obs_meta), obs_key)

    ! Set the time
    time_obs = set_date(r%year, r%month, r%day, r%hour, r%minute, r%second)
    call get_time(time_obs, osec, oday)

    ! Set the observation value and error
    call create_3d_obs(r%latitude, r%longitude, r%vert, vert_loc_kind, r%obs_value, obs_type_id, r%obs_error, oday, osec, qc, obs, obs_key)
    call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    first_obs = .false.
  end do obsloop
  write(*,*) 'Number of observations in obs_seq: ', get_num_obs(obs_seq)

  ! Write the output file
  if (get_num_obs(obs_seq) > 0) then
    print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
    call write_obs_seq(obs_seq, obsseq_output_path)
  endif

  call finalize_utilities()

contains

  ! Parses a line from the CSV file into `data`
  ! Also does some basic validation on the fields
  subroutine parse_csv_line(line, data, is_valid)
    character(len=*), intent(in) :: line
    type(csv_line), intent(out) :: data
    logical, intent(out) :: is_valid

    character(len=max_line_length) :: fields(13)
    integer :: field_count, ios


    is_valid = .false.
    call split_csv_line(line, fields, field_count)

    ! Check if we have the expected number of fields (12 or 13, allowing empty
    ! last field because some observations don't have metadata)
    if (field_count < 12 .or. field_count > 13) then
      write(*,*) 'Error: Expected 12-13 fields, found ', field_count
      return
    endif

    ! Field 1: obs_type (string)
    data%obs_type = trim(adjustl(fields(1)))
    if (len_trim(data%obs_type) == 0) then
      write(*,*) 'Error: Empty obs_type field'
      return
    endif

    ! Field 2: longitude (real)
    read(fields(2), *, iostat=ios) data%longitude
    if (ios /= 0) then
      write(*,*) 'Error: Invalid longitude: ', trim(fields(2))
      return
    endif
    ! Ensure longitude is in the range [0, 360]
    data%longitude = mod(data%longitude, 360.0_r8)
    if (data%longitude < 0.0_r8) data%longitude = data%longitude + 360.0_r8
    if (data%longitude < 0.0_r8 .or. data%longitude > 360.0_r8) then
      write(*,*) 'Error: Longitude out of range: ', data%longitude
      return
    endif

    ! Field 3: latitude (real)
    read(fields(3), *, iostat=ios) data%latitude
    if (ios /= 0) then
      write(*,*) 'Error: Invalid latitude: ', trim(fields(3))
      return
    endif
    if (abs(data%latitude) > 90.0_r8) then
      write(*,*) 'Error: Latitude out of range: ', data%latitude
      return
    endif

    ! Field 4: vert (real)
    read(fields(4), *, iostat=ios) data%vert
    if (ios /= 0) then
      write(*,*) 'Error: Invalid vertical coordinate: ', trim(fields(4))
      return
    endif

    ! Field 5: year (integer)
    read(fields(5), *, iostat=ios) data%year
    if (ios /= 0) then
      write(*,*) 'Error: Invalid year: ', trim(fields(5))
      return
    endif
    if (data%year < 1900 .or. data%year > 2100) then
      write(*,*) 'Error: Year out of reasonable range: ', data%year
      return
    endif

    ! Field 6: month (integer)
    read(fields(6), *, iostat=ios) data%month
    if (ios /= 0) then
      write(*,*) 'Error: Invalid month: ', trim(fields(6))
      return
    endif
    if (data%month < 1 .or. data%month > 12) then
      write(*,*) 'Error: Month out of range: ', data%month
      return
    endif

    ! Field 7: day (integer)
    read(fields(7), *, iostat=ios) data%day
    if (ios /= 0) then
      write(*,*) 'Error: Invalid day: ', trim(fields(7))
      return
    endif
    if (data%day < 1 .or. data%day > 31) then
      write(*,*) 'Error: Day out of range: ', data%day
      return
    endif

    ! Field 8: hour (integer)
    read(fields(8), *, iostat=ios) data%hour
    if (ios /= 0) then
      write(*,*) 'Error: Invalid hour: ', trim(fields(8))
      return
    endif
    if (data%hour < 0 .or. data%hour > 23) then
      write(*,*) 'Error: Hour out of range: ', data%hour
      return
    endif

    ! Field 9: minute (integer)
    read(fields(9), *, iostat=ios) data%minute
    if (ios /= 0) then
      write(*,*) 'Error: Invalid minute: ', trim(fields(9))
      return
    endif
    if (data%minute < 0 .or. data%minute > 59) then
      write(*,*) 'Error: Minute out of range: ', data%minute
      return
    endif

    ! Field 10: second (integer)
    read(fields(10), *, iostat=ios) data%second
    if (ios /= 0) then
      write(*,*) 'Error: Invalid second: ', trim(fields(10))
      return
    endif
    if (data%second < 0 .or. data%second > 59) then
      write(*,*) 'Error: Second out of range: ', data%second
      return
    endif

    ! Field 11: obs_value (real)
    read(fields(11), *, iostat=ios) data%obs_value
    if (ios /= 0) then
      write(*,*) 'Error: Invalid observation value: ', trim(fields(11))
      return
    endif

    ! Field 12: obs_error (real)
    read(fields(12), *, iostat=ios) data%obs_error
    if (ios /= 0) then
      write(*,*) 'Error: Invalid observation error: ', trim(fields(12))
      return
    endif
    if (data%obs_error < 0.0_r8) then
      write(*,*) 'Error: Observation error cannot be negative: ', data%obs_error
      return
    endif

    ! Field 13: obs_meta (string) - can be empty
    if (field_count >= 13) then
      data%obs_meta = trim(adjustl(fields(13)))
    else
      data%obs_meta = ''  ! Empty if field not present
    endif

    ! If we get here, all fields were parsed successfully
    is_valid = .true.

  end subroutine parse_csv_line

  ! Splits a CSV line into fields by finding the commas
  subroutine split_csv_line(line, fields, field_count)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: fields(:)
    integer, intent(out) :: field_count

    integer :: i, start_pos, end_pos, line_len
    logical :: in_field

    field_count = 0
    line_len = len_trim(line)
    start_pos = 1
    in_field = .false.

    ! Handle empty line
    if (line_len == 0) return

    do i = 1, line_len + 1
      if (i <= line_len) then
        if (line(i:i) == ',') then
          ! Found delimiter
          if (in_field) then
            end_pos = i - 1
            field_count = field_count + 1
            if (field_count <= size(fields)) then
              if (end_pos >= start_pos) then
                fields(field_count) = line(start_pos:end_pos)
              else
                fields(field_count) = ''
              endif
            endif
            in_field = .false.
          else
            ! Empty field
            field_count = field_count + 1
            if (field_count <= size(fields)) then
              fields(field_count) = ''
            endif
          endif
          start_pos = i + 1
        else
          ! Regular character
          if (.not. in_field) then
            in_field = .true.
            start_pos = i
          endif
        endif
      else
        ! End of line
        if (in_field) then
          end_pos = line_len
          field_count = field_count + 1
          if (field_count <= size(fields)) then
            fields(field_count) = line(start_pos:end_pos)
          endif
        endif
      endif
    enddo

  end subroutine split_csv_line

  ! Reads a CSV file line by line and populates the data_array with the parsed
  ! observations. The number of records read is returned in num_records.
  ! The maximum number of records to read is specified by max_records.
  ! If filename is empty, reads from stdin instead.
  subroutine read_csv_file(filename, data_array, num_records, max_records)
    character(len=*), intent(in) :: filename
    type(csv_line), intent(out) :: data_array(:)
    integer, intent(out) :: num_records
    integer, intent(in) :: max_records

    integer :: unit_num, ios, line_num
    character(len=max_line_length) :: line
    type(csv_line) :: temp_data
    logical :: is_valid, use_stdin

    num_records = 0
    line_num = 0
    use_stdin = (len_trim(filename) == 0)

    if (use_stdin) then
      ! Use stdin (unit 5)
      unit_num = 5
      write(*,*) 'Reading from stdin...'
      ios = 0  ! Initialize ios for stdin reading
    else
      ! Open file
      open(newunit=unit_num, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
        write(*,*) 'Error: Cannot open file: ', trim(filename)
        return
      endif
      write(*,*) 'Reading from file: ', trim(filename)
    endif

    ! Read file/stdin line by line
    do while (ios == 0 .and. num_records < max_records)
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit

      line_num = line_num + 1

      ! Skip empty lines
      if (len_trim(line) == 0) cycle

      ! Skip header line (first line)
      if (line_num == 1) then
        write(*,*) 'Skipping header line: ', trim(line)
        cycle
      endif

      ! Parse the line
      call parse_csv_line(line, temp_data, is_valid)

      if (is_valid) then
        num_records = num_records + 1
        data_array(num_records) = temp_data
      else
        write(*,*) 'Warning: Skipping invalid line ', line_num, ': ', trim(line)
      endif

    enddo

    if (.not. use_stdin) then
      close(unit_num)
    endif

    write(*,*) 'Successfully read ', num_records, ' records from ', line_num, ' lines'

  end subroutine read_csv_file

  ! Maps the string representation of observation types from the CSV file to the
  ! corresponding DART obs_type and vertical type.
  subroutine map_observation_type(string_obs_type, dart_obs_type, dart_vert_type)
    use location_mod, only : VERTISUNDEF, VERTISHEIGHT
    use obs_kind_mod
    implicit none

    character(len=*), intent(in) :: string_obs_type
    integer, intent(out) :: dart_obs_type, dart_vert_type

    ! Map the string to the corresponding DART obs_type
    select case (trim(string_obs_type))
     case ('AIRSENSE_AOD')
      dart_obs_type = AIRSENSE_AOD
      dart_vert_type = VERTISUNDEF
     case ('EC_ATL_EXT355')
      dart_obs_type = LIDAR_EXTINCTION
      dart_vert_type = VERTISHEIGHT
     case ('HLOS_WIND')
      dart_obs_type = AEOLUS_L2B_HLOS
      dart_vert_type = VERTISHEIGHT
     case default
      dart_obs_type = -1 ! Error value
      dart_vert_type = VERTISUNDEF
      write(*,*) 'Unknown observation type: ', string_obs_type
      stop
    end select
  end subroutine map_observation_type

  ! For obs_types that require metadata (such as AEOLUS HLOS winds), this
  ! subroutine handles parsing the metadata string and using the appropriate
  ! "set metadata" routine. The format of the metadata string is vertical bar-separated
  ! key-value pairs, e.g. "key1=value1|key2=value2".
  ! If the obs_type does not require metadata, this subroutine does nothing.
  ! The obs_key is returned as -1 if no metadata was set, or the actual
  ! observation key if metadata was set.
  subroutine handle_obs_metadata(obs_type_id, obs_meta_str, obs_key)
    use obs_kind_mod, only : AEOLUS_L2B_HLOS
    use obs_def_aeolus_hlos_mod, only: set_aeolus_metadata
    implicit none

    integer, intent(in) :: obs_type_id
    character(len=*), intent(in) :: obs_meta_str
    integer, intent(out) :: obs_key

    ! Max 10 metadata pairs
    type(metadata_pair) :: meta_pairs(10)
    integer :: num_meta_pairs

    ! Storage for the various fields we are searching, might get messy when
    ! we implement more metadata-requiring obs_types
    character(len=100) :: azimuth_str, windresultid_str
    logical :: azimuth_found, windresultid_found
    real(r8) :: azimuth
    integer :: windresultid

    call parse_metadata_string(trim(obs_meta_str), meta_pairs, num_meta_pairs)

    if (obs_type_id == AEOLUS_L2B_HLOS) then
      ! Search for azimuth and wind_result_id in metadata
      call get_metadata_value(meta_pairs, num_meta_pairs, 'azimuth', azimuth_str, azimuth_found)
      call get_metadata_value(meta_pairs, num_meta_pairs, 'wind_result_id', windresultid_str, windresultid_found)

      if (.not. azimuth_found) then
        write(*,*) 'Error: Missing required metadata "azimuth" for AEOLUS_L2B_HLOS observation'
        stop
      endif
      if (.not. windresultid_found) then
        write(*,*) 'Error: Missing required metadata "wind_result_id" for AEOLUS_L2B_HLOS observation'
        stop
      endif

      ! Remember to convert the azimuth to real and wind_result_id to integer
      read(azimuth_str, *) azimuth
      read(windresultid_str, *) windresultid
      print*, 'Setting AEOLUS metadata: wind_result_id=', windresultid, ' azimuth=', azimuth
      call set_aeolus_metadata(obs_key, windresultid, azimuth)

      return
    end if

    ! If we reach here, no metadata was set
    obs_key = -1
  end subroutine handle_obs_metadata

  subroutine parse_metadata_string(meta_str, metadata, num_metadata)
    character(len=*), intent(in) :: meta_str
    type(metadata_pair), intent(out) :: metadata(:)
    integer, intent(out) :: num_metadata

    character(len=200) :: work_str
    integer :: pos, next_pos, eq_pos

    num_metadata = 0
    work_str = trim(meta_str) // '|'  ! Append vertical bar for easier parsing
    pos = 1

    do while (pos < len_trim(work_str) .and. num_metadata < size(metadata))
      next_pos = index(work_str(pos:), '|')
      if (next_pos == 0) exit
      next_pos = pos + next_pos - 1

      eq_pos = index(work_str(pos:next_pos-1), '=')
      if (eq_pos > 0) then
        num_metadata = num_metadata + 1
        metadata(num_metadata)%key = adjustl(trim(work_str(pos:pos+eq_pos-2)))
        metadata(num_metadata)%value = adjustl(trim(work_str(pos+eq_pos:next_pos-1)))
      endif

      pos = next_pos + 1
    enddo
  end subroutine

! Helper subroutine to get a value by key
  subroutine get_metadata_value(metadata, num_metadata, key, value, found)
    type(metadata_pair), intent(in) :: metadata(:)
    integer, intent(in) :: num_metadata
    character(len=*), intent(in) :: key
    character(len=*), intent(out) :: value
    logical, intent(out) :: found
    integer :: i

    found = .false.
    do i = 1, num_metadata
      if (trim(metadata(i)%key) == trim(key)) then
        value = trim(metadata(i)%value)
        found = .true.
        return
      endif
    enddo
  end subroutine
end program convert_universal_csv

! ! <next few lines under version control, do not edit>
! ! $URL$
! ! $Id$
! ! $Revision$
! ! $Date$
