! filepath: /home/thgeorgiou/work/DART_AIRSENSE/observations/obs_converters/earthcare/ec_atl_ebd.f90
! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program ec_atl_ebd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! ec_atl_ebd - program that reads EarthCARE L2 ATL EBD netCDF files
  !              and writes a DART obs_seq file containing particle
  !              extinction coefficient observations.
  !
  ! The netCDF files contain observations in the ScienceData group with:
  ! - particle_extinction_coefficient_355nm (along_track, JSG_height)
  ! - particle_extinction_coefficient_355nm_error (along_track, JSG_height)
  ! - latitude, longitude, time (along_track)
  ! - height (along_track, JSG_height)
  ! - quality_status, simple_classification (along_track, JSG_height)
  !
  ! Filters applied:
  ! - quality_status <= 1 (good/likely good observations)
  ! - simple_classification == 3 (aerosol observations only)
  !
  ! Usage: ec_atl_ebd <input_netcdf_file> <output_obs_seq_file>
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use types_mod, only : r8, digits12, i2
  use utilities_mod, only : initialize_utilities, register_module, finalize_utilities, &
    error_handler, E_ERR, E_MSG, E_WARN
  use time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
    GREGORIAN, get_time
  use obs_sequence_mod, only : obs_sequence_type, obs_type, &
    static_init_obs_sequence, init_obs, init_obs_sequence, &
    set_copy_meta_data, set_qc_meta_data, write_obs_seq, get_num_obs
  use location_mod, only : VERTISHEIGHT
  use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq
  use obs_kind_mod, only : LIDAR_EXTINCTION
  use netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, &
    nc_get_dimension_size, nc_get_variable, nc_get_variable_info
  use netcdf

  implicit none

  ! version controlled file description for error handling, do not edit
  character(len=256), parameter :: source   = "$URL$"
  character(len=32 ), parameter :: revision = "$Revision$"
  character(len=128), parameter :: revdate  = "$Date$"

  ! File paths
  character(len=255) :: input_file = ""
  character(len=255) :: output_file = ""
  character(len=255) :: coord_file = ""

  ! NetCDF file ID and group ID
  integer :: ncid, grp_id

  ! Dimensions
  integer :: nalong, nheight

  ! Data arrays
  real(r8), allocatable :: lat(:), lon(:)
  real(digits12), allocatable :: time_var(:)
  real(r8), allocatable :: height(:,:)
  real(r8), allocatable :: extinction(:,:), extinction_error(:,:)
  integer(i2), allocatable :: quality(:,:), classification(:,:)

  ! DART observation sequence variables
  type(obs_sequence_type) :: obs_seq
  type(obs_type) :: obs, prev_obs
  type(time_type) :: obs_time, prev_time
  integer :: num_copies, num_qc, max_obs
  real(r8) :: qc = 0.0_r8
  logical :: first_obs

  ! Time conversion variables
  integer :: days, seconds, year, month, day, hour, minute, sec
  real(digits12) :: time_origin_days, obs_time_days
  type(time_type) :: time_origin

  ! Loop variables and counters
  integer :: i, j, num_good_obs
  integer :: coord_unit

  logical :: debug = .true.
  integer :: ndims
  integer :: dimlens(10)
  character(len=64) :: dimnames(10)

  ! Initialize DART
  call initialize_utilities('ec_atl_ebd')
  call register_module(source, revision, revdate)

  ! Handle command line arguments
  if (command_argument_count() .ne. 2) then
    write(*,*) 'Usage: ec_atl_ebd <input_netcdf_file> <output_obs_seq_file>'
    write(*,*) 'Converts EarthCARE L2 ATL EBD netCDF files to DART obs_seq format'
    stop
  end if

  call get_command_argument(1, value=input_file)
  call get_command_argument(2, value=output_file)

  ! Create coordinate file name by replacing .out with .coords
  coord_file = trim(output_file)
  if (index(coord_file, '.out') > 0) then
    coord_file = coord_file(1:index(coord_file, '.out')-1) // '.coords'
  else
    coord_file = trim(coord_file) // '.coords'
  endif

  ! Set calendar type
  call set_calendar_type(GREGORIAN)

  ! Open netCDF file and get ScienceData group
  ncid = nc_open_file_readonly(input_file, 'ec_atl_ebd')
  call check_nc_error(nf90_inq_ncid(ncid, 'ScienceData', grp_id), 'finding ScienceData group')

  ! Get dimensions
  nalong = nc_get_dimension_size(grp_id, 'along_track')
  nheight = nc_get_dimension_size(grp_id, 'JSG_height')

  write(*,*) 'Reading ', nalong, ' along-track points with ', nheight, ' height levels'

  ! Allocate arrays
  allocate(lat(nalong), lon(nalong), time_var(nalong))
  allocate(height(nheight, nalong))
  allocate(extinction(nheight, nalong), extinction_error(nheight, nalong))
  allocate(quality(nheight, nalong), classification(nheight, nalong))

  ! Read data from netCDF file
  if (debug) write(*,*) 'Reading latitude'
  call nc_get_variable(grp_id, 'latitude', lat)
  if (debug) write(*,*) 'Reading longitude'
  call nc_get_variable(grp_id, 'longitude', lon)
  if (debug) write(*,*) 'Reading time'
  call nc_get_variable(grp_id, 'time', time_var)
  if (debug) write(*,*) 'Reading height'
  call nc_get_variable(grp_id, 'height', height)
  if (debug) write(*,*) 'Reading extinction coefficient'
  call nc_get_variable(grp_id, 'particle_extinction_coefficient_355nm', extinction)
  if (debug) write(*,*) 'Reading extinction error'
  call nc_get_variable(grp_id, 'particle_extinction_coefficient_355nm_error', extinction_error)
  if (debug) write(*,*) 'Reading quality status'
  call nc_get_variable(grp_id, 'quality_status', quality)
  if (debug) write(*,*) 'Reading simple classification'
  call nc_get_variable(grp_id, 'simple_classification', classification)

  ! Close netCDF file
  call nc_close_file(ncid)

  ! Convert longitude to [0, 360] range
  where(lon < 0.0_r8) lon = lon + 360.0_r8

  ! Set up time origin (2000-01-01 00:00:00)
  time_origin = set_date(2000, 1, 1, 0, 0, 0)
  call get_time(time_origin, seconds, days)
  time_origin_days = days + seconds/86400.0_r8

  ! Initialize observation sequence
  num_copies = 1
  num_qc = 1
  max_obs = nalong * nheight ! Conservative estimate

  call static_init_obs_sequence()
  call init_obs(obs, num_copies, num_qc)
  call init_obs(prev_obs, num_copies, num_qc)
  call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)
  call set_copy_meta_data(obs_seq, 1, 'observation')
  call set_qc_meta_data(obs_seq, 1, 'Data QC')

  first_obs = .true.
  num_good_obs = 0

  ! Open coordinate file for writing
  open(newunit=coord_unit, file=coord_file, status='replace', action='write')
  write(coord_unit, '(A)') 'height_index,along_track_index,longitude,latitude,time'

  ! Process observations
  do i = 1, nheight
    do j = 1, nalong

      ! Skip if missing data (check for fill values)
      if (abs(lat(j)) > 90.0_r8) cycle
      if (abs(lon(j)) > 360.0_r8) cycle
      if (abs(time_var(j)) > 1e30_digits12) cycle
      if (abs(height(i,j)) > 1e30_r8) cycle
      if (abs(extinction(i,j)) > 1e30_r8) cycle
      if (abs(extinction_error(i,j)) > 1e30_r8) cycle

      ! Apply quality filters
      if (quality(i,j) > 0) cycle  ! Only good/likely good observations
      if ((classification(i,j) /= 3) .and. (classification(i,j) /= 0)) cycle  ! Only aerosol observations and clear!

      ! Skip if extinction error is non-positive
      if (extinction_error(i,j) <= 0.0_r8) cycle

      ! Skip too large extinctions, might be clouds
      if (extinction(i,j) > 0.00075_r8) cycle  ! Adjust threshold as needed

      ! Convert time to DART format
      obs_time_days = time_origin_days + time_var(j) / 86400.0_digits12
      days = floor(obs_time_days)
      seconds = nint((obs_time_days - real(days,digits12)) * 86400.0_digits12)
      obs_time = set_time(seconds, days)

      ! Create observation
      call create_3d_obs(lat(j), lon(j), height(i,j), VERTISHEIGHT, &
        extinction(i,j), LIDAR_EXTINCTION, &
        extinction_error(i,j), days, seconds, qc, obs)

      ! Add to sequence
      call add_obs_to_seq(obs_seq, obs, obs_time, prev_obs, prev_time, first_obs)
      first_obs = .false.
      num_good_obs = num_good_obs + 1

      ! Write coordinates to file
      write(coord_unit, '(I0,A,I0,A,F0.6,A,F0.6,A,F0.3)') i, ',', j, ',', lon(j), ',', lat(j), ',', time_var(j)

    end do
  end do

  ! Close coordinate file
  close(coord_unit)

  write(*,*) 'Total observations processed: ', num_good_obs
  write(*,*) 'Number of observations in obs_seq: ', get_num_obs(obs_seq)
  write(*,*) 'Coordinate file written: ', trim(coord_file)

  ! Write output file
  if (get_num_obs(obs_seq) > 0) then
    write(*,*) 'Writing obs_seq file: ', trim(output_file)
    call write_obs_seq(obs_seq, output_file)
  else
    write(*,*) 'No valid observations to write, not writing obs_seq file!'
    ! call error_handler(E_WARN, 'ec_atl_ebd', 'No valid observations found', source)
  endif

  ! Clean up
  deallocate(lat, lon, time_var, height, extinction, extinction_error, quality, classification)
  call finalize_utilities()

contains

  ! Simple wrapper for netCDF error checking
  subroutine check_nc_error(status, context)
    integer, intent(in) :: status
    character(len=*), intent(in) :: context

    if (status /= nf90_noerr) then
      call error_handler(E_ERR, 'ec_atl_ebd', &
        'NetCDF error in '//trim(context)//': '//trim(nf90_strerror(status)), source)
    endif
  end subroutine check_nc_error

end program ec_atl_ebd
