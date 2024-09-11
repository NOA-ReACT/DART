! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_aeolus_l2b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_aeolus_l2b - program that reads an AEOLUS L2B wind profile
!                        file and writes a DART obs_seq file using the
!                        DART library routines.
!
!     created Sep. 2023 Thanasis Georgiou, National Observatory of Athens
!             Used the COSMOS converters as a guideline
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use     iso_c_binding, only : c_ptr, c_long
   use         types_mod, only : r8, missing_r8
   use     utilities_mod, only : initialize_utilities,  register_module, finalize_utilities
   use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file
   use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
      increment_time, set_time, get_time, GREGORIAN, print_date
   use      location_mod, only : VERTISHEIGHT
   use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
      static_init_obs_sequence, init_obs, write_obs_seq, &
      init_obs_sequence, get_num_obs, &
      set_copy_meta_data, set_qc_meta_data
   use      obs_kind_mod, only : AEOLUS_L2B_HLOS
   use          sort_mod, only : index_sort
   use obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
      create_3d_obs, getvar_int, getdimlen, getvar_real_2d, &
      getvar_int_2d, query_varname
   use obs_def_aeolus_hlos_mod, only: set_aeolus_metadata

   implicit none

   ! version controlled file description for error handling, do not edit
   character(len=256), parameter :: source   = &
      "$URL$"
   character(len=32 ), parameter :: revision = "$Revision$"
   character(len=128), parameter :: revdate  = "$Date$"

   ! TODO Make dynamic?
   character(len=255) :: obs_in_file = ""
   character(len=255) :: obs_out_file = "obs_seq.aeolus"

   ! Metadata per HLOS observation
   type hlos_metadata
      real(r8) :: azimuth_a ! Azimuth angle
      integer :: obs_id ! Unique (per file) ID of observation
   end type hlos_metadata
   type(hlos_metadata), allocatable, dimension(:) :: obs_metadata

   ! Count of "WindResults" in the input file
   integer*4 :: num_mie = 0
   integer*4 :: num_rayleigh = 0

   ! For looping
   integer :: i

   ! Per observation: unique id, coordinates, HLOS, error, validity flag
   real(r8) :: latitude, longitude, altitude, azimuth_a, hlos_m, hlos_err_m
   integer :: obs_id, validity_flag, hlos, hlos_err, observation_type
   real(r8) :: obs_seconds_from_base
   type(time_type) :: obs_timestamp, prev_time
   integer :: days, seconds

   type(obs_sequence_type) :: obs_seq
   type(obs_type)          :: obs, prev_obs
   integer                 :: num_copies, num_qc, max_obs, obskey
   logical                 :: first_obs

   ! For CODA interface
   include "coda.inc"
   integer :: coda_result ! Result of call, non-zero indicates a failure
   integer*8 :: coda_pf ! TODO Why doesn't c_ptr work?
   ! type(c_ptr) :: coda_pf ! Pointer for file

   ! CODA file product class and type, for sanity checking
   character*32 product_class
   character*32 product_type


!------------
! start of executable code
!------------

   call initialize_utilities('convert_aeolus_l2b')
   call register_module(source, revision, revdate)

   ! Read input/output paths from args
   if (command_argument_count() .eq. 0) then
      write(*,*) 'Usage: convert_aeolus_l2b <input_file> (<output_file>)'
      write(*,*) 'input_file should point to an AEOLUS L2B wind profile file (.DBL)'
      write(*,*) 'If output file is omitted, then obs_seq.aeolus is used'
      stop 2
   end if
   if (command_argument_count() .ge. 1) then
      call get_command_argument(1, value=obs_in_file)
      if (command_argument_count() .ge. 2) then
         call get_command_argument(2, value=obs_out_file)
      end if
   end if

   ! time setup
   call set_calendar_type(GREGORIAN)

   ! Initialize coda library
   coda_result = coda_init()
   call handle_coda_error()

   coda_result = coda_open(obs_in_file, coda_pf)
   call handle_coda_error()

   call check_product_type_and_class()

   ! Get number of windresults
   num_mie = coda_fetch_uint32(coda_pf, [character(len=40) :: "sph", "NumMieWindResults"])
   write(*,*) 'Num of mie windresults', num_mie
   num_rayleigh = coda_fetch_uint32(coda_pf, [character(len=40) :: "sph", "NumRayleighWindResults"])
   write(*,*) 'Num of rayleigh windresults', num_rayleigh

   write(*,*) 'Total windresults' , num_mie + num_rayleigh

   ! If we found no observations abort now
   max_obs = num_mie + num_rayleigh + 1
   if (max_obs == 0) then
      write(*,*) 'Error: no windresults found'
      stop 1
   end if

   ! Start writing observations
   num_copies = 1
   num_qc = 1

   call static_init_obs_sequence()
   call init_obs(obs, num_copies, num_qc)
   call init_obs(prev_obs, num_copies, num_qc)
   call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

   call set_copy_meta_data(obs_seq, 1, 'observation')
   call set_qc_meta_data(obs_seq, 1, 'AEOLUS QC')

   ! We do two loops here, one for mie measurements and one for rayleigh measurements
   ! The code duplication is not particularly nice but I couldn't think of a way to
   ! generalize this with the data structures available in stock FORTRAN. Contributions
   ! welcome.
   do i = 0, num_mie - 1
      ! Time
      obs_seconds_from_base = coda_fetch_double(coda_pf, [character(len=40) :: "mie_geolocation", "?", "windresult_geolocation", "datetime_cog"], i)
      obs_timestamp = set_date(2000, 1, 1, 0, 0, 0) ! Base time
      obs_timestamp = increment_time(obs_timestamp, INT(obs_seconds_from_base))
      call get_time(obs_timestamp, seconds, days)

      ! Geolocation
      obs_id = coda_fetch_uint32(coda_pf, [character(len=40) :: "mie_geolocation", "?", "wind_result_id"], i)
      altitude = coda_fetch_int32(coda_pf, [character(len=40) :: "mie_geolocation", "?", "windresult_geolocation", "altitude_vcog"], i)
      azimuth_a = coda_fetch_double(coda_pf, [character(len=40) :: "mie_geolocation", "?", "windresult_geolocation", "los_azimuth"], i)

      ! For some reason, lat/lon are not correctly converted by CODA
      ! We disable automatic conversion here and do it manually
      coda_result = coda_set_option_perform_conversions(0)

      latitude = coda_fetch_int32(coda_pf, [character(len=40) :: "mie_geolocation", "?", "windresult_geolocation", "latitude_cog"], i) / 1000000.0
      longitude = coda_fetch_int32(coda_pf, [character(len=40) :: "mie_geolocation", "?", "windresult_geolocation", "longitude_cog"], i) / 1000000.0

      coda_result = coda_set_option_perform_conversions(1)

      ! HLOS & QC
      hlos = coda_fetch_int32(coda_pf, [character(len=40) :: "mie_hloswind", "?", "windresult", "mie_wind_velocity"], i)
      hlos_err = coda_fetch_uint32(coda_pf, [character(len=40) :: "mie_wind_prod_conf_data", "?", "mie_wind_qc", "hlos_error_estimate"], i)
      validity_flag = coda_fetch_uint8(coda_pf, [character(len=40) :: "mie_hloswind", "?", "windresult", "validity_flag"], i)
      observation_type = coda_fetch_uint8(coda_pf, [character(len=40) :: "mie_hloswind", "?", "windresult", "observation_type"], i)

      ! Convert HLOS from cm/s to m/s
      hlos_m = real(hlos, r8) / 100.0
      hlos_err_m = real(hlos_err, r8) / 100.0

      ! Reject flagged observations (flag == 0) and Mie-clear retrievals (type != 1)
      if (validity_flag == 0 .or. observation_type /= 1) then
         cycle
      end if

      ! Reject observations w/ instrument error larger than 12m/s
      if (hlos_err_m > 4.0) then
         cycle
      end if

      call set_aeolus_metadata(obskey, obs_id, azimuth_a)

      ! We use a hardcoded 2.5m/s value for the observation error here
      hlos_err_m = 2.5
      ! TODO Dynamically compute perhaps?
      call create_3d_obs(latitude, longitude, altitude, VERTISHEIGHT, hlos_m, AEOLUS_L2B_HLOS, hlos_err_m, days, seconds, real(validity_flag, r8), obs, obskey)

      call add_obs_to_seq(obs_seq, obs, obs_timestamp, prev_obs, prev_time, first_obs)
   end do
   do i = 0, num_rayleigh - 1
      ! Time
      obs_seconds_from_base = coda_fetch_double(coda_pf, [character(len=40) :: "rayleigh_geolocation", "?", "windresult_geolocation", "datetime_cog"], i)
      obs_timestamp = set_date(2000, 1, 1, 0, 0, 0) ! Base time
      obs_timestamp = increment_time(obs_timestamp, INT(obs_seconds_from_base))
      call get_time(obs_timestamp, seconds, days)

      ! Geolocation
      obs_id = coda_fetch_uint32(coda_pf, [character(len=40) :: "rayleigh_geolocation", "?", "wind_result_id"], i)
      altitude = coda_fetch_int32(coda_pf, [character(len=40) :: "rayleigh_geolocation", "?", "windresult_geolocation", "altitude_vcog"], i)
      azimuth_a = coda_fetch_double(coda_pf, [character(len=40) :: "rayleigh_geolocation", "?", "windresult_geolocation", "los_azimuth"], i)

      ! For some reason, lat/lon are not correctly converted by CODA
      ! We disable automatic conversion here and do it manually
      coda_result = coda_set_option_perform_conversions(0)

      latitude = coda_fetch_int32(coda_pf, [character(len=40) :: "rayleigh_geolocation", "?", "windresult_geolocation", "latitude_cog"], i) / 1000000.0
      longitude = coda_fetch_int32(coda_pf, [character(len=40) :: "rayleigh_geolocation", "?", "windresult_geolocation", "longitude_cog"], i) / 1000000.0

      coda_result = coda_set_option_perform_conversions(1)

      ! HLOS & QC
      hlos = coda_fetch_int32(coda_pf, [character(len=40) :: "rayleigh_hloswind", "?", "windresult", "rayleigh_wind_velocity"], i)
      hlos_err = coda_fetch_uint32(coda_pf, [character(len=40) :: "rayleigh_wind_prod_conf_data", "?", "rayleigh_wind_qc", "hlos_error_estimate"], i)
      validity_flag = coda_fetch_uint8(coda_pf, [character(len=40) :: "rayleigh_hloswind", "?", "windresult", "validity_flag"], i)
      observation_type = coda_fetch_uint8(coda_pf, [character(len=40) :: "rayleigh_hloswind", "?", "windresult", "observation_type"], i)

      ! Convert HLOS from cm/s to m/s
      hlos_m = real(hlos, r8) / 100.0
      hlos_err_m = real(hlos_err, r8) / 100.0

      ! Reject flagged observations (flag == 0) and Rayleigh-cloudy retrievals (type != 2)
      if (validity_flag == 0 .or. observation_type /= 2) then
         cycle
      end if

      ! Reject observations w/ instrument error larger than 12m/s
      if (hlos_err_m > 4.0) then
         cycle
      end if

      call set_aeolus_metadata(obskey, obs_id, azimuth_a)

      ! We use a hardcoded 4.0m/s value for the observation error here
      hlos_err_m = 4.0
      ! TODO Dynamically compute perhaps?
      call create_3d_obs(latitude, longitude, altitude, VERTISHEIGHT, hlos_m, AEOLUS_L2B_HLOS, hlos_err_m, days, seconds, real(validity_flag, r8), obs, obskey)

      call add_obs_to_seq(obs_seq, obs, obs_timestamp, prev_obs, prev_time, first_obs)
   end do

   call write_obs_seq(obs_seq, obs_out_file)

   ! Clean-up coda
   coda_result = coda_close(coda_pf)
   call handle_coda_error()
   call coda_done()

   call finalize_utilities()

contains
   ! Checks coda status to see if the last call succeded, otherwise prints
   ! the error message.
   ! TODO Use DART's error handler
   subroutine handle_coda_error
      implicit none
      include "coda.inc"

      integer err
      character*75 errstr

      if (coda_result .eq. 0) return

      err = coda_get_errno()
      call coda_errno_to_string(err, errstr)
      write(*,*) 'Error: ' // errstr
      stop 1
   end subroutine

   ! Verifies that the opened DBL file has the correct product type and class
   subroutine check_product_type_and_class
      implicit none
      include "coda.inc"

      ! Check if product class is AEOLUS
      coda_result = coda_get_product_class(coda_pf, product_class)
      call handle_coda_error()
      write(*,*) 'Product class = ' // product_class

      if (product_class .ne. 'AEOLUS') then
         write(*,*) 'Error: product class is not AEOLUS'
         stop 1
      end if

      ! Check if product type is ALD_U_N_2B
      coda_result = coda_get_product_type(coda_pf, product_type)
      call handle_coda_error()
      write(*,*) 'Product type = ' // product_type

      if (product_type .ne. 'ALD_U_N_2B') then
         write(*,*) 'Error: product type is not ALD_U_N_2B'
         stop 1
      end if
   end subroutine

   ! Points a cursor at the given path
   ! At any point of the path, you can use ? to indicate an array index, which is supplied using the `index` argument
   subroutine coda_point_cursor(cursor, path, index)
      integer*8, intent(in) :: cursor
      character(len=*), dimension(:), intent(in) :: path
      integer, intent(in), optional :: index

      include "coda.inc"

      integer(kind=c_long) :: index_c
      integer :: idx

      do idx = 1, size(path)
         ! Use index if path node is ?
         if (present(index) .and. trim(path(idx)) == '?') then
            index_c = index
            coda_result = coda_cursor_goto_array_element_by_index(cursor, index_c)
         else
            coda_result = coda_cursor_goto_record_field_by_name(cursor, trim(path(idx)))
         end if
         call handle_coda_error()
      end do
   end subroutine

   ! Fetch a uint32 using a coda cursor
   function coda_fetch_uint32(product, path, index) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path
      integer, intent(in), optional :: index

      include "coda.inc"

      integer*8 :: cursor
      integer*4 :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path, index)

      ! Read value
      coda_result = coda_cursor_read_uint32(cursor, r)
      call handle_coda_error()

      ! Clean-up
      call coda_cursor_delete(cursor)
      call handle_coda_error()
   end function

   ! Fetch a uint8 using a coda cursor
   function coda_fetch_uint8(product, path, index) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path
      integer, intent(in), optional :: index

      include "coda.inc"

      integer*8 :: cursor
      integer*1 :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path, index)

      ! Read value
      coda_result = coda_cursor_read_uint8(cursor, r)
      call handle_coda_error()

      ! Clean-up
      call coda_cursor_delete(cursor)
      call handle_coda_error()
   end function

   ! Fetch a int32 using a coda cursor
   function coda_fetch_int32(product, path, index) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path
      integer, intent(in), optional :: index

      include "coda.inc"

      integer*8 :: cursor
      integer*4 :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path, index)

      ! Read value
      coda_result = coda_cursor_read_int32(cursor, r)
      call handle_coda_error()

      ! Clean-up
      call coda_cursor_delete(cursor)
      call handle_coda_error()
   end function

   ! Fetch a double using a coda cursor
   function coda_fetch_double(product, path, index) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path
      integer, intent(in), optional :: index

      include "coda.inc"

      integer*8 :: cursor
      real(r8) :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path, index)

      ! Read value
      coda_result = coda_cursor_read_double(cursor, r)
      call handle_coda_error()

      ! Clean-up
      call coda_cursor_delete(cursor)
      call handle_coda_error()
   end function

   ! Fetch a float using a coda cursor
   function coda_fetch_float(product, path, index) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path
      integer, intent(in), optional :: index

      include "coda.inc"

      integer*8 :: cursor
      real(r8) :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path, index)

      ! Read value
      coda_result = coda_cursor_read_float(cursor, r)
      call handle_coda_error()

      ! Clean-up
      call coda_cursor_delete(cursor)
      call handle_coda_error()
   end function

end program

! ! <next few lines under version control, do not edit>
! ! $URL$
! ! $Id$
! ! $Revision$
! ! $Date$
