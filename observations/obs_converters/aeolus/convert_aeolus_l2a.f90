! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_aeolus_l2a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_aeolus_l2a - program that reads an AEOLUS L2A aerosol profile
!                        file and writes a DART obs_seq file using the
!                        DsART library routines.
!
!     created April. 2024 Thanasis Georgiou, National Observatory of Athens
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
   use      obs_kind_mod, only : LIDAR_EXTINCTION
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

   ! I/O paths
   character(len=255) :: obs_in_file = ""
   character(len=255) :: obs_out_file = "obs_seq.aeolus"

   ! Count of profiles in the input file
   integer*4 :: num_mle_profiles = 0

   ! For looping
   integer :: i, meas, height

   ! Path for indexing into coda arrays
   character(len=40), allocatable :: path(:)

   ! Observation geolocation
   real(r8), dimension(30, 25) :: latitude, longitude, altitude
   real(r8), dimension(30, 24) :: latitude_mid, longitude_mid, altitude_mid
   real(r8), dimension(24) :: latitude_mean, longitude_mean, altitude_mean
   real(r8):: centroid_time

   ! Feature mask
   integer, dimension(30, 24) :: feature_mask
   integer, dimension(24) :: feature_mask_sum
   integer fm

   ! Observation data
   real(r8), dimension(24) :: extinction, extinction_stdev
   real(r8) :: extinction_err, max_extinction = -1000.0

   ! Observation time
   real(r8) :: obs_seconds_from_base
   type(time_type) :: obs_timestamp, prev_time
   integer :: days, seconds

   ! Counts for rejected observations
   integer :: rejected_feature_mask = 0, rejected_no_retrieval = 0, rejected_below_ground = 0

   ! Obs. sequence and related metadata
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

   call initialize_utilities('convert_aeolus_l2a')
   call register_module(source, revision, revdate)

   ! Read input/output paths from args
   if (command_argument_count() .eq. 0) then
      write(*,*) 'Usage: convert_aeolus_l2a <input_file> (<output_file>)'
      write(*,*) 'input_file should point to an AEOLUS L2A aerosol profile file (.DBL)'
      write(*,*) 'If output file is omitted, then obs_seq.aeolus is used'
      stop
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
   num_mle_profiles = coda_fetch_int32(coda_pf, [character(len=12) :: "sph", "num_prof_mle"])
   max_obs = num_mle_profiles * 24
   write(*,*) 'Num of MLE profiles', num_mle_profiles
   write(*,*) 'Num of MLE bins', max_obs

   ! If we found no observations abort now
   if (max_obs == 0) then
      write(*,*) 'Error: no data found'
      stop
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

   ! The outer loop targets each MLE profile, while the inner loops iterates over the vertial bins
   do i = 0, num_mle_profiles - 1
      ! Read geolocation array for all measurements
      ! The path is /geolocation[i]/measurement_geolocation[meas]/rayleigh_geolocation_height_bin[height]/VARIABLE
      if (allocated(path)) deallocate(path)
      allocate(path(7))
      path(1) = "geolocation"
      write(path(2), '(I0)') i
      path(2) = "$" // trim(path(2))
      path(3) = "measurement_geolocation"
      path(5) = "rayleigh_geolocation_height_bin"

      do meas = 1, 30
         write(path(4), '(I0)') meas - 1
         path(4) = "$" // trim(path(4))

         do height = 1, 25
            write(path(6), '(I0)') height - 1
            path(6) = "$" // trim(path(6))

            coda_result = coda_set_option_perform_conversions(0)

            path(7) = "latitude_of_height_bin"
            latitude(meas, height) = coda_fetch_int32(coda_pf, path) / 1000000.0

            path(7) = "longitude_of_height_bin"
            longitude(meas, height) = coda_fetch_int32(coda_pf, path) / 1000000.0

            coda_result = coda_set_option_perform_conversions(1)

            path(7) = "altitude_of_height_bin"
            altitude(meas, height) = coda_fetch_double(coda_pf, path)
         end do
      end do

      ! Read centroid time
      deallocate(path)
      allocate(path(3))

      path(1) = "sca_mle_opt_properties"
      write(path(2), '(I0)') i
      path(2) = "$" // trim(path(2))
      path(3) = "starttime"
      centroid_time = coda_fetch_double(coda_pf, path)

      ! Read feature mask
      deallocate(path)
      allocate(path(5))
      path(1) = "feature_mask"
      write(path(2), '(I0)') i
      path(2) = "$" // trim(path(2))
      path(3) = "feature_mask_indices"
      path(5) = "feature_mask_index"
      do meas = 1, 30
         ! write(path(4), '(I0)') meas - 1
         ! path(4) = "$" // trim(path(4))

         do height = 1, 24
            write(path(4), '(I0)') (meas - 1) * 24 + (height - 1)
            path(4) = "$" // trim(path(4))

            fm = coda_fetch_int8(coda_pf, path)
            ! Feature masks are:
            ! - Negative values are surface, attenuated, ...
            ! - 0 is clear sky
            ! - 0-4 are mostly molecular
            ! - 10 is certainy cloud
            ! We reject anything outside the [0, 8]
            if (fm < 0 .or. fm >= 6) then
               feature_mask(meas, height) = 0
            else
               feature_mask(meas, height) = 1
            end if
         end do
      end do


      ! Unstagger the arrays
      latitude_mid = (latitude(:, 1:24) + latitude(:, 2:25)) / 2.0
      longitude_mid = (longitude(:, 1:24) + longitude(:, 2:25)) / 2.0
      altitude_mid = (altitude(:, 1:24) + altitude(:, 2:25)) / 2.0

      ! Compute average of all measurements
      latitude_mean = sum(latitude_mid, dim=1) / 30.0
      longitude_mean = sum(longitude_mid, dim=1) / 30.0
      altitude_mean = sum(altitude_mid, dim=1) / 30.0

      ! Sum the feature mask
      feature_mask_sum = sum(feature_mask, dim=1)

      ! Read extinction
      extinction(:) = missing_r8

      deallocate(path)
      allocate(path(5)) ! Reuse path
      path(1) = "sca_mle_opt_properties"
      write(path(2), '(I0)') i
      path(2) = "$" // trim(path(2))
      path(3) = "sca_mle_optical_properties"
      path(5) = "extinction"
      do height = 1, 24
         write(path(4), '(I0)') height - 1
         path(4) = "$" // trim(path(4))

         ! Documentation says that Extinction is given in [1e-6/m] and that the missing value
         ! is -1e6, but in practice it seems the missing value is -1. I still take it that
         ! the value is in [1e-6/m] and convert it to [1/m]
         extinction(height) = coda_fetch_double(coda_pf, path)
         if (extinction(height) == -1.0) then
            extinction(height) = missing_r8
         else
            extinction(height) = extinction(height) / 1e6
         end if
      end do

      ! Read extinction variance
      path(1) = "sca_mle_pcd"
      path(3) = "profile_pcd_bins"
      path(5) = "extinction_variance"
      do height = 1, 24
         write(path(4), '(I0)') height - 1
         path(4) = "$" // trim(path(4))

         ! Variance is in [1/m^2] units, square root of it is the standard deviation [1/m]
         extinction_stdev(height) = sqrt(coda_fetch_double(coda_pf, path))
      end do

      ! Parse time
      obs_timestamp = set_date(2000, 1, 1, 0, 0, 0) ! Base time
      obs_timestamp = increment_time(obs_timestamp, INT(centroid_time))
      call get_time(obs_timestamp, seconds, days)

      if (i .eq. 0) then
         write(*,*) "First timestamp:"
         call print_date(obs_timestamp)
      end if


      ! Write observations
      do height = 1, 24
         ! Skip missing values
         if (extinction(height) == missing_r8) then
            rejected_no_retrieval = rejected_no_retrieval + 1
            cycle
         end if
         if (altitude_mean(height) < 0.0) then
            rejected_below_ground = rejected_below_ground + 1
            cycle
         end if

         ! If more than half of the measurement feature masks in this BRC are rejected, skip
         ! Feature mask of 1 means we accept, so less than 15 accepted bins means we skip
         if (feature_mask_sum(height) .lt. 15) then
            rejected_feature_mask = rejected_feature_mask + 1
            cycle
         end if

         ! Set error to 25% of the extinction value, unless extinction is 0
         extinction_err = extinction(height) * 0.25
         if (extinction_err .eq. 0.0) extinction_err = 0.05

         ! write(*,*) 'Extinction', extinction(height), 'STDEV ', extinction_stdev(height), 'STDEV (%)', (extinction_stdev(height) / extinction(height))*100
         ! write(*,*) 'Extinction', extinction(height), 'Error ', extinction_err, 'Error (%)', (extinction_err / extinction(height))*100
         if (extinction(height) > max_extinction) max_extinction = extinction(height)

         ! Set observation error
         call create_3d_obs(latitude_mean(height), longitude_mean(height), altitude_mean(height), VERTISHEIGHT, extinction(height), LIDAR_EXTINCTION, extinction_stdev(height), days, seconds, real(1.0, r8), obs)

         call add_obs_to_seq(obs_seq, obs, obs_timestamp, prev_obs, prev_time, first_obs)
      end do
   end do

   call write_obs_seq(obs_seq, obs_out_file)

   ! Write some statistics
   write(*,*) 'Total obs in file:', num_mle_profiles * 24
   write(*,*) 'Rejected due to no retrieval:', rejected_no_retrieval
   write(*,*) 'Rejected due to below ground:', rejected_below_ground
   write(*,*) 'Rejected due to feature mask:', rejected_feature_mask
   write(*,*) 'Total obs written:', get_num_obs(obs_seq)
   write(*,*) 'Max ext:', max_extinction

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
      stop
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
         stop
      end if

      ! Check if product type is ALD_U_N_2B
      coda_result = coda_get_product_type(coda_pf, product_type)
      call handle_coda_error()
      write(*,*) 'Product type = ' // product_type

      if (product_type .ne. 'ALD_U_N_2A') then
         write(*,*) 'Error: product type is not ALD_U_N_2A'
         stop
      end if
   end subroutine

   ! Points a cursor at the given path
   ! At any point of the path, you can use a integer prefixed by a dollar sign (e.g. $12) to indicate and array index
   subroutine coda_point_cursor(cursor, path)
      integer*8, intent(in) :: cursor
      character(len=*), dimension(:), intent(in) :: path

      include "coda.inc"

      integer(kind=c_long) :: index_c
      character(len=len(path(1))) :: node
      integer :: idx, node_len
      integer :: coda_idx

      do idx = 1, size(path)
         ! Use index if path node starts with $
         node = path(idx)

         if (node(1:1) == '$') then
            read(node(2:len(node)), *) coda_idx
            index_c = int(coda_idx, kind=c_long)
            coda_result = coda_cursor_goto_array_element_by_index(cursor, index_c)
         else
            coda_result = coda_cursor_goto_record_field_by_name(cursor, trim(path(idx)))
         end if
         call handle_coda_error()
      end do
   end subroutine

   ! Fetch a uint32 using a coda cursor
   function coda_fetch_uint32(product, path) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path

      include "coda.inc"

      integer*8 :: cursor
      integer*4 :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path)

      ! Read value
      coda_result = coda_cursor_read_uint32(cursor, r)
      call handle_coda_error()

      ! Clean-up
      call coda_cursor_delete(cursor)
      call handle_coda_error()
   end function

   ! Fetch a uint8 using a coda cursor
   function coda_fetch_uint8(product, path) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path

      include "coda.inc"

      integer*8 :: cursor
      integer*1 :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path)

      ! Read value
      coda_result = coda_cursor_read_uint8(cursor, r)
      call handle_coda_error()

      ! Clean-up
      call coda_cursor_delete(cursor)
      call handle_coda_error()
   end function

   ! Fetch a int8 using a coda cursor
   function coda_fetch_int8(product, path) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path

      include "coda.inc"

      integer*8 :: cursor
      integer*1 :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path)

      ! Read value
      coda_result = coda_cursor_read_int8(cursor, r)
      call handle_coda_error()

      ! Clean-up
      call coda_cursor_delete(cursor)
      call handle_coda_error()
   end function


   ! Fetch a int32 using a coda cursor
   function coda_fetch_int32(product, path) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path

      include "coda.inc"

      integer*8 :: cursor
      integer*4 :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path)

      ! Read value
      coda_result = coda_cursor_read_int32(cursor, r)
      call handle_coda_error()

      ! Clean-up
      call coda_cursor_delete(cursor)
      call handle_coda_error()
   end function

   ! Fetch a double using a coda cursor
   function coda_fetch_double(product, path) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path

      include "coda.inc"

      integer*8 :: cursor
      real(r8) :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path)

      ! Read value
      coda_result = coda_cursor_read_double(cursor, r)
      call handle_coda_error()

      ! Clean-up
      call coda_cursor_delete(cursor)
      call handle_coda_error()
   end function

   ! Fetch a float using a coda cursor
   function coda_fetch_float(product, path) result(r)
      integer*8, intent(in) :: product
      character(len=*), dimension(:), intent(in) :: path

      include "coda.inc"

      integer*8 :: cursor
      real(r8) :: r

      ! Prepare a cursor for reading data
      cursor = coda_cursor_new()
      call handle_coda_error()
      coda_result = coda_cursor_set_product(cursor, coda_pf)
      call handle_coda_error()

      call coda_point_cursor(cursor, path)

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
