! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! AEOLUS_L2B_HLOS,  QTY_HLOS_WIND
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_aeolus_hlos_mod, only : set_aeolus_metadata, &
!                                           read_aeolus_metadata, &
!                                           write_aeolus_metadata, &
!                                           interactive_aeolus_metadata, &
!                                           get_hlos
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(AEOLUS_L2B_HLOS)
!         call get_hlos(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(AEOLUS_L2B_HLOS)
!      call read_aeolus_metadata(obs_def%key, key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(AEOLUS_L2B_HLOS)
!      call write_aeolus_metadata(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(AEOLUS_L2B_HLOS)
!      call interactive_aeolus_metadata(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_aeolus_hlos_mod

   use        types_mod, only : r8, PI, metadatalength, MISSING_R8, DEG2RAD
   use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
      logfileunit, get_unit, open_file, close_file, &
      file_exist, ascii_file_format
   use     location_mod, only : location_type, set_location, get_location, &
      VERTISHEIGHT, VERTISLEVEL, set_location_missing
   use     obs_kind_mod, only : QTY_HLOS_WIND, QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT
   use  assim_model_mod, only : interpolate

   use obs_def_utilities_mod, only : track_status
   use ensemble_manager_mod,  only : ensemble_type

   implicit none
   private

   public :: set_aeolus_metadata, read_aeolus_metadata, write_aeolus_metadata, &
      interactive_aeolus_metadata, get_hlos

! version controlled file description for error handling, do not edit
   character(len=256), parameter :: source   = "$URL$"
   character(len=32 ), parameter :: revision = "$Revision$"
   character(len=128), parameter :: revdate  = "$Date$"

   character(len=256) :: string1, string2
   logical, save      :: module_initialized = .false.

! Metadata for each HLOS observation
! The azimuth angle is required to compute the HLOS wind vector from the U
! and V wind components
! These slides have some info about the geometry of Aeolus' measurements
! https://confluence.ecmwf.int/download/attachments/46596815/aeolus_obs_operator.pdf?version=1&modificationDate=1524579985341&api=v2

   type hlos_metadata
      private
      real(r8) :: azimuth_a ! Azimuth angle
      integer :: obs_id ! Unique (per file) ID of observation
   end type hlos_metadata

   type(hlos_metadata), allocatable, dimension(:) :: obs_metadata
   type(hlos_metadata) :: missing_metadata
   character(len=6), parameter :: HLOSSTRING = 'hlos'

   logical :: debug = .FALSE.
   integer :: maxkey = 30*50000  ! 30 vertical bins, 50000 profiles
   integer :: lastKey = 0

contains

! Allocates the array for storing the observation metadata
   subroutine initialize_module
      if (module_initialized) return
      call register_module(source, revision, revdate)

      missing_metadata%azimuth_a = MISSING_R8

      allocate(obs_metadata(maxkey))
      obs_metadata(:) = missing_metadata

      module_initialized = .true.
   end subroutine initialize_module

   subroutine set_aeolus_metadata(key, obs_id, azimuth_a)
      integer, intent(out) :: key
      integer, intent(in) :: obs_id
      real(r8), intent(in) :: azimuth_a

      if (.not. module_initialized) call initialize_module

      ! TODO check if key is in range/grow array maybe?

      lastKey = lastKey + 1
      key = lastKey

      obs_metadata(key)%azimuth_a = azimuth_a
      obs_metadata(key)%obs_id = obs_id
   end subroutine set_aeolus_metadata

   ! Reads metadata for a single observation and stores into `obs_metadata`
   subroutine read_aeolus_metadata(key, obsId, ifile, fform)
      integer, intent(out) :: key    ! index into local metadata
      integer, intent(in) :: obsID
      integer, intent(in) :: ifile
      character(len=*), intent(in), optional :: fform

      logical :: is_ascii
      integer :: ierr
      character(len=6) :: header

      if (.not. module_initialized) call initialize_module

      lastKey = lastKey + 1
      key = lastKey

      is_ascii = ascii_file_format(fform)
      write(string2, *)'observation #', obsID

      if (is_ascii) then
         read(ifile, *, iostat=ierr) header
         call check_iostat(ierr, 'read_aeolus_metadata', 'header', string2)
         if (trim(header) /= trim(HLOSSTRING)) then
            write(string1, *) 'Expected HLOS header [', trim(HLOSSTRING), '] but got [', trim(header), ']'
            call error_handler(E_ERR, 'read_aeolus_metadata', string1, source, &
               revision, revdate, text2=string2)
         endif

         read(ifile, *, iostat=ierr) obs_metadata(key)%azimuth_a
         read(ifile, *, iostat=ierr) obs_metadata(key)%obs_id
         call check_iostat(ierr, 'read_aeolus_metadata', 'azimuth_a', string2)
      else
         read(ifile, iostat=ierr) header
         call check_iostat(ierr, 'read_aeolus_metadata', 'header', string2)
         if (trim(header) /= trim(HLOSSTRING)) then
            write(string1, *) 'Expected HLOS header [', trim(HLOSSTRING), '] but got [', trim(header), ']'
            call error_handler(E_ERR, 'read_aeolus_metadata', string1, source, &
               revision, revdate, text2=string2)
         endif

         read(ifile, iostat=ierr) obs_metadata(key)%azimuth_a
         read(ifile, *, iostat=ierr) obs_metadata(key)%obs_id
         call check_iostat(ierr, 'read_aeolus_metadata', 'azimuth_a', string2)
      endif
   end subroutine read_aeolus_metadata

   ! Writes metadata for one observation into the file
   subroutine write_aeolus_metadata(key, ifile, fform)
      integer, intent(in) :: key
      integer, intent(in) :: ifile
      character(len=*), intent(in), optional :: fform

      logical :: is_ascii
      integer :: ierr

      if (.not. module_initialized) call initialize_module

      is_ascii = ascii_file_format(fform)

      if (is_ascii) then
         write(ifile, *) HLOSSTRING
         write(ifile, *) obs_metadata(key)%azimuth_a
         write(ifile, *) obs_metadata(key)%obs_id
      else
         write(ifile) HLOSSTRING
         write(ifile) obs_metadata(key)%azimuth_a
         write(ifile) obs_metadata(key)%obs_id
      endif
   end subroutine write_aeolus_metadata

   ! For interactively creating an HLOS measurement
   ! TODO maybe add some error checking
   subroutine interactive_aeolus_metadata(key)
      integer, intent(in) :: key

      if (.not. module_initialized) call initialize_module

      write(string1, *) 'Azimuth angle (degrees)'
      read(*, *) obs_metadata(key)%azimuth_a
   end subroutine interactive_aeolus_metadata

   ! Gets the HLOS wind vector from the U and V wind components
   subroutine get_hlos(state_handle, ens_size, location, key, val, istatus)
      type(ensemble_type), intent(in) :: state_handle
      integer, intent(in) :: ens_size
      type(location_type), intent(in) :: location
      integer, intent(in) :: key
      real(r8), intent(out) :: val(ens_size)
      integer, intent(out) :: istatus(ens_size)

      real(r8) :: u(ens_size), v(ens_size), az_rad(ens_size)

      ! First grab U and V from the model state
      call interpolate(state_handle, ens_size, location, QTY_U_WIND_COMPONENT, u, istatus)
      call interpolate(state_handle, ens_size, location, QTY_V_WIND_COMPONENT, v, istatus)

      ! Then compute the HLOS wind vector
      az_rad = obs_metadata(key)%azimuth_a * DEG2RAD
      val = -u * sin(az_rad) - v * cos(az_rad)
   end subroutine get_hlos


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


end module obs_def_aeolus_hlos_mod
! END DART PREPROCESS MODULE CODE
