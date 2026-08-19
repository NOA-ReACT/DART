! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

!> Optical properties of the GOCART dust and sea salt bins.
!>
!> Reads the look-up table shared by the GOCART AOD (obs_def_GOCART_AOD_mod) and the
!> atmospheric lidar (obs_def_atm_lidar_mod) forward operators, and interpolates it to the
!> ambient hygroscopic growth factor. It mirrors the offline python operator
!> (`extinction-operator/extinction_operator.py`) so that the two produce the same numbers.
!>
!> The table is a CSV with one row per (species, wavelength) for humidity-independent species
!> (dust) and one row per (species, wavelength, growth factor) for hygroscopic ones (sea salt):
!>
!>   wavelength[nm],species,ext_eff[meter ** 2 / kilogram],fine_fraction[1],growth_factor[1]
!>   532,DUST_1,960.7737123287676,0.49405568193826266,
!>   532,SEAS_1,3007.590681818182,1.0,1.0
!>
!> Only the first three columns are mandatory. `growth_factor[1]` is blank for the species that
!> do not take up water, and a table without the column at all reproduces the old, purely dry
!> behaviour. `fine_fraction[1]` is the share of a bin's extinction coming from particles below
!> the ambient-size cutoff the table was built with, and is only needed for the fine/coarse mode
!> split.

module gocart_optics_mod

   use types_mod, only : r8, MISSING_R8, MISSING_I, &
      t_kelvin, es_alpha, es_beta, es_gamma
   use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
      find_namelist_in_file, check_namelist_read, NAMELIST_NOT_PRESENT, &
      do_nml_file, do_nml_term, nmlfileunit
   use read_csv_mod, only : csv_file_type, csv_open, csv_close, csv_get_nrows, &
      csv_get_field, csv_field_exists
   use sort_mod, only : index_sort
   use location_mod, only : location_type
   use obs_kind_mod, only : QTY_TEMPERATURE, QTY_PRESSURE, QTY_VAPOR_MIXING_RATIO, &
      QTY_GC_DUST_BIN1, QTY_GC_DUST_BIN2, QTY_GC_DUST_BIN3, &
      QTY_GC_DUST_BIN4, QTY_GC_DUST_BIN5, &
      QTY_GC_SEAS_BIN1, QTY_GC_SEAS_BIN2, QTY_GC_SEAS_BIN3, &
      QTY_GC_SEAS_BIN4
   use assim_model_mod, only : interpolate
   use obs_def_utilities_mod, only : track_status
   use ensemble_manager_mod, only : ensemble_type

   implicit none
   private

   public :: gocart_optics_init, &
      N_GOCART_BINS, GOCART_BIN_QTYS, GOCART_BIN_NAMES, &
      bin_optics_type, get_bin_optics, bin_ext_eff, bin_fine_fraction, &
      optics_are_hygroscopic, optics_have_fine_fraction, &
      get_expected_growth_factor

   character(len=*), parameter :: source   = 'gocart_optics_mod.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''

   ! The nine GOCART aerosol bins, in the order the tables and the operators use them
   integer, parameter :: N_GOCART_BINS = 9
   integer, parameter :: GOCART_BIN_QTYS(N_GOCART_BINS) = &
      [QTY_GC_DUST_BIN1, QTY_GC_DUST_BIN2, QTY_GC_DUST_BIN3, QTY_GC_DUST_BIN4, QTY_GC_DUST_BIN5, &
       QTY_GC_SEAS_BIN1, QTY_GC_SEAS_BIN2, QTY_GC_SEAS_BIN3, QTY_GC_SEAS_BIN4]
   character(len=6), parameter :: GOCART_BIN_NAMES(N_GOCART_BINS) = &
      [character(len=6) :: 'DUST_1', 'DUST_2', 'DUST_3', 'DUST_4', 'DUST_5', &
                           'SEAS_1', 'SEAS_2', 'SEAS_3', 'SEAS_4']

   ! Column headers of the optical properties CSV
   character(len=*), parameter :: WAVELENGTH_COLUMN    = 'wavelength[nm]'
   character(len=*), parameter :: SPECIES_COLUMN       = 'species'
   character(len=*), parameter :: EXT_EFF_COLUMN       = 'ext_eff[meter ** 2 / kilogram]'
   character(len=*), parameter :: FINE_FRACTION_COLUMN = 'fine_fraction[1]'
   character(len=*), parameter :: GROWTH_FACTOR_COLUMN = 'growth_factor[1]'

   ! Ratio of the molar masses of water vapour and dry air. Spelled out rather than taken from
   ! `gas_constant / gas_constant_v` so that it matches the python operator digit for digit.
   real(r8), parameter :: EPSILON_WATER = 0.622_r8

   !> The optical properties of one aerosol bin at one wavelength. Humidity-independent species
   !> have a single node; hygroscopic ones carry the table's growth factor axis, in ascending
   !> order. `fine_fraction` is only allocated when the table has the column.
   type :: bin_optics_type
      integer :: wavelength = MISSING_I  ! [nm]
      integer :: qty        = MISSING_I  ! DART quantity of the bin
      integer :: nnodes     = 0          ! 1 when the bin does not take up water
      real(r8), allocatable :: growth_factor(:)  ! [1], ascending
      real(r8), allocatable :: ext_eff(:)        ! [meter ** 2 / kilogram], per dry mass
      real(r8), allocatable :: fine_fraction(:)  ! [1]
   end type bin_optics_type

   type(bin_optics_type), allocatable :: optics(:)
   integer :: n_optics = 0

   logical :: module_initialized = .false.
   logical :: have_growth_factor = .false.
   logical :: have_fine_fraction = .false.
   logical :: warned_about_humidity = .false.

   character(len=512) :: string1, string2

   ! Namelist
   character(len=256) :: optical_properties_file = 'GOCART_AOD_optical_properties.csv'
   ! Hygroscopicity parameter (kappa-Koehler) of sea salt. Zieger et al. (2017,
   ! doi:10.1038/ncomms15883): sea salt grows less than pure NaCl (kappa 1.28) due to hydrates.
   real(r8) :: kappa  = 1.1_r8
   ! The growth factor diverges as RH -> 1, so RH is clamped here. This is also the ceiling of
   ! the growth factor axis of the shipped tables.
   real(r8) :: max_rh = 0.99_r8
   logical  :: debug  = .false.

   namelist /gocart_optics_nml/ optical_properties_file, kappa, max_rh, debug

contains

!----------------------------------------------------------------------
!> Reads the namelist and the optical properties table. Safe to call repeatedly.

subroutine gocart_optics_init()

integer :: iunit, io

if (module_initialized) return
module_initialized = .true.

call register_module(source, revision, revdate)

call find_namelist_in_file('input.nml', 'gocart_optics_nml', iunit, optional_nml = .true.)
if (iunit /= NAMELIST_NOT_PRESENT) then
   read(iunit, nml = gocart_optics_nml, iostat = io)
   call check_namelist_read(iunit, io, 'gocart_optics_nml')
endif
if (do_nml_file()) write(nmlfileunit, nml = gocart_optics_nml)
if (do_nml_term()) write(     *     , nml = gocart_optics_nml)

if (kappa < 0.0_r8) then
   write(string1, *) 'kappa must not be negative, got ', kappa
   call error_handler(E_ERR, 'gocart_optics_init', string1, source, revision, revdate)
endif
if (max_rh <= 0.0_r8 .or. max_rh >= 1.0_r8) then
   write(string1, *) 'max_rh must be in (0, 1), got ', max_rh
   call error_handler(E_ERR, 'gocart_optics_init', string1, source, revision, revdate)
endif

call read_optical_properties()

end subroutine gocart_optics_init

!----------------------------------------------------------------------
!> Reads the CSV into one `bin_optics_type` per (wavelength, species) present in it.

subroutine read_optical_properties()

character(len=*), parameter :: routine = 'read_optical_properties'

type(csv_file_type) :: cf
integer  :: nrows, nwaves, nentries, irow, iwave, ibin, ientry, nnodes
integer,  allocatable :: row_wavelength(:), row_qty(:), wavelengths(:), rows(:), order(:)
real(r8), allocatable :: row_wavelength_real(:), row_ext_eff(:), row_fine_fraction(:), &
                         row_growth_factor(:)
character(len=64), allocatable :: row_species(:)

call csv_open(optical_properties_file, cf, forced_delim = ',', context = routine)

nrows = csv_get_nrows(cf)
if (nrows <= 0) then
   write(string1, *) 'No optical properties in ' // trim(optical_properties_file)
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
endif

allocate(row_wavelength_real(nrows), row_species(nrows), row_ext_eff(nrows), &
         row_fine_fraction(nrows), row_growth_factor(nrows), &
         row_wavelength(nrows), row_qty(nrows))

! The wavelength is read as a real and rounded: the tables have been written both as '532'
! and as '532.0', and string_to_integer() rejects the latter.
call csv_get_field(cf, WAVELENGTH_COLUMN, row_wavelength_real, context = routine)
call csv_get_field(cf, SPECIES_COLUMN,    row_species,         context = routine)
call csv_get_field(cf, EXT_EFF_COLUMN,    row_ext_eff,         context = routine)

! Blank fields come back as MISSING_R8, which is how the dust rows of a growth-factor-resolved
! table say that they do not take up water.
have_fine_fraction = csv_field_exists(cf, FINE_FRACTION_COLUMN)
if (have_fine_fraction) then
   call csv_get_field(cf, FINE_FRACTION_COLUMN, row_fine_fraction, context = routine)
else
   row_fine_fraction(:) = MISSING_R8
endif

row_growth_factor(:) = MISSING_R8
if (csv_field_exists(cf, GROWTH_FACTOR_COLUMN)) then
   call csv_get_field(cf, GROWTH_FACTOR_COLUMN, row_growth_factor, context = routine)
endif
have_growth_factor = any(row_growth_factor /= MISSING_R8)

call csv_close(cf)

do irow = 1, nrows
   row_wavelength(irow) = nint(row_wavelength_real(irow))
   row_qty(irow) = species_to_qty(row_species(irow))
   if (row_ext_eff(irow) == MISSING_R8 .and. row_qty(irow) /= MISSING_I) then
      write(string1, *) 'Missing ' // EXT_EFF_COLUMN // ' on row ', irow, &
                        ' of ' // trim(optical_properties_file)
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   endif
end do

call unique_wavelengths(row_wavelength, row_qty, wavelengths, nwaves)

! One entry per (wavelength, bin) that has at least one row
allocate(optics(nwaves * N_GOCART_BINS))
allocate(rows(nrows), order(nrows))
nentries = 0

do iwave = 1, nwaves
   do ibin = 1, N_GOCART_BINS

      nnodes = 0
      do irow = 1, nrows
         if (row_wavelength(irow) /= wavelengths(iwave)) cycle
         if (row_qty(irow) /= GOCART_BIN_QTYS(ibin)) cycle
         nnodes = nnodes + 1
         rows(nnodes) = irow
      end do
      if (nnodes == 0) cycle

      nentries = nentries + 1
      ientry = nentries
      optics(ientry)%wavelength = wavelengths(iwave)
      optics(ientry)%qty        = GOCART_BIN_QTYS(ibin)
      optics(ientry)%nnodes     = nnodes

      allocate(optics(ientry)%growth_factor(nnodes))
      allocate(optics(ientry)%ext_eff(nnodes))
      if (have_fine_fraction) allocate(optics(ientry)%fine_fraction(nnodes))

      if (nnodes == 1) then
         ! Humidity-independent species: a single node, pinned at the dry radius
         optics(ientry)%growth_factor(1) = 1.0_r8
      else
         ! Several rows for one (species, wavelength) only make sense if they sit on a growth
         ! factor axis we can interpolate along
         do irow = 1, nnodes
            if (row_growth_factor(rows(irow)) /= MISSING_R8) cycle
            write(string1, *) trim(GOCART_BIN_NAMES(ibin)) // ' has ', nnodes, ' rows at ', &
                              wavelengths(iwave), ' nm but no ' // GROWTH_FACTOR_COLUMN
            call error_handler(E_ERR, routine, string1, source, revision, revdate)
         end do
         call index_sort(row_growth_factor(rows(1:nnodes)), order(1:nnodes), nnodes)
         rows(1:nnodes) = rows(order(1:nnodes))
         optics(ientry)%growth_factor(:) = row_growth_factor(rows(1:nnodes))
      endif

      optics(ientry)%ext_eff(:) = row_ext_eff(rows(1:nnodes))
      if (have_fine_fraction) then
         do irow = 1, nnodes
            if (row_fine_fraction(rows(irow)) /= MISSING_R8) cycle
            write(string1, *) 'Missing ' // FINE_FRACTION_COLUMN // ' for ' // &
                              trim(GOCART_BIN_NAMES(ibin)) // ' at ', wavelengths(iwave), ' nm'
            call error_handler(E_ERR, routine, string1, source, revision, revdate)
         end do
         optics(ientry)%fine_fraction(:) = row_fine_fraction(rows(1:nnodes))
      endif

   end do
end do

if (nentries == 0) then
   write(string1, *) 'No known GOCART species found in ' // trim(optical_properties_file)
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
endif
n_optics = nentries

write(string1, '(A,I0,A,I0,A)') 'Read ', nentries, ' optical property entries over ', &
                                nwaves, ' wavelengths'
write(string2, '(A,L1,A,L1)') 'from ' // trim(optical_properties_file) // &
                              ': hygroscopic = ', have_growth_factor, &
                              ', fine/coarse split available = ', have_fine_fraction
call error_handler(E_MSG, routine, string1, source, revision, revdate, text2=string2)

if (debug) then
   do ientry = 1, nentries
      write(*, *) 'wavelength ', optics(ientry)%wavelength, ' nm, quantity ', &
                  optics(ientry)%qty, ', ', optics(ientry)%nnodes, ' growth factor node(s), ', &
                  'ext_eff ', optics(ientry)%ext_eff(1), ' .. ', &
                  optics(ientry)%ext_eff(optics(ientry)%nnodes)
   end do
endif

deallocate(row_wavelength_real, row_species, row_ext_eff, row_fine_fraction, &
           row_growth_factor, row_wavelength, row_qty, wavelengths, rows, order)

end subroutine read_optical_properties

!----------------------------------------------------------------------
!> The DART quantity of a GOCART species name, MISSING_I if we do not know the species.

function species_to_qty(species) result(qty)

character(len=*), intent(in) :: species
integer                      :: qty

integer :: ibin

qty = MISSING_I
do ibin = 1, N_GOCART_BINS
   if (trim(species) /= trim(GOCART_BIN_NAMES(ibin))) cycle
   qty = GOCART_BIN_QTYS(ibin)
   return
end do

end function species_to_qty

!----------------------------------------------------------------------
!> The distinct wavelengths of the rows belonging to a known species, in the order they first
!> appear in the file.

subroutine unique_wavelengths(row_wavelength, row_qty, wavelengths, nwaves)

integer,              intent(in)  :: row_wavelength(:)
integer,              intent(in)  :: row_qty(:)
integer, allocatable, intent(out) :: wavelengths(:)
integer,              intent(out) :: nwaves

integer :: irow, iwave

allocate(wavelengths(size(row_wavelength)))
nwaves = 0

ROWS: do irow = 1, size(row_wavelength)
   if (row_qty(irow) == MISSING_I) cycle ROWS
   do iwave = 1, nwaves
      if (wavelengths(iwave) == row_wavelength(irow)) cycle ROWS
   end do
   nwaves = nwaves + 1
   wavelengths(nwaves) = row_wavelength(irow)
end do ROWS

end subroutine unique_wavelengths

!----------------------------------------------------------------------
!> Whether the table resolves the optical properties over a growth factor axis, i.e. whether
!> the operators have to work out the ambient humidity. A table without the column reproduces
!> the old dry behaviour exactly.

function optics_are_hygroscopic() result(hygroscopic)

logical :: hygroscopic

call gocart_optics_init()
hygroscopic = have_growth_factor

end function optics_are_hygroscopic

!----------------------------------------------------------------------
!> Whether the table carries the fine fraction needed for the fine/coarse mode split.

function optics_have_fine_fraction() result(available)

logical :: available

call gocart_optics_init()
available = have_fine_fraction

end function optics_have_fine_fraction

!----------------------------------------------------------------------
!> The optical properties of one aerosol bin at one wavelength. Resolve this once per bin,
!> outside of any loop over model levels or ensemble members.

subroutine get_bin_optics(wavelength, qty, props)

integer,               intent(in)  :: wavelength
integer,               intent(in)  :: qty
type(bin_optics_type), intent(out) :: props

integer :: ientry

call gocart_optics_init()

do ientry = 1, n_optics
   if (optics(ientry)%wavelength /= wavelength) cycle
   if (optics(ientry)%qty /= qty) cycle
   props = optics(ientry)
   return
end do

write(string1, *) 'No optical properties for quantity ', qty, ' at ', wavelength, ' nm in ' // &
                  trim(optical_properties_file)
call error_handler(E_ERR, 'get_bin_optics', string1, source, revision, revdate)

end subroutine get_bin_optics

!----------------------------------------------------------------------
!> Mass extinction efficiency [meter ** 2 / kilogram], per dry mass, at the given growth
!> factors. Pass a growth factor of 1 for a table that is not growth-factor-resolved.

subroutine bin_ext_eff(props, growth_factor, ext_eff)

type(bin_optics_type), intent(in)  :: props
real(r8),              intent(in)  :: growth_factor(:)
real(r8),              intent(out) :: ext_eff(:)

integer :: i

do i = 1, size(growth_factor)
   ext_eff(i) = interpolate_on_growth(props%growth_factor, props%ext_eff, growth_factor(i))
end do

end subroutine bin_ext_eff

!----------------------------------------------------------------------
!> The share of the bin's extinction that belongs to the fine mode, at the given growth
!> factors. Constant per wavelength for the species that do not take up water, and following
!> hygroscopic growth for the ones that do.

subroutine bin_fine_fraction(props, growth_factor, fine_fraction)

type(bin_optics_type), intent(in)  :: props
real(r8),              intent(in)  :: growth_factor(:)
real(r8),              intent(out) :: fine_fraction(:)

integer :: i

if (.not. allocated(props%fine_fraction)) then
   write(string1, *) 'The fine/coarse mode split needs the ' // FINE_FRACTION_COLUMN // &
                     ' column, which ' // trim(optical_properties_file) // ' does not have'
   call error_handler(E_ERR, 'bin_fine_fraction', string1, source, revision, revdate)
endif

do i = 1, size(growth_factor)
   fine_fraction(i) = interpolate_on_growth(props%growth_factor, props%fine_fraction, &
                                            growth_factor(i))
end do

end subroutine bin_fine_fraction

!----------------------------------------------------------------------
!> Linear interpolation of one table column onto a growth factor, clamped to the endpoints
!> outside of the table's range (the behaviour of numpy's interp, which the python operator
!> relies on).

function interpolate_on_growth(x, y, xi) result(yi)

real(r8), intent(in) :: x(:)
real(r8), intent(in) :: y(:)
real(r8), intent(in) :: xi
real(r8)             :: yi

integer :: n, lo, hi, mid

n = size(x)

if (n == 1 .or. xi <= x(1)) then
   yi = y(1)
   return
else if (xi >= x(n)) then
   yi = y(n)
   return
endif

! Bracket xi between two nodes. The shipped tables are on a uniform axis but do not rely on it.
lo = 1
hi = n
do while (hi - lo > 1)
   mid = (lo + hi) / 2
   if (xi < x(mid)) then
      hi = mid
   else
      lo = mid
   endif
end do

yi = y(lo) + (y(hi) - y(lo)) * (xi - x(lo)) / (x(hi) - x(lo))

end function interpolate_on_growth

!----------------------------------------------------------------------
!> The kappa-Koehler hygroscopic growth factor g = r_wet / r_dry at a model location, from the
!> relative humidity of the model state.
!>
!> The saturation vapour pressure follows Bolton (1980) and the Kelvin term is neglected:
!>
!>    es = es_alpha * exp(es_beta * Tc / (Tc + es_gamma))
!>    g  = (1 + kappa * RH / (1 - RH)) ** (1/3)
!>
!> RH is clamped to `max_rh` because g diverges as RH -> 1. Members whose interpolation fails
!> come back as MISSING_R8 with a non-zero istatus, in the usual forward operator convention.

subroutine get_expected_growth_factor(state_handle, ens_size, location, growth_factor, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: growth_factor(ens_size)
integer,             intent(out) :: istatus(ens_size)

real(r8), dimension(ens_size) :: qvap, tmpk, pres, tc, es, e, rh
integer,  dimension(ens_size) :: this_istatus
logical :: return_now

call gocart_optics_init()

growth_factor(:) = 1.0_r8
istatus(:) = 0

call interpolate(state_handle, ens_size, location, QTY_VAPOR_MIXING_RATIO, qvap, this_istatus)
call track_status(ens_size, this_istatus, growth_factor, istatus, return_now)
if (return_now) then
   call warn_about_humidity()
   return
endif

call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, tmpk, this_istatus)
call track_status(ens_size, this_istatus, growth_factor, istatus, return_now)
if (return_now) then
   call warn_about_humidity()
   return
endif

call interpolate(state_handle, ens_size, location, QTY_PRESSURE, pres, this_istatus)
call track_status(ens_size, this_istatus, growth_factor, istatus, return_now)
if (return_now) then
   call warn_about_humidity()
   return
endif

where (istatus == 0 .and. (tmpk <= 0.0_r8 .or. pres <= 0.0_r8))
   growth_factor = MISSING_R8
   istatus = 99
end where
if (all(istatus /= 0)) return

tc = 0.0_r8

where (istatus == 0)
   tc = tmpk - t_kelvin
   ! Bolton (1980) saturation vapour pressure over liquid water [Pa]
   es = es_alpha * exp(es_beta * tc / (tc + es_gamma))
   ! Vapour pressure from the mixing ratio [Pa]
   e  = max(qvap, 0.0_r8) * pres / (EPSILON_WATER + max(qvap, 0.0_r8))
   rh = min(max(e / es, 0.0_r8), max_rh)
   growth_factor = (1.0_r8 + kappa * rh / (1.0_r8 - rh)) ** (1.0_r8 / 3.0_r8)
end where

end subroutine get_expected_growth_factor

!----------------------------------------------------------------------
!> The table is growth-factor-resolved but the humidity is not available: say once what is
!> most likely missing, then let the usual istatus machinery reject the observations.

subroutine warn_about_humidity()

if (warned_about_humidity) return
warned_about_humidity = .true.

write(string1, *) 'Could not compute the humidity for ' // trim(optical_properties_file) // &
                  ', observations will be rejected'
write(string2, *) 'QVAPOR, THM/T, MU and PH all have to be in the model state for ' // &
                  'QTY_VAPOR_MIXING_RATIO, QTY_TEMPERATURE and QTY_PRESSURE to interpolate'
call error_handler(E_MSG, 'get_expected_growth_factor', string1, source, revision, revdate, &
                   text2=string2, text3='this message will only print once')

end subroutine warn_about_humidity

end module gocart_optics_mod
