ICON
====

Overview
--------

DART interface to the `ICON <https://www.icon-model.org>`__ model: icosahedral
triangular grid, terrain following vertical coordinate in height (SLEVE). Written
for limited-area runs with the NWP physics, and driven by
`WRF-Ensembly <https://github.com/NOA-ReACT/wrf_ensembly>`__. Global grids work as
well, nests do not.

The state is read from ICON output files on the native grid, so every state variable
is on cell centres with the netCDF dimensions ``(time, height, ncells)`` (full
levels), ``(time, height_2, ncells)`` (half levels, e.g. ``w``) or ``(time, ncells)``.
The names of the dimensions do not matter, their order and sizes do.

How the analysis gets back into ICON
------------------------------------

ICON's prognostic wind is the component normal to the cell edges (``vn``), and its
thermodynamic variables are ``rho`` and ``theta_v``. None of these go into the DART
state. Instead, use the diagnostic ``u``, ``v``, ``temp`` and ``pres`` on cells, which
ICON writes with its first guess (output group ``dwd_fg_atm_vars``), and start the
next forecast from the first guess plus an analysis file with these four variables
(``initicon_nml``: ``init_mode = 1``, or ``5`` for an incremental analysis update).
ICON then turns the difference in ``u`` and ``v`` into an increment of ``vn`` and
rebuilds ``rho`` and ``theta_v`` from the analysed temperature and pressure. The
moisture variables can go into the same analysis file.

Files
-----

``template_file``
   An ICON output file with all state variables, e.g. a first guess file.

``grid_file``
   The ICON grid file. ``clon``, ``clat``, ``vlon``, ``vlat``, ``vertex_of_cell``,
   ``cells_of_vertex`` and, for limited-area grids, ``refin_c_ctrl`` are read.

``height_file``
   An ICON output file with the half level heights ``z_ifc``. It can be the same file
   as ``template_file``. Full level heights are the mean of the two half levels
   around them.

Interpolation
-------------

Horizontal
   The cell centre nearest to the location is found first. The cells around each of
   its three vertices, together with the cell itself, make a fan of triangles of cell
   centres. The location is interpolated linearly (barycentric weights) in the
   triangle that contains it. On the rim of a limited-area grid, where a location can
   be inside a cell but not between three cell centres, the value of the cell is used.

Vertical
   Linear in height between the two full levels around the location, in each of the
   three columns. For locations in pressure, the levels are found in each member's own
   pressure column (linear in the logarithm of pressure), which needs a state variable
   with ``QTY_PRESSURE``. Locations in the lowest half layer get the value of the
   lowest level. Locations below the model surface or above the highest full level
   are rejected.

Model levels
   ICON counts levels from the model top. Locations with ``VERTISLEVEL`` count from
   the ground, with the lowest level being both 0 and 1, and return a status of 2 at
   or above the top level. This is what the column operators (AOD) expect, as in the
   WRF interface.

Surface
   Locations with ``VERTISSURFACE`` use the screen level diagnostics if they are in
   the state (``QTY_2M_TEMPERATURE``, ``QTY_2M_SPECIFIC_HUMIDITY``,
   ``QTY_10M_U_WIND_COMPONENT``, ``QTY_10M_V_WIND_COMPONENT``,
   ``QTY_SURFACE_PRESSURE``), otherwise the lowest model level.

Quantities that are computed when they are not in the state:

==============================  =====================================================
Quantity                        From
==============================  =====================================================
``QTY_GEOPOTENTIAL_HEIGHT``     height of the levels (also ``QTY_GEOMETRIC_HEIGHT``)
``QTY_SURFACE_ELEVATION``       lowest half level
``QTY_VAPOR_MIXING_RATIO``      ``QTY_SPECIFIC_HUMIDITY`` (ICON's ``qv``)
``QTY_SPECIFIC_HUMIDITY``       ``QTY_VAPOR_MIXING_RATIO``
``QTY_DENSITY``                 pressure, temperature and specific humidity
``QTY_POTENTIAL_TEMPERATURE``   pressure and temperature
==============================  =====================================================

Localisation
------------

Every state element has a height, and the vertical localisation is always in height
(``VERTISHEIGHT``). Observation locations are converted to height once, before the
assimilation: keep ``convert_all_obs_verticals_first`` and
``convert_all_state_verticals_first`` at ``.true.`` in ``&assim_tools_nml``.

Namelist
--------

.. code-block:: fortran

   &model_nml
      template_file               = 'icon_template.nc'
      grid_file                   = 'icon_grid.nc'
      height_file                 = 'icon_heights.nc'
      assimilation_period_days    = 0
      assimilation_period_seconds = 21600
      boundary_width              = 4
      sfc_elev_max_diff           = -1.0
      debug                       = 0
      variables = 'u',        'QTY_U_WIND_COMPONENT',  'NA',  'NA', 'UPDATE',
                  'v',        'QTY_V_WIND_COMPONENT',  'NA',  'NA', 'UPDATE',
                  'temp',     'QTY_TEMPERATURE',       'NA',  'NA', 'UPDATE',
                  'pres',     'QTY_PRESSURE',          'NA',  'NA', 'UPDATE',
                  'qv',       'QTY_SPECIFIC_HUMIDITY', '0.0', 'NA', 'UPDATE',
                  'pres_sfc', 'QTY_SURFACE_PRESSURE',  'NA',  'NA', 'NO_COPY_BACK',
   /

``variables``
   Five columns per state variable: netCDF name, DART quantity, minimum, maximum
   (``NA`` for none) and ``UPDATE`` or ``NO_COPY_BACK``.

``boundary_width``
   Observations whose nearest cell is within this many rows of cells from the lateral
   boundary of a limited-area grid are rejected (``refin_c_ctrl`` of the grid file).
   0 keeps all of them.

``sfc_elev_max_diff``
   Surface observations are rejected if the model surface is more than this many
   metres away from the height of the observation. Negative: no check.

Return codes of ``model_interpolate``
-------------------------------------

==  ===========================================================
 2  model level at or above the top (``VERTISLEVEL``)
11  outside of the grid
12  in the boundary zone of a limited-area grid
15  below the model surface
16  above the highest full level
17  surface elevation too different (``sfc_elev_max_diff``)
18  vertical coordinate does not fit the variable
19  location in pressure, but no pressure in the state
20  quantity is not in the state and cannot be computed
21  unknown vertical coordinate
22  the interpolated value is not a number
==  ===========================================================
