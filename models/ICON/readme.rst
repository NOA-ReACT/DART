ICON
====

The interface implements static geometry, configurable state registration and
state locations. Pressure, observation interpolation and native history time I/O
remain unfinished. This is not yet a working assimilation interface.

Geometry initialization
-----------------------

``icon_geometry_mod`` owns a read-only geometry object and provides idempotent
``init_icon_geometry`` and ``end_icon_geometry`` routines. ``model_mod`` invokes
them during initialization and shutdown. A separate ``geometry_check`` program
exercises geometry without requiring history fields or state registration.

Place this namelist in the experiment working directory's ``input.nml``::

   &icon_geometry_nml
      horizontal_grid_file = 'icon_grid.nc',
      vertical_geometry_file = 'icon_geometry.nc',
      nudge_zone_width = 10,
   /

Use paths to the actual files or local symlinks. The filenames above are defaults;
10 is an example, not a default. For a limited-area domain the user must supply
the effective ICON nudging width in cell rows. Global grids do not use this value.
No width is inferred from metadata. The specified zone is fixed at four cell
rows, matching ICON 2026.04 ``mo_impl_constants_grf``. The observation exclusion
mask is ``1 <= refin_c_ctrl <= 4 + nudge_zone_width``. It is computed once in
memory; no ``boundary_mask.nc`` is read or written. This mask must not suppress
filter updates to boundary state elements.

The native horizontal grid supplies ``clon``, ``clat``, ``vlon``, ``vlat`` in
radians, ``vertex_of_cell``, ``neighbor_cell_index``, and (for LAM)
``refin_c_ctrl``. Required global metadata are ``uuidOfHGrid``, ``global_grid``
and positive ``semi_major_axis``. Native dimensions are ``cell``, ``vertex``,
``nv=3``; connectivity is one-based and missing neighbors are preserved as -1.
Nested-grid boundary mappings are outside this first implementation.

The complementary static file supplies ``z_mc(full_level,cell)`` and
``z_ifc(interface_level,cell)`` in CDL order, in metres, with no time dimension.
The Fortran arrays use (cell,level) order. Axis names are taken from each variable,
not inferred from a ``height`` suffix; levels must number K and K+1. Native
single/double precision is accepted. Horizontal UUID and cell count must match;
``uuidOfVGrid`` is retained for future state-file checks. Original cell ordering
is required; UUID comparison cannot detect arbitrary external permutations.

Checks cover rank/order/sizes, units, identities, finite/nonmissing heights,
connectivity bounds and valid boundary-row coverage. Full levels must lie
strictly between their interfaces, ordered top to bottom. Cartesian cell centers
are derived as unit vectors. Flat columns are allowed: validation does not
establish that a file contains the intended experiment's terrain.

Building and testing
--------------------

``work/quickbuild.sh`` includes ``geometry_check`` as a serial model program.
The local HPC wrapper also supports it, preserving separate build products::

   sbatch ~/software/config/DART/build-dart.sbatch icon-oneapi-geometry-02

Use a fresh build label. After a successful build, load its retained
``environment.sh`` and run ``bin/geometry_check`` from a separate experiment
directory containing the geometry namelist. Use the site's compute allocation
for execution. The program reports geometry counts/ranges, checks repeat
initialization and cleanup/reinitialization, then prints ``geometry_check PASS``.

Synthetic integration tests require Python 3 and ``ncgen``, with no Python
scientific packages. They create their own temporary files::

   python3 ~/software/src/DART/models/ICON/tests/test_geometry.py \
      /absolute/build/path/bin/geometry_check

An optional ``--work-dir /new/experiment/test-directory`` retains fixtures and
logs. Tests cover LAM/global masks, zero nudging, axis names, native float heights,
and rejection of incompatible UUIDs, dimensions, indices, units, missing/nonfinite
heights, invalid vertical bracketing and insufficient boundary-row coverage.
Scientific terrain verification and interpolation remain separate milestones.

State indexing and metadata
---------------------------

``icon_state_mod`` validates ICON field semantics and supplies locations.
DART's ``parse_variables`` and ``add_domain`` own state registration and file
metadata; ``get_model_variable_indices`` supplies the cell, level, variable and
quantity. No separate flattened-index map or custom state reader is introduced.
This follows the MPAS separation between DART indexing and model locations;
no MPAS source was copied.

``model_nml`` uses the standard DART name/quantity/copy-back table::

   &model_nml
      template_file = 'icon_history.nc',
      state_variables = 'u',    'QTY_U_WIND_COMPONENT',       'UPDATE',
                        'v',    'QTY_V_WIND_COMPONENT',       'UPDATE',
                        'w',    'QTY_VERTICAL_VELOCITY',      'UPDATE',
                        'temp', 'QTY_TEMPERATURE',            'UPDATE',
                        'qv',   'QTY_SPECIFIC_HUMIDITY',       'UPDATE',
                        'qc',   'QTY_CLOUDWATER_MIXING_RATIO', 'UPDATE',
                        'qi',   'QTY_ICE_MIXING_RATIO',        'UPDATE',
   /

An omitted/all-blank table selects these seven defaults. An explicitly supplied
nonempty table replaces the defaults completely: deleting a row removes that
field, including a trailing field. Rows may be reordered. Each row is a complete
triple with no gaps. Native NetCDF field names are used, not U/T aliases.
``NO_COPY_BACK`` is DART's output-write flag, not a promise to freeze that field
inside the filter. Pressure fields are deliberately rejected for this milestone.

Additional supported triples use these names and quantities (all ``UPDATE``
unless explicitly changed)::

   qr   QTY_RAINWATER_MIXING_RATIO
   qs   QTY_SNOW_MIXING_RATIO
   qg   QTY_GRAUPEL_MIXING_RATIO
   qh   QTY_HAIL_MIXING_RATIO
   qnc  QTY_DROPLET_NUMBER_CONCENTR
   qni  QTY_ICE_NUMBER_CONCENTRATION
   qnr  QTY_RAIN_NUMBER_CONCENTR
   qns  QTY_SNOW_NUMBER_CONCENTR
   qng  QTY_GRAUPEL_NUMBER_CONCENTR
   qnh  QTY_HAIL_NUMBER_CONCENTR

These use DART's existing hydrometeor quantity identifiers. Stored values retain
ICON conventions: hydrometeor mass fractions in kg kg-1, number fields in kg-1,
and qv as specific humidity. There is no implicit conversion to dry-air mixing
ratios or number per volume. These quantity labels do not by themselves establish
compatibility with a forward operator; operator-specific definitions/conversions
must be checked when interpolation/operators are implemented.

The template must have matching grid UUIDs, complete cell/level dimensions, and
one unlimited ``time`` record. Selected fields must be floating-point arrays
with CDL order (time,level,cell). W has K+1 interfaces; the other supported fields
have K full levels. Units are checked (surrounding whitespace is ignored).
Unknown/duplicate fields, inconsistent quantities, incorrect dimensions or units
fail clearly. Extra unselected fields are ignored. Metadata initialization does
not scan all state values or yet decode the time value. Use DART state ingestion
and the forthcoming ICON time reader for those responsibilities.

State longitude/latitude use native cell centers converted to degrees, with
longitude wrapped into [0,360). W uses ``z_ifc(cell,level)``; all other fields use
``z_mc(cell,level)``. Every location uses ``VERTISHEIGHT``. Boundary cells retain
valid state locations. W metadata is available even though W interpolation is
deferred. ``end_model`` releases model arrays; DART owns the domain registry,
so model shutdown is terminal for this instance.

``state_metadata_check`` is another serial model target. Its checks exercise
``model_mod`` initialization and metadata, index inversion, field ordering,
quantities, full/interface heights and negative identity-observation indices.
Small fixtures check every element; large real inputs sample first/last cells
on every level of every selected field. Run the synthetic tests with::

   python3 ~/software/src/DART/models/ICON/tests/test_state_metadata.py \
      /absolute/build/path/bin/state_metadata_check

A standalone diagnostic input also needs the standard DART ``utilities_nml``,
``obs_kind_nml`` and ``location_nml`` entries (empty entries with the slash on a
separate line are sufficient), plus ``icon_geometry_nml`` and ``model_nml``.
The supplied ``work/input.nml`` includes the normal DART namelists.
