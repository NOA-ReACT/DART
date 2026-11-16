ICON
====

The interface implements static geometry, configurable state registration and
state locations, and full-level height observation interpolation with ensemble
evaluation and native-domain/boundary acceptance. Pressure and native history
time I/O remain unfinished. This is not yet a complete assimilation interface.

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
Scientific verification of the supplied terrain remains a separate task.

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
must be checked when implementing observation-type forward operators.

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

Auxiliary triangulation construction
-----------------------------------

``icon_triangulation_mod`` provides a serial constructor independent of model
state and MPI. It adapts ICON 2026.04's repartition/global-cap numerical path,
preserving perturbation, sorting, precision predicates and triangle order.
Required real128 precision is explicit; no lower-precision fallback is supplied.

Call ``build_auxiliary_triangulation(lon, lat, result)`` after geometry loading.
The ``auxiliary_triangulation_type`` result contains one-based ordered
``vertex(3,ntriangles)`` indices into sorted points, ``point_cell`` mapping those
points to original native cells, and unperturbed ``xyz`` coordinates. Native
cell triples are ``point_cell(vertex(:,t))``. Call ``result%clear()`` to release
storage. This constructor includes all centers, regardless of the boundary mask.
The horizontal-search lifecycle now invokes it automatically when its cache is absent.

``triangulation_check`` is a serial build target. In a compute allocation, with
``utilities_nml`` in ``input.nml``, invoke::

   triangulation_check points.txt diagnostic.txt

``points.txt`` contains the number of points, then one longitude/latitude pair
in radians per line. Alternatively::

   triangulation_check --geometry diagnostic.txt

This reads ``icon_geometry_nml`` and exercises the initialized geometry path.
The diagnostic checks repeated construction and writes the sorted mapping,
original coordinates and ordered triangles. This text format is NOT an ICON
``cc_delaunay`` file or a production cache. The production search/cache lifecycle
is described below.

Imported numerical code is separated into ``icon_delaunay_mod``,
``icon_delaunay_types_mod`` and ``icon_spherical_coordinates_mod``, with preserved
ICON attribution and the BSD-3-Clause license/author list under ``licenses``.
The serial kernel is extracted from ``triangulate_mthreaded``, not ICON's
separate ``triangulate`` implementation. A native edge-buffer bound check and
thin output structural checks reject invalid products; they do not repair them
or alter ICON's geometric thresholds.

A tiny regional synthetic grid (0.0001-degree spacing near latitude 38 degrees)
exposes invalid topology in the pinned ICON reference itself. DART rejects that
result by its structural checks. This is not a demonstrated grid-resolution
limit and no new geometric rejection threshold is introduced.

The independent oracle builder ``tests/triangulation/build_reference.py`` builds
unmodified pinned ICON Delaunay sources with serial infrastructure adapters.
It does not link DART's adapted numerical modules. Its companion
``reference_driver.f90.in`` is a test template, excluded from ordinary Fortran
source discovery. Run it on a compute node with the Intel/NetCDF environment::

   python3 tests/triangulation/build_reference.py /path/to/ICON /new/reference/build

Validation targets the controlled serial reference. Equality with multi-rank
ICON output is not established. The separate search oracle now checks the
controlled serial exact-edge selection and weights.

Horizontal search and cache
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``static_init_model`` initializes search after geometry and state registration.
Set ``model_nml:horizontal_search_cache_file`` (default
``icon_horizontal_search.nc`` in the run directory) to a shared absolute path
when cycling uses separate run directories. Task 0 loads or constructs the
complete product and broadcasts it with DART utilities. Existing compatible
caches bypass both triangulation and octree construction.

The NetCDF cache includes ordered connectivity/mapping/physical coordinates and
all native octree arrays, not just triangle bounding boxes. Loading performs
thin identity, coordinate-order, version, shape and index checks. Incompatible
files fail clearly; remove them explicitly or select another path to rebuild.
It is a DART-specific format, not native ICON ``cc_delaunay`` input.

A unique temporary sibling is closed before atomic, no-replace POSIX hard-link
publication. A competing publisher's file is loaded and checked. Save failure
warns and preserves the valid in-memory search. The cache excludes boundary
masks, vertical geometry and observation-specific data.

``icon_horizontal_search_mod`` provides ``build_horizontal_search``,
``search_candidates`` and ``query_horizontal_search``. The query returns ordered
native cell IDs, barycentric weights and a status. Its success does not imply
native-domain or LAM-boundary acceptance. The public model operator now applies
those checks and height interpolation; pressure and model time remain deferred.

For geometry-only preparation, use the MPI-capable diagnostic with one process
inside a compute allocation and the usual geometry namelist::

   printf '0\n' > no_queries.txt
   horizontal_search_check --geometry /shared/cache/search.nc no_queries.txt preparation

Nonempty query files contain a count and longitude/latitude pairs in radians.
Results are written as ``<output>.<rank>``. The synthetic single-process mode is::

   horizontal_search_check points.txt cache.nc queries.txt results [roundtrip.nc]

The independent ``tests/search/build_reference.py`` compiles pinned ICON search
routines with test infrastructure only. ``tests/search/verify_search.py`` accepts
DART build, reference build, retained triangulation fixtures, and new output
directory arguments. It checks native candidate order/selection/weights,
independent affine reproduction, exact cache round trips, failures and MPI.
See the project ``docs/horizontal_search_implementation.md`` for schema,
provenance, launch instructions and retained validation results.


Height observation operator
---------------------------

``model_interpolate`` supports configured full-level scalar quantities and
cell-centered earth-relative U/V at ``VERTISHEIGHT`` locations (metres, geometric
height). It uses independent per-column linear interpolation, then the unchanged
ICON horizontal weights. Exact full-level endpoints are accepted; there is no
extrapolation or additional upper-atmosphere cutoff. At least two full levels
are required. W remains registered at interface heights but its observation
interpolation is unsupported, as are pressure and true surface diagnostics.

Include DART's standard nearest-cell bin configuration in ``input.nml``::

   &xyz_location_nml
      nx=20, ny=20, nz=20
   /

These are performance parameters, not boundary widths. The fixed-sphere
Cartesian chord metric has the same nearest-cell ordering as arc distance.
The Cartesian nearest search now uses a distance bound to unsearched bins before
stopping; an existing premature-stop defect was exposed by the M1 global fixture.
The shared utility's great-circle and nearest-N paths are unchanged.

Native LAM membership uses a separate in-memory index over native vertex
triangles, with ICON's oriented-edge sign tolerance and no enlarged triangles.
The existing exclusion mask applies to the nearest cell and all three auxiliary
corners, even zero-weight corners. Every column must bracket the height. Neither
stencil substitution nor boundary-state freezing is performed. This index is
built once on each task and is not added to the auxiliary cache.

``icon_vertical_mod`` owns bracketing, ``icon_domain_mod`` owns domain/nearest
indices, and ``icon_interpolate_mod`` builds reusable stencils and obtains member
values through DART ``get_dart_vector_index`` / ``get_state``. DART's caller must
establish the state window. Ordinary initialized queries contain no collectives.
Missing or nonfinite values are ignored only if their final horizontal times
vertical coefficient is exactly zero. Otherwise only affected members fail;
there is no weight threshold or renormalization.

Status codes (failed member values are always ``MISSING_R8``):

* 0: success.
* 1: unsupported vertical coordinate, W, or negative quantity request.
* 2: quantity not registered in this state (including deferred pressure).
* 3: outside native LAM coverage, or no resolvable containing native triangle.
* 4: excluded nearest cell or auxiliary corner.
* 5: nearest/auxiliary search failure, including invalid auxiliary weights.
* 6: height outside any selected column, or fewer than two full levels.
* 7: required missing/nonfinite member data, or nonfinite arithmetic result.
* 8: nonfinite location coordinates.

``height_interpolation_check`` is an MPI-capable public-operator diagnostic.
``tests/height/verify_height.py BUILD OUTPUT`` generates coherent global/LAM
fixtures with sloping columns and runs single/two-task and cold/warm-cache tests.
Run inside a compute allocation with two MPI slots. See the project report
``docs/height_interpolation_implementation.md`` for equations and evidence.
