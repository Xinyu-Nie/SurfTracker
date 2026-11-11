End-to-end U-fiber workflow
===========================

The default command is organized as::

    FreeSurfer subject + supplied FOD
        -> explicit alignment choice
           already aligned / supplied RAS affine / optional ANTs rigid fit
        -> import white/pial surfaces and original annotation IDs
        -> recommended curvature, sulcus/gyrus segmentation and HJ cleanup
        -> anatomical ROI union, seed and endpoint masks
        -> 0.5-mm inward surface and vertex-interpolated projected FODs
        -> constrained surface tracking
        -> validated tracks.tck, tracks.trk and tracks.vtk

Required inputs
---------------

* A FreeSurfer subject directory with ``mri/orig.mgz``, ``surf/lh.white``,
  ``surf/lh.pial``, ``label/lh.aparc.annot`` and corresponding RH files when
  both hemispheres are selected. ``--annotation`` selects another annotation.
  FreeSurfer executables are not needed to read these existing outputs.
* A 4D NIfTI FOD coefficient image with the declared real even-degree SH
  convention, selected with ``--sh-basis legacy-real-even``. FOD estimation
  from raw DWI is not performed. ``--sh-to-world`` supplies a separate SH
  direction transform if its axes differ from physical world RAS.
* ROI names or packed annotation IDs for targeted modes, individually or as
  explicit lists. Whole-hemisphere mode selects all cortical parcels.
  The original annotation is preserved; a union is a separate binary mask.
* Exactly one alignment choice: ``--assume-scanner-aligned``, a supplied
  ``--scanner-to-world`` forward physical RAS 4x4 affine, or ``--register-ants``.

Current recommended defaults
----------------------------

================================== =======================================
Stage / parameter                  Default
================================== =======================================
Backend                            Independent Python; no legacy executables
Workflow preprocessing profile     recommended (agreement-candidate alias)
Geometry and feature smoothing     0.5, consistent FEM mass for both
Curvature                          implicit-grid; no coordinate rounding
Grid / derivative sigma / band     1 / 1.5 / 6 mm
Gradient fallback threshold        0.2
Segmentation energy                robust-unit, regularization 1
Skeleton / label cleanup           hj-flux / pruned-branches
Endpoint flux percentile           5
Gyral / sulcal pruning             25 / 15 mm
Unsupported reconstruction limit   5 mm2
Between-gyri seeds                 connected-intersection
Within-region seeds                full intersection with merged ROI
Inward surface offset              0.5 mm
FOD interpolation                  vertex
Seeding                            area-weighted, 20 attempts/mm2, RNG 42
Tracking step / turn limit         0.5 mm / 14 degrees per mm
Angular resolution / minimum FOD   2 degrees / 0.01
Length range                       8--150 mm
Tracking threads                   4 (capped to Numba pool)
Stop at endpoint targets           enabled for both constrained modes
Additional registration            off unless explicitly selected
ANTs fit when selected             rigid, Mattes MI, one thread, COM init
================================== =======================================

The requested thread count is an upper bound: hosts with a smaller Numba pool
use that pool limit. Tracking reports retain the requested count in config
and record effective_threads and thread_pool_limit. This does not change seeded
trajectory sampling or integration parameters.

``--preprocess-profile previous-defaults`` selects the older cotan/lumped/
area-weighted/distance-skeleton baseline. ``recommended``, ``default`` and
``agreement-candidate`` select the same validated profile. Explicit flags
override profile values. Python no-config ``run_pipeline`` and ``preprocess``
use ``PreprocessConfig.recommended()``. The bare ``PreprocessConfig()``
constructor remains an explicit historical baseline for compatibility.

Already aligned FreeSurfer and FOD inputs
-----------------------------------------

Between precentral and postcentral gyri::

    swmtrack run --subject-dir /data/subject/freesurfer \
        --fod /data/subject/fod.nii.gz --output-dir output/central \
        --assume-scanner-aligned --sh-basis legacy-real-even \
        --gyrus-a precentral --gyrus-b postcentral

The temporal test uses ``--gyrus-a superiortemporal --gyrus-b middletemporal``.
Each side can also contain several disjoint ROI names/IDs. Seeds are sulcal
vertices inside their union connected through sulcal vertices to the shared
group interface. ``--between-gyri-seed-policy intersection`` retains the simple
reference option. Both endpoints must reach distinct group targets.

Within a merged parietal region::

    swmtrack run --subject-dir /data/subject/freesurfer \
        --fod /data/subject/fod.nii.gz --output-dir output/parietal \
        --assume-scanner-aligned --sh-basis legacy-real-even \
        --endpoint-policy within_gyrus \
        --gyrus-a superiorparietal inferiorparietal postcentral supramarginal precuneus

For the two-parcel intraparietal-region test use only ``superiorparietal
inferiorparietal``, or IDs ``9221140 14433500`` in standard aparc. Seeds are
all sulcal vertices inside the union, target vertices are its gyral class,
and the full union is traversable. Both same-original-parcel and cross-parcel
endpoints are accepted. No named-sulcus atlas or branch/component gate is
applied. This union does not isolate the anatomically named intraparietal sulcus.
Do not supply ``--gyrus-b`` in within-region mode.

Optional ANTs registration
--------------------------

The earlier HCP comparison did not support mandatory additional registration
for already aligned data. Use this stage only when a fit is requested. It
requires ``mri/brain.mgz`` and ANTsPy >=0.6.3, optionally in an isolated Python
environment::

    python -m venv ants-env
    ants-env/bin/python -m pip install 'antspyx>=0.6.3'

Use an existing 3D mean b0 image::

    swmtrack run --subject-dir /data/subject/freesurfer \
        --fod /data/subject/fod.nii.gz --output-dir output/registered \
        --register-ants --b0 /data/subject/mean_b0.nii.gz \
        --dwi-mask /data/subject/dwi_mask.nii.gz \
        --registration-python ants-env/bin/python \
        --sh-basis legacy-real-even --gyrus-a precentral --gyrus-b postcentral

Or replace ``--b0 ...`` with::

    --dwi /data/subject/corrected_dwi.nii.gz --bvals /data/subject/bvals

The workflow averages volumes with b-value <100 s/mm2. b-values must be finite,
nonnegative and match the volume count. DWI must already have the relevant
motion/eddy/distortion corrections; this stage does not perform them or
reconstruct FODs. The b0/DWI and optional mask must share the supplied FOD's
spatial grid and affine. A supplied mean b0 avoids scanning the 4D series.

ANTs registers FreeSurfer ``brain.mgz`` to the mean b0. It uses a supplied mask
or ANTs' intensity-based ``get_mask`` when omitted; inspect the saved mask.
``--registration-initialization header`` selects identity initialization
instead of centers of mass. ``--registration-threads`` defaults to one for
repeatability. The explicit random seed is recorded. Only rigid fitting is
exposed; no nonlinear deformation is applied.

The image pull affine is explicitly inverted and converted from LPS to a
forward RAS point transform for the surface. A native ANTs landmark check
validates the conversion. The FOD image and SH coefficients are never resampled
or rotated by this registration stage. ``--sh-to-world`` remains a separate
input convention. ANTs transform conventions are documented at
https://antspy.readthedocs.io/en/stable/registration.html .

Inspect ``registration/t1_in_b0.nii.gz`` against ``registration/mean_b0.nii.gz``
and the transformed surfaces. Successful optimization, coordinate conversion,
image bounds and surface membership do not establish anatomical alignment.

Outputs and reuse
-----------------

``run.json`` records the resolved settings, inputs/source hashes, alignment
matrix, ROI IDs, stages and final outputs. Optional registration also saves
the b0, masks, native affine, forward ``scanner_ras_to_fod_world.txt``, runtime
version and log. Each hemisphere has original/imported and preprocessed VTK
surfaces, seed/endpoint/ROI arrays, projected-FOD cache, trial diagnostics,
TCK/TRK/VTK tracks and export validation. Empty accepted outputs are reported.

Repeat the identical command with ``--resume`` to reuse hash-verified stages,
including registration. Changed inputs, settings, source code or registration
runtime require a new output directory. Existing standalone segmentation
caches can be reused by ``preprocess --reuse-segmentation`` for new ROI lists.
The low-level ``prepare``/``track`` commands continue to respect explicit masks.

Whole-hemisphere option
-------------------------------

Add ``--whole-hemisphere`` and omit ROI selection. The public workflow seeds
all cortical sulcal segmentation regions, stops both ends in the binary union
of gyral segmentation regions, applies the length limits, and classifies the
accepted trajectories by their two original endpoint parcels. Unknown and
medial-wall annotation regions and noncortex vertices are excluded. All other
recommended defaults remain unchanged. Use ``tracks.vtk`` cell labels and the
per-track endpoint tables for selection; see README.rst for label definitions.
