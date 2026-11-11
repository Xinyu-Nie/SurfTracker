SurfTracker: superficial white matter tracking
==============================================

SurfTracker uses the PyPI distribution name ``SurfTracker`` and provides the Python
import package and command ``swmtrack``. It runs from existing FreeSurfer
outputs and a supplied FOD image to sulcus/gyrus segmentation, surface tracking
and TCK/TRK/VTK tractograms in physical RAS+ millimeters. The default pipeline
uses independent NumPy/SciPy/Numba implementations. It requires no MATLAB,
FreeSurfer executable, VTK runtime or separately supplied legacy executables.
Optional rigid registration uses ANTsPy. FOD estimation is outside this pipeline.

This is the 0.3.0.dev17 release-preparation tree. The numerical implementation
and recommended defaults match dev11; this revision separates portable tests
and documentation from the institution-specific research workspace. The
legal copyright-holder and license authorization still need to be finalized
before publication. Human contributor credits and both citations are prepared. See ``RELEASE_PREPARATION.rst``.

Install
-------

After publication on PyPI, install in a Python 3.10+ environment::

    python -m pip install SurfTracker
    swmtrack doctor

For this release-preparation source tree, install from this directory::

    python -m pip install .
    swmtrack doctor
    swmtrack run --help

For optional ANTs registration, install ``.[registration]`` from source
(or ``SurfTracker[registration]`` after publication), or select a separate
ANTsPy interpreter with ``--registration-python``. Linux/Python 3.12 is the
environment used for scientific cohort validation; GitHub Actions checks the
configured compatibility matrix. The pure-Python wheel tag alone is not evidence
that all dependency/platform combinations have been tested.

Inputs
------

* A FreeSurfer subject containing ``mri/orig.mgz``, white and pial surfaces,
  and the requested hemisphere's annotation file (``aparc`` by default).
* A 4D NIfTI FOD with the supported even-degree real-SH convention, explicitly
  declared through ``--sh-basis legacy-real-even``. Matching channel counts
  alone do not establish SH basis compatibility. See ``docs/input_conventions.rst``.
* An explicit alignment choice: verified existing scanner/FOD alignment,
  a supplied physical RAS affine, or optional ANTs rigid registration using
  an appropriate b0/DWI reference. Headers alone do not establish registration.

MRtrix3 FODs and orientation
----------------------------

Standard MRtrix3 WM FODs from dwi2fod csd or msmt_csd use the supported
orthonormal real even-degree SH basis. Select --sh-basis legacy-real-even;
this existing option name refers to the compatible basis, not the older
MRtrix 0.2 normalization. Use the WM FOD rather than GM/CSF tissue maps.
SurfTracker reads NIfTI, so export a .mif FOD with MRtrix3 first::

    mrconvert wmfod.mif wmfod.nii.gz -datatype float32

An existing MRtrix3 NIfTI FOD can be supplied directly. Keep all SH channels;
format conversion does not register anatomy or change an incompatible SH basis.
See the `MRtrix3 SH convention <https://mrtrix.readthedocs.io/en/latest/concepts/spherical_harmonics.html>`_
and `mrconvert documentation <https://mrtrix.readthedocs.io/en/latest/reference/commands/mrconvert.html>`_.

An existing 101107 LH experiment tested MRtrix3 3.0.8 MSMT-CSD WM FODs
(lmax=8), raw and mtnormalise-normalized. All 45 basis functions matched
MRtrix3 sh2amp on 187 directions (maximum absolute error 5.70e-8). Tracking
and export checks passed using existing surfaces/masks; this was a separate
tracker compatibility experiment. See docs/validation.rst for its exact scope.

No particular voxel storage orientation is required. The full NIfTI affine
supports RAS/LAS/LPS storage and oblique grids. Surface coordinates are converted
from FreeSurfer surface RAS via orig.mgz to scanner RAS and then into FOD-world
RAS using the declared alignment. Physical output coordinates are RAS+ mm;
the FOD image is not resampled into a canonical RAS voxel array.

Standard MRtrix3 FOD directions use scanner/world RAS, so normally omit
--sh-to-world (identity), regardless of voxel storage order. This separate
orthogonal 3x3 matrix is only needed for a known different SH direction frame;
it does not convert a basis or register anatomy. Matching orientation labels
are not anatomical registration: choose verified --assume-scanner-aligned,
--scanner-to-world or optional --register-ants. Registration moves surfaces
into FOD space and leaves the FOD unchanged. A separately rotated/warped FOD
requires appropriate FOD reorientation. See docs/input_conventions.rst and the
`MRtrix3 coordinate convention <https://mrtrix.readthedocs.io/en/latest/getting_started/image_data.html#coordinate-system>`_.

Tracking examples
-----------------

Between precentral and postcentral gyri::

    swmtrack run --subject-dir /data/subject/freesurfer \
        --fod /data/subject/fod.nii.gz --output-dir output/central \
        --assume-scanner-aligned --sh-basis legacy-real-even \
        --gyrus-a precentral --gyrus-b postcentral

For superior--middle temporal tracking, use ``--gyrus-a superiortemporal
--gyrus-b middletemporal``. Between-gyri seeding defaults to the
connectivity-filtered sulcal intersection; ``--between-gyri-seed-policy
intersection`` retains the simple intersection option.

Within a merged parietal region::

    swmtrack run --subject-dir /data/subject/freesurfer \
        --fod /data/subject/fod.nii.gz --output-dir output/parietal \
        --assume-scanner-aligned --sh-basis legacy-real-even \
        --endpoint-policy within_gyrus \
        --gyrus-a superiorparietal inferiorparietal postcentral supramarginal precuneus

Names, original packed annotation IDs and mixed lists are accepted. For example,
``--gyrus-a 9221140 14433500 --endpoint-policy within_gyrus`` merges the standard
aparc superiorparietal and inferiorparietal parcels. These are annotation IDs,
not FreeSurfer volume LUT IDs. Seeds are sulcal segmentation intersected with
the union; endpoints are gyral segmentation intersected with the same union.
Original parcel IDs remain available. No named-sulcus atlas is applied, so
this two-parcel union does not isolate the anatomically named intraparietal sulcus.

Whole-hemisphere tracking with endpoint classification::

    swmtrack run --subject-dir /data/subject/freesurfer \
        --fod /data/subject/fod.nii.gz --output-dir output/whole \
        --assume-scanner-aligned --sh-basis legacy-real-even \
        --whole-hemisphere

Whole-hemisphere mode seeds all cortical sulcal regions and uses one binary
union of gyral regions as the termination target. Both ends must belong to
that target. Unknown and medial-wall regions are excluded, with the FreeSurfer
cortex mask applied when present. Omit ``--gyrus-a/b``. Default retained lengths
are 8--150 mm; use ``--min-length-mm`` and ``--max-length-mm`` to change them.
Endpoint-selected subsets can leave and re-enter their selected annotation
regions; targeted modes instead restrict the propagation domain.

For all modes, ``--hemisphere lh|rh|both`` selects the hemispheres (default both).
The complete defaults and optional ANTs/b0/DWI workflow are documented in
``docs/end_to_end.rst``. More ready-to-adapt commands are in ``examples/README.rst``.

Outputs and QC
--------------

Each hemisphere's ``tractography`` directory contains ``tracks.tck``,
``tracks.trk``, ``tracks.vtk``, a JSON report and per-attempt diagnostics.
``run.json`` records input/source hashes, parameters, transforms and stages.
The preprocessing directory retains the sulcal/gyral segmentation, original
annotation, source vertex IDs and binary seed/ROI/endpoint masks on the full
white surface. Repeat the identical command with ``--resume`` to reuse
hash-verified completed stages. A changed input/configuration/code requires
a new output directory.

Whole-hemisphere mode adds per-track endpoint CSV/NPZ tables, unordered
parcel-pair counts and a symmetric count matrix. Tract VTK cell scalars include
``endpoint_a_id``, ``endpoint_b_id``, ``pair_low_id``, ``pair_high_id``,
``length_mm``, ``accepted_trial_id`` and ``endpoint_class`` (1 same parcel,
2 cross parcel, 3 ambiguous). A terminal face whose gyral vertices disagree
on parcel ID receives ID 0 (ambiguous); the full output retains that track.
TCK/TRK carry geometry; the accompanying tables carry endpoint metadata.

For additional full white/inflated label VTKs, this source tree includes::

    python scripts/export_qc.py --run-dir output/whole \
        --output-dir qc/whole --inflated

The helper uses the recorded FreeSurfer subject directory, or an explicit
``--subject-dir`` override, to locate optional inflated surfaces. It preserves
source vertex identity and checks every exported point/face/array by readback.
Point data include ``candidate_segmentation`` (0 gyrus, 1 sulcus, 2 outside
tracking domain), ``seed``, ``roi``, ``endpoint_region`` and original
``anatomical_roi`` IDs, with RGB arrays. See ``docs/qc.rst`` for viewing details.

Tracks lie on a 0.5-mm inward surface. Overlay them with the saved
``projection/tracking_surface.world.vtk`` or a transparent white surface.
Inflated surfaces are for label QC only. Computational checks and a nonempty
bundle do not establish anatomical accuracy or correct input registration.

Build and test
--------------

The portable suite generates its own fixtures; it requires no HCP subjects,
FOD downloads, institutional directories, legacy binaries or research drivers::

    python -m pip install '.[test,release]'
    OPENBLAS_NUM_THREADS=1 python -m pytest -q
    python -m build
    python scripts/check_distributions.py dist/*.whl dist/*.tar.gz
    python -m twine check dist/*

The archive checker rejects native/data/cache payloads and mismatched source
versions. ``--require-license`` additionally checks for project license text
and metadata once the approved license is added. The helpers do not publish
anything. ``docs/testing.rst`` describes tests of extracted sdists and installed
wheels outside the source checkout. ``docs/validation.rst`` summarizes the
existing numerical/cohort evidence and its limits.

Primary citations
-----------------

When using the pipeline in scientific work, please cite BOTH articles and report
the software version and processing settings:

* Yonggang Shi, Paul M. Thompson, Ivo Dinov and Arthur W. Toga.
  Hamilton-Jacobi Skeleton on Cortical Surfaces. IEEE Transactions on Medical
  Imaging 27(5), 664--673 (2008). https://doi.org/10.1109/TMI.2007.913279 .
* Xinyu Nie, Jialiang Ruan, Maria Concepción Garcia Otaduy, Lea Tenenholz Grinberg,
  John Ringman and Yonggang Shi. Surface-Based Probabilistic Fiber Tracking in
  Superficial White Matter. IEEE Transactions on Medical Imaging 43(3),
  1113--1124 (2024). https://doi.org/10.1109/TMI.2023.3329451 .

CITATION.cff records both as primary references; CITATIONS.bib supplies BibTeX.
The software citation credits Xinyu Nie and Yonggang Shi; AUTHORS.rst records
Jialiang Ruan's C++ contributions and AI-assisted development separately. No
software DOI or publication date is invented. The journal years above follow the
final issues; the earlier DOI years are retained exactly.

Licensing preparation
---------------------

BSD-3-Clause is recommended. Its proposed standard text is in the local licensing/
review folder; the legal holder and license authorization remain unresolved.
No active grant is declared yet. Research credits are in AUTHORS.rst and direct
dependency terms are summarized in THIRD_PARTY_NOTICES.rst. See docs/licensing.rst
and docs/ownership_review.rst for the rationale and institutional context.

Example data and publishing
---------------------------

Run ``python scripts/synthetic_demo.py --output-dir output/synthetic`` for an
analytic prepare/track/export example with no human data. It does not exercise
FreeSurfer import or cortical segmentation. ``docs/example_data.rst`` explains
conditional HCP S900 redistribution and separate data-use terms. Human example
data remain outside the Python distributions.

The prepared GitHub Actions workflow checks build artifacts on a proposed
Linux/Python 3.10--3.13 and macOS/Windows/Python 3.12 matrix. See the repository Actions page for current results. ``docs/publishing.rst`` covers TestPyPI rehearsal and Trusted
Publishing setup; no automatic publishing job is enabled.

Repository layout
-----------------

README.md is the GitHub landing page. The main branch contains the Python
pipeline, portable tests, documentation and synthetic examples. Historical
implementations, their original notices and the old example archive remain on
the `legacy branch <https://github.com/Xinyu-Nie/SurfTracker/tree/legacy>`_ and
`dev-c++ branch <https://github.com/Xinyu-Nie/SurfTracker/tree/dev-c%2B%2B>`_,
with earlier revisions in Git history. No proposed Python license applies to
those historical assets. The external native-kernel backend is no longer
supported; the independent Python pipeline and its recommended defaults remain.
