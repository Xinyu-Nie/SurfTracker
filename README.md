# SurfTracker

Surface-based probabilistic tractography for superficial white matter (SWM).
The Python pipeline starts from existing FreeSurfer outputs and an FOD image,
computes sulcus/gyrus segmentation and tracking masks, and exports U-fiber
tractograms in physical RAS+ millimeters.

| Interface | Name |
| --- | --- |
| PyPI distribution | `SurfTracker` |
| Python import | `swmtrack` |
| Command | `swmtrack` |

This branch prepares the Python package for release. It is not yet published
on PyPI. The project copyright holder remains unspecified and BSD-3-Clause is
the proposed Python license; no active project license grant is declared here.
See [licensing status](docs/licensing.rst) and [contributor credits](AUTHORS.rst).

## Install

Use a Python 3.10+ environment and install from this checkout:

```bash
python -m pip install .
swmtrack doctor
swmtrack run --help
```

After a PyPI release, the installation command will be:

```bash
python -m pip install SurfTracker
```

The Python pipeline uses NumPy, SciPy, Numba, NiBabel and meshio. It does not
require the legacy C++ binaries, MATLAB, FreeSurfer executables or a VTK runtime.
Use a fresh environment when migrating from the old development distribution
named `swmtrack`; both names install the same import package and command.

## Inputs and alignment

- A FreeSurfer subject directory with `mri/orig.mgz`, white/pial surfaces and
  the requested hemisphere's annotation (`aparc` by default).
- A supplied 4D NIfTI FOD using the supported real even-degree spherical
  harmonic convention, explicitly declared with `--sh-basis legacy-real-even`.
  FOD estimation and DWI correction are outside this pipeline. Channel counts
  alone do not establish basis compatibility.
- Exactly one alignment choice: verified existing alignment using
  `--assume-scanner-aligned`, a forward physical RAS affine using
  `--scanner-to-world`, or optional ANTs rigid registration.

See [input conventions](docs/input_conventions.rst) and the
[complete workflow and defaults](docs/end_to_end.rst). Image bounds and headers
alone do not establish anatomical registration.

## Between-gyri tracking

Precentral to postcentral:

```bash
swmtrack run --subject-dir /data/subject/freesurfer \
  --fod /data/subject/fod.nii.gz --output-dir output/central \
  --assume-scanner-aligned --sh-basis legacy-real-even \
  --gyrus-a precentral --gyrus-b postcentral
```

For superior/middle temporal U-fibers, change the selection to:

```bash
--gyrus-a superiortemporal --gyrus-b middletemporal
```

The default seeds use the connectivity-filtered intersection of sulcal
segmentation with the two ROI groups. Both endpoints must reach the distinct
group targets. Use `--between-gyri-seed-policy intersection` for the simple
intersection option. ROI names, original packed annotation IDs and lists are
accepted; packed annotation IDs are different from FreeSurfer volume LUT IDs.

## Within a merged region

For a whole-parietal ROI union:

```bash
swmtrack run --subject-dir /data/subject/freesurfer \
  --fod /data/subject/fod.nii.gz --output-dir output/parietal \
  --assume-scanner-aligned --sh-basis legacy-real-even \
  --endpoint-policy within_gyrus \
  --gyrus-a superiorparietal inferiorparietal postcentral supramarginal precuneus
```

The pipeline merges the requested parcels automatically. Seeds are the merged
ROI intersected with sulcal segmentation; endpoints are its gyral segmentation.
Both same-parcel and cross-parcel endpoints inside the union are allowed.
Original parcel IDs remain available for reporting. A two-parcel
`superiorparietal inferiorparietal` union provides an intraparietal-region
example, but does not anatomically isolate a named sulcus.

## Whole hemisphere and endpoint classification

```bash
swmtrack run --subject-dir /data/subject/freesurfer \
  --fod /data/subject/fod.nii.gz --output-dir output/whole \
  --assume-scanner-aligned --sh-basis legacy-real-even \
  --whole-hemisphere
```

This seeds cortical sulcal regions, terminates both ends in the binary union
of gyral regions, filters lengths and classifies accepted tracks by their
original endpoint parcels. Ambiguous boundary faces retain endpoint ID 0.
They remain in the full tractogram and are excluded from named parcel subsets.
Endpoint-selected subsets can leave and re-enter a selected ROI; targeted
tracking restricts the propagation domain to the requested region.

For all modes, `--hemisphere lh|rh|both` selects hemispheres (default `both`).
The default length range is 8–150 mm. The recommended curvature, seed and
tracking settings are documented in [the workflow guide](docs/end_to_end.rst).

## Optional ANTs registration

Install `.[registration]` from source. Replace `--assume-scanner-aligned` with
`--register-ants --b0 /data/subject/mean_b0.nii.gz`, or supply corrected 4D DWI
and b-values using `--register-ants --dwi ... --bvals ...`.

Registration requires `mri/brain.mgz` and an appropriate b0/DWI reference in
the FOD grid. `--registration-python` can select an isolated ANTsPy environment.
Inspect the transformed anatomy and surfaces before interpreting the tracks.
See [registration details](docs/end_to_end.rst).

## Outputs and QC

Each hemisphere has TCK/TRK/VTK tractograms, masks, reports and processing
provenance. Whole-hemisphere output additionally includes endpoint tables,
parcel-pair counts and VTK cell labels. `run.json` records inputs, parameters,
transforms and stage hashes. `--resume` reuses verified stages only for an
identical configuration and matching inputs/source.

Export full-resolution white and optional inflated label VTKs:

```bash
python scripts/export_qc.py --run-dir output/whole \
  --output-dir qc/whole --inflated
```

QC arrays include `candidate_segmentation` (0 gyrus, 1 sulcus, 2 excluded),
`seed`, `roi`, `endpoint_region`, original `anatomical_roi` IDs and source
vertex IDs. Tracks lie on a 0.5-mm inward surface: use the saved
`projection/tracking_surface.world.vtk` or transparent white surface for
physical overlays. Inflated geometry is for label QC. See [QC guidance](docs/qc.rst).

An analytic example needs no human data:

```bash
python scripts/synthetic_demo.py --output-dir output/synthetic
```

This exercises preparation, tracking and validated exports on a flat mesh;
it does not test cortical segmentation. HCP S900 examples require separate
data-use terms and recipient acceptance. See [example-data guidance](docs/example_data.rst).

## Validation and development

The portable suite creates synthetic fixtures. The scientific validation covers
four subjects, both hemispheres and the central, temporal, merged-parietal and
whole-hemisphere workflows. A small sensitivity screen found subject-dependent
temporal responses to offset and FOD threshold. Counts and numerical checks do
not establish anatomical correctness or complete sulcal coverage.
See [validation scope](docs/validation.rst).

```bash
python -m pip install '.[test,release]'
python -m build
python -m twine check --strict dist/*
python scripts/check_distributions.py dist/*.whl dist/*.tar.gz
python scripts/verify_artifacts.py --wheel dist/*.whl --sdist dist/*.tar.gz \
  --output-dir release-audit/artifact-tests
```

[GitHub Actions](.github/workflows/ci.yml) checks source/wheel builds, installed
artifacts outside the checkout and the synthetic example across the configured
Python/platform matrix. It does not publish packages. See [testing](docs/testing.rst)
and [release preparation](docs/publishing.rst).

## Legacy C++ implementation

The original [C++ implementation](Surface_tracking_c++/README.md) is retained
unchanged for reference, including its contributor and third-party notices.
It is independent of the Python installation and excluded from Python wheel
and source distributions. Historical C++ documentation may reference the MATLAB
implementation, which has been removed from the current repository tree.
Earlier revisions remain in Git history.

The existing `101107.tar.gz` Git LFS example pointer is retained for legacy
compatibility and is excluded from the Python distributions. Consult the
applicable HCP data-use terms before obtaining or using human example data.
The proposed Python license does not relicense legacy source or example data.

## Cite both methods

Please cite **both** articles when using this pipeline. If you prepare the ROIs and WM surface by yourself, you can only cite 2:

1. Yonggang Shi, Paul M. Thompson, Ivo Dinov and Arthur W. Toga.
   *Hamilton-Jacobi Skeleton on Cortical Surfaces.* IEEE Transactions on Medical
   Imaging 27(5):664–673 (2008). [DOI: 10.1109/TMI.2007.913279](https://doi.org/10.1109/TMI.2007.913279).
2. Xinyu Nie, Jialiang Ruan, Maria Concepción Garcia Otaduy, Lea Tenenholz Grinberg,
   John Ringman and Yonggang Shi. *Surface-Based Probabilistic Fiber Tracking in
   Superficial White Matter.* IEEE Transactions on Medical Imaging
   43(3):1113–1124 (2024). [DOI: 10.1109/TMI.2023.3329451](https://doi.org/10.1109/TMI.2023.3329451).

[CITATION.cff](CITATION.cff) and [CITATIONS.bib](CITATIONS.bib) provide citation
metadata. [AUTHORS.rst](AUTHORS.rst) records contributor roles and AI-assisted
development separately. [THIRD_PARTY_NOTICES.rst](THIRD_PARTY_NOTICES.rst)
describes the scope of dependency and legacy notices.
