# SurfTracker

Surface-based probabilistic tractography for superficial white matter (SWM).
The Python pipeline starts from existing FreeSurfer outputs and an FOD image,
computes sulcus/gyrus segmentation and tracking masks, and exports U-fiber
tractograms in physical RAS+ millimeters.

| Interface | Name |
| --- | --- |
| Distribution | `SurfTracker` |
| Python import | `swmtrack` |
| Command | `swmtrack` |

This public repository provides documentation, usage examples and compiled
binary downloads. The scientific core is maintained privately. The `swmtrack`
command, Python imports and recommended numerical settings are preserved.

## Download and install

The current build is **0.3.0.dev18 for CPython 3.12 on Linux x86_64**.
Windows, macOS, other Python versions and other architectures need separate
builds. The native wheel has manylinux2014 / glibc 2.17 compatibility tags;
execution was tested on a glibc 2.31 host.

- [Download the complete binary kit](binaries/0.3.0.dev18/SurfTracker-0.3.0.dev18-linux-cp312-binary.zip): wheel, pinned validation dependencies, examples and documentation.
- [Download the wheel](binaries/0.3.0.dev18/surftracker-0.3.0.dev18-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.whl).
- [Checksums](binaries/0.3.0.dev18/SHA256SUMS) · [Download manifest](binaries/0.3.0.dev18/manifest.json).

After unzipping the kit, create a fresh Python 3.12 environment and follow its
README. To install directly from GitHub:

```bash
python3.12 -m venv surftracker-env
source surftracker-env/bin/activate
python -m pip install \
  -c https://raw.githubusercontent.com/Xinyu-Nie/SurfTracker/main/binaries/0.3.0.dev18/requirements-validation.txt \
  https://raw.githubusercontent.com/Xinyu-Nie/SurfTracker/main/binaries/0.3.0.dev18/surftracker-0.3.0.dev18-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.whl
swmtrack doctor
swmtrack run --help
```

The constraints reproduce the numerical dependency versions used for validation.
No compiler or Cython runtime is required. SurfTracker uses NumPy, SciPy,
NiBabel and meshio; Numba remains a dependency for thread-configuration
compatibility. Core algorithms execute in native modules. Existing FreeSurfer
outputs can be read without FreeSurfer executables, MATLAB or a VTK runtime.

Use the supplied wheel for installation. This documentation/download repository
contains no source package to build with `pip install .`. Use a fresh environment
when migrating from the older distribution named `swmtrack`, which uses the same
import package and command. See [download details](docs/downloads.rst).

## Inputs and alignment

- A FreeSurfer subject directory with `mri/orig.mgz`, white/pial surfaces and
  the requested hemisphere's annotation (`aparc` by default).
- A supplied 4D NIfTI WM FOD, including standard MRtrix3 FODs, using the
  supported real even-degree spherical-harmonic convention. Declare it with
  `--sh-basis legacy-real-even`; see the export example below. FOD estimation
  and DWI correction are outside this pipeline. 
- Exactly one alignment choice: verified existing alignment using
  `--assume-scanner-aligned`, a forward physical RAS affine using
  `--scanner-to-world`, or optional ANTs rigid registration.

See [input conventions](docs/input_conventions.rst) and the
[complete workflow and defaults](docs/end_to_end.rst). Image bounds and headers
alone do not establish anatomical registration.

### MRtrix3 FODs

WM FODs produced by MRtrix3 `dwi2fod csd` or `dwi2fod msmt_csd` use the
supported SH basis: even degrees, ascending `m = -l ... l`, negative-`m` sine
and positive-`m` cosine terms, with orthonormal scaling and the Condon–Shortley
phase. Common channel counts are 6, 15, 28 and 45 for `lmax` 2, 4, 6 and 8.
The existing option name `legacy-real-even` selects this compatible basis.

SurfTracker reads NIfTI FODs (`.nii` or `.nii.gz`), not `.mif` files directly.
Export the **WM FOD**, keeping all SH channels; for multi-tissue CSD, use the WM
output rather than the GM/CSF maps:

```bash
# Format conversion with MRtrix3; no SH-basis conversion is needed.
mrconvert wmfod.mif wmfod.nii.gz -datatype float32

# Use --assume-scanner-aligned only after verifying anatomical alignment.
swmtrack run --subject-dir /data/subject/freesurfer \
  --fod wmfod.nii.gz --output-dir output/central \
  --assume-scanner-aligned --sh-basis legacy-real-even \
  --gyrus-a precentral --gyrus-b postcentral
```

### Orientation and physical coordinates

**No specific voxel storage orientation is required.** RAS-, LAS-, LPS-ordered
and oblique NIfTI grids are sampled using the full header affine. SurfTracker
works in physical **RAS+ millimeters** (right, anterior, superior), without
resampling or rewriting the FOD into a canonical RAS voxel array.

- **FreeSurfer surfaces:** the importer converts surface/tkregister RAS to
  scanner RAS using `mri/orig.mgz`, then applies the selected scanner-to-FOD
  transform. Keep the original FreeSurfer subject files together.
- **Anatomical alignment:** use `--assume-scanner-aligned` only when scanner
  RAS already aligns with FOD-world RAS; otherwise supply `--scanner-to-world`
  (a forward 4×4 physical RAS point transform) or request `--register-ants`
  with a b0/DWI reference in the FOD grid. RAS coordinates alone do not imply
  that two datasets are anatomically registered.
- **FOD directions:** standard MRtrix3 FOD directions use scanner/world RAS
  axes, so normally omit `--sh-to-world` (identity). Voxel axis order does not
  change that convention: do not apply the image-affine rotation a second time.
  For externally generated coefficients in a different direction frame,
  `--sh-to-world` accepts an orthogonal 3×3 matrix mapping those directions
  into FOD-world RAS. This matrix does not convert the SH basis or register
  the anatomy. See [MRtrix3 coordinates](https://mrtrix.readthedocs.io/en/latest/getting_started/image_data.html#coordinate-system).

The pipeline transforms the surfaces into the supplied FOD space, keeping the
FOD/SH image unchanged, including during optional ANTs registration. For a
separately rotated or warped FOD image, use FOD-aware reorientation before
tracking; editing a header or resampling coefficient channels as scalar images
alone does not rotate the lobes. See
[MRtrix3 FOD reorientation](https://mrtrix.readthedocs.io/en/latest/reference/commands/mrtransform.html).

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

Install ANTsPy in the current or a separate environment as described in the
[workflow guide](docs/end_to_end.rst). Replace `--assume-scanner-aligned` with
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

## Python interface

The same Python calls are available after installing the wheel:

```python
from swmtrack.tracking import TrackingConfig, prepare, track
from swmtrack.preprocessing import PreprocessConfig, preprocess
from swmtrack.workflow import run_pipeline

config = PreprocessConfig.recommended()
```

[API examples](docs/python_interface.rst) cover end-to-end and prepared-surface
usage. Dataclass fields, defaults and ordinary call signatures are preserved.
Core source inspection and Numba Dispatcher internals are outside the compiled
interface.

## Validation and distribution

The installed compiled wheel passed **155 portable tests**. All **107 checked
API signatures/dataclass definitions**, 5,000 RNG draws, **70 prepared tracking
runs** and **four fresh FreeSurfer/FOD cases** matched the Python reference.
Observed curvature-feature and streamline-coordinate differences were zero in
those comparisons. The synthetic example produces **98 accepted tracks from
256 attempts**. See [validation scope](docs/validation.rst).

[Public CI](.github/workflows/ci.yml) audits repository contents, wheel integrity,
installation and the synthetic example. The full scientific tests and core
implementation are maintained privately. Human datasets and cohort outputs are
excluded from downloads. Numerical agreement does not establish anatomical
accuracy, complete sulcal coverage or correct registration.

Current `main` contains documentation, approved general helpers and compiled
artifacts. Earlier public commits and the preserved historical branches remain
accessible. Compilation excludes core source files from new downloads; it does
not prevent reverse engineering. See [distribution scope](docs/downloads.rst).

## Cite both methods

Please cite **both** articles when using this pipeline. If you prepare the ROIs and WM surface yourself, you can cite only 2: 

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
describes dependency terms and redistribution notices. Project-wide licensing
terms remain unspecified; these dependency notices do not license SurfTracker
itself.
