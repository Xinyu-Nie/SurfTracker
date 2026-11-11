Changes
========

0.3.0.dev16
-----------

* Prepare the GitHub Python release tree and a new main README with installation,
  inputs, between-gyri, merged-region and whole-hemisphere workflows.
* Preserve the legacy C++ directory and notices unchanged; remove MATLAB from
  the current tree while retaining repository history.
* Exclude legacy source and the existing example archive from Python packages.
* Configure PR/main CI checks without enabling package publication.
* Make rotation/native-interface fixtures portable across current SciPy and
  operating systems; cap requested tracking workers to the available Numba pool
  and record the effective count without changing seeded trajectories.

0.3.0.dev15
-----------

* Rename the PyPI distribution to SurfTracker, retaining the swmtrack import
  package, python -m swmtrack and swmtrack console command.
* Update archive checks, metadata fixtures and publication/install examples.
* Verify installed distribution metadata and the registered console entry point
  outside the checkout. Retain scientific algorithms and default settings.

0.3.0.dev14
-----------

* Prepare build/artifact CI for Linux/Python 3.10--3.13 and macOS/Windows 3.12.
* Add a self-contained analytic prepare/track/export example with QC VTKs.
* Document HCP S900 redistribution conditions and a separate data route.
* Document TestPyPI and Trusted Publishing preparation without enabling uploads.
* Keep copyright holder unspecified and scientific defaults unchanged.

0.3.0.dev13
-----------

* Add two verified primary paper citations in CFF, BibTeX and the README.
* Add human software-author/contributor roles, a separate AI-assistance
  acknowledgment, project URLs and scoped dependency notices, with research
  affiliations separate from copyright ownership.
* Prepare an unmodified BSD-3-Clause proposal without asserting an unconfirmed
  copyright holder or activating a license grant.
* Preserve the runtime algorithms/defaults; only the version string changes.

0.3.0.dev12
-----------

* Prepare a clean Python release tree with a portable synthetic test suite.
* Preserve the recommended-default regression checks without importing
  historical experiment scripts; retain those experiments in the original
  research workspace.
* Explicitly include the portable tests, documentation and helper scripts in
  the source archive. Add a source/archive content and version checker.
* Add a portable white/inflated QC-label exporter for completed pipeline runs.
* Retain the dev11 numerical algorithms and defaults; only the package version
  string changes in the runtime source modules.

0.3.0.dev11
-----------

* Add optional whole-hemisphere sulcal seeding with a merged gyral target,
  physical length limits and original endpoint-parcel classification.
* Export endpoint/length cell scalars, per-track tables and parcel-pair counts.
* Batch VTK connectivity writes while preserving points and streamline order.

Earlier Python development
--------------------------

* Recommended implicit-grid curvature with consistent FEM mass smoothing,
  robust-unit segmentation energy and existing Python HJ/pruning/cleanup.
* Connectivity-filtered intersection seeds between gyri; full sulcal
  intersection within an explicitly merged region.
* ROI-name/packed-ID lists, optional ANTs rigid registration and hash-verified
  checkpoints from FreeSurfer and supplied FODs to validated tractograms.
