SurfTracker: compiled superficial white matter tracking
=======================================================

This repository provides documentation, examples and compiled binary downloads.
The current build is 0.3.0.dev18 for Linux x86_64 / CPython 3.12. The Python import
and console command remain swmtrack. The scientific core is maintained privately.

See README.md for installation, MRtrix3 FOD compatibility, physical RAS coordinate
handling, between-gyri, merged-region and whole-hemisphere examples, optional
ANTs registration and QC. Download the complete binary kit from
binaries/0.3.0.dev18/ and follow its README. This repository contains no source
package to build with pip install . .

Documentation
-------------

* docs/downloads.rst: supported binary, checksums and distribution scope.
* docs/end_to_end.rst: input requirements, commands and recommended defaults.
* docs/input_conventions.rst: MRtrix3 SH coefficients and coordinate frames.
* docs/python_interface.rst: preserved Python calls.
* docs/qc.rst: segmentation, seed/ROI and tractogram QC.
* docs/validation.rst: compiled/source agreement and scientific limits.
* docs/example_data.rst: synthetic example and separate HCP data terms.

Citations
---------

Please cite both primary methods for use of the full pipeline:

* Shi et al. Hamilton-Jacobi Skeleton on Cortical Surfaces (2008).
  https://doi.org/10.1109/TMI.2007.913279 .
* Nie et al. Surface-Based Probabilistic Fiber Tracking in Superficial White
  Matter (2024). https://doi.org/10.1109/TMI.2023.3329451 .

CITATION.cff and CITATIONS.bib provide complete credits. AUTHORS.rst records
contributions. THIRD_PARTY_NOTICES.rst records dependency notices. No project-wide
license grant or copyright-holder assignment is added by this update.
