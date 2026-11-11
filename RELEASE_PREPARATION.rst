Release preparation
===================

The Python package in this repository is prepared for SurfTracker (PyPI
distribution: SurfTracker; Python import and command: swmtrack). The main branch
contains the Python pipeline, portable tests, public workflow documentation and
four generic helper scripts. Historical implementations, their original notices
and the old example archive are retained on the legacy and dev-c++ branches and
in Git history. Subject-specific experiments remain in the research workspace.

The main tree and Python distributions contain no historical source/executable/
library payload, cohort dataset, large result directory, virtual environment or
research-test dependency. The old 101107.tar.gz example pointer and its sole
Git LFS rule have been removed from main. MANIFEST.in continues to guard against
accidentally including historical source/data in Python distributions.

The development version is 0.3.0.dev17. The optional external native-executable
backend and --tools-json interface have been removed. The independent Python
algorithms, FOD convention identifier and recommended numerical defaults remain
unchanged. Portable unit/integration tests and clean archive/install checks are
rerun. The copyright holder remains unspecified and no active project license
has been added.

Before publication
------------------

* Confirm the project copyright holder/license authorization and review the
  prepared author/maintainer metadata before adding an active LICENSE.
  BSD-3-Clause is proposed, with the holder still unresolved; this preparation
  does not activate a grant or assign rights to historical third-party code. scripts/check_distributions.py --require-license will fail
  until license text and distribution metadata are provided.
* Review the prepared citation/contributor metadata and PyPI naming, the scientific QC and
  supported platform/Python range, and configure repository CI/publishing.
* Run python -m twine check and the isolated wheel/sdist test procedure against
  the exact final artifacts; rehearse on TestPyPI before a public release.
* Preserve the existing GitHub history and scope the tagged source tree to
  the Python release. Repository preparation does not publish a PyPI release.

Additional preparation in dev14
-------------------------------

Prepared build/artifact CI, an analytic synthetic tracking example, HCP S900
data-sharing guidance and a TestPyPI/Trusted Publishing runbook. The holder
remains unspecified; no active license or publisher workflow was added. The CI
matrix runs on PRs/main; PyPI account/contact setup remains outstanding,
and a separate human-data example needs confirmed source-release provenance.
See docs/example_data.rst and docs/publishing.rst.

Distribution naming in dev15
----------------------------

The PyPI distribution is SurfTracker. Python imports, module execution and
the console command remain swmtrack. Build artifacts use the normalized
surftracker name; publishing examples and archive verification follow that
name. Both installed distribution metadata and the swmtrack entry point are
checked outside the source checkout. No package upload is performed.

GitHub preparation in dev16
---------------------------

README.md now describes the Python package and links to the preserved C++
implementation. Legacy C++ source and notices remain byte-identical to the
base revision. MATLAB is removed only in the proposed current tree; history
is preserved. Explicit manifest exclusions keep legacy source and the existing
Git LFS example pointer out of Python distributions. CI checks pull requests
and main without enabling automatic package uploads.

Python-only main in dev17
-------------------------

Remove Surface_tracking_c++, 101107.tar.gz and .gitattributes from main while
preserving the legacy and dev-c++ branch tips. The attribute file contained only
an LFS rule for *.gz, and the deleted example was its only tracked match.
Remove the external native backend, executable overrides and payload tests;
update public documentation to link to the historical branches. Retain current
Python computations and the supported FOD convention.
