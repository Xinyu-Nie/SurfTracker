Release preparation
===================

The Python package in this repository is prepared for SurfTracker (PyPI
distribution: SurfTracker; Python import and command: swmtrack). The historical C++ source is retained separately in Surface_tracking_c++/;
MATLAB is removed from the current tree and remains in Git history. Runtime modules, portable tests, public workflow documentation and
four generic helper scripts are included. Subject-specific benchmark, license
recovery, external-tool bundling, dataset reconstruction and cohort plotting
scripts remain in the original research workspace.

The Python distributions contain no legacy source/executable/library payload,
cohort dataset, large scientific result directory, virtual environment, research
test dependency or original Git metadata. MANIFEST.in controls their contents.
A GitHub source archive instead contains its tagged tree, including the
preserved C++ source and the existing example-data pointer.

The development version is 0.3.0.dev16. Dev12 separated portable tests; dev13
adds verified citations, human contributor roles and a scoped BSD-3-Clause
proposal. The frozen dev11 numerical validation artifacts remain intact. Numerical algorithms and defaults are unchanged. The runtime version string
and optional native-payload reinstall message reflect the new distribution. No cohort computation
needs to be repeated to validate this test/documentation reorganization.
Portable unit/integration tests and clean archive/install checks are rerun.

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
* Preserve the existing GitHub history and scope the new tagged source tree
  to the intended Python release. This preparation does not initialize Git,
  change the remote repository, upload packages or publish a release.

Additional preparation in dev14
-------------------------------

Prepared build/artifact CI, an analytic synthetic tracking example, HCP S900
data-sharing guidance and a TestPyPI/Trusted Publishing runbook. The holder
remains unspecified; no active license or publisher workflow was added. The CI
matrix needs remote execution, account/name/contact setup remains outstanding,
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
