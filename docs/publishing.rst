Release rehearsal and publication
=================================

Current status
--------------

The PyPI distribution and repository/project name are SurfTracker. The Python
import package and console command remain swmtrack. The manual TestPyPI workflow
is .github/workflows/publish-testpypi.yml, with GitHub environment testpypi.
Adding the pending publisher on TestPyPI is a separate account-configuration
step. No production PyPI publisher is configured by this setup. The legal
copyright holder remains unspecified. BSD-3-Clause remains a proposal.

The --require-license archive check is a project release-readiness choice,
not a PyPI requirement to prove copyright ownership. Default CI deliberately
accepts license-pending development artifacts and reports that status. Complete
the project's license/authorization review before a public licensed release.

Automated checks
----------------

.github/workflows/ci.yml builds wheel/sdist, validates archive contents and
strict README rendering, tests extracted source and installed wheels outside
the checkout, rebuilds from the sdist, and runs the synthetic example. It uses
read-only repository permissions and has no publishing step. Its configured
matrix is Linux/Python 3.10--3.13 plus macOS/Windows on Python 3.12. Inspect
the current PR checks or repository Actions page for the exact run results.
These synthetic/package checks do not replace scientific cohort validation.

From a fresh environment::

    python -m pip install '.[test,release]'
    python -m build
    python -m twine check --strict dist/*
    python scripts/check_distributions.py dist/*.whl dist/*.tar.gz
    python scripts/verify_artifacts.py --wheel dist/*.whl --sdist dist/*.tar.gz \
        --output-dir release-audit/artifact-tests

For a final licensed release, additionally use --require-license. Build in a
clean directory so dist contains exactly one wheel and one source archive.
Use a new version for each changed public artifact; existing PyPI files cannot
be replaced. Preserve the exact checked artifacts and their SHA256 checksums.

TestPyPI rehearsal
------------------

The distribution name is SurfTracker (normalized to surftracker). Confirm
the maintainer email and the actual PyPI/TestPyPI account administrators.
Repository ownership is not software copyright
ownership. TestPyPI and PyPI have separate accounts and publisher settings.
An available name is not a reservation until the service accepts the project.

GitHub Trusted Publishing setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a GitHub environment named testpypi under repository Settings ->
Environments, restricted to deployments from the main branch. It needs no
package-index password or API-token secret. In the GitHub tab of the TestPyPI
account's Publishing page, add this pending publisher:

===================== ========================
Field                 Value
===================== ========================
TestPyPI Project Name SurfTracker
Owner                 Xinyu-Nie
Repository name       SurfTracker
Workflow name         publish-testpypi.yml
Environment name      testpypi
===================== ========================

Owner is the GitHub username, not the TestPyPI account username. The workflow
field is the filename, without .github/workflows/. All values must match the
workflow/repository exactly. Adding a pending publisher does not upload a
package or reserve its name. The first successful upload creates the project.

Manual validation and upload
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wait for Package checks to pass for the main commit to publish. In repository
Actions, select Publish to TestPyPI -> Run workflow -> main. Leave the
"Upload the verified distributions to TestPyPI" checkbox unchecked for a
validation-only rehearsal. The workflow finds a successful push-triggered
ci.yml run for that exact commit and retrieves its Linux/Python 3.12 artifacts.
It verifies their recorded SHA256 hashes, package-source hashes, archive
contents and README rendering. It does not rebuild distributions. Missing or
expired CI artifacts require rerunning Package checks first.

After the pending publisher is configured and distribution is authorized, run
the workflow again on main with the upload checkbox checked. Only this explicit
choice enables the publishing job, which uses the testpypi environment and
job-scoped id-token: write. The validated archives are passed to the publishing
job as an artifact and uploaded only to https://test.pypi.org/legacy/.
The workflow has no push, tag or pull-request publishing trigger and no
production PyPI endpoint. Existing filenames are not silently skipped: use a
new package version for changed distributions after an upload.

The current package version remains the version in pyproject.toml; this setup
does not automatically change it to a release candidate or activate a license.
TestPyPI is publicly accessible, so it is a distribution rehearsal rather than
private storage. After upload, install and test from the index as shown below.

Alternative manual upload
~~~~~~~~~~~~~~~~~~~~~~~~~

After the release is authorized, upload only the checked artifacts, for example::

    python -m twine upload --repository testpypi dist/*

This is a manual publishing instruction, not executed by preparation scripts.
Use a project-scoped credential or a configured Trusted Publisher; do not put
tokens in source files. Then install from TestPyPI in a fresh environment.
Install known runtime/test dependencies from PyPI first, and obtain just the
project from TestPyPI with --no-deps, avoiding mixed-index dependency selection::

    python -m pip install 'numpy>=1.24' 'scipy>=1.11' 'nibabel>=5' 'numba>=0.61' 'meshio>=5' pytest
    python -m pip install --no-deps --index-url https://test.pypi.org/simple/ SurfTracker==VERSION
    swmtrack doctor

Replace VERSION with the exact rehearsal version. When migrating from the
older distribution named swmtrack, use a fresh environment: both distributions
install the same swmtrack import package and command. Do not install both
distributions together. Query installed version metadata with
importlib.metadata.version("SurfTracker"); import swmtrack remains unchanged.

Compare downloaded archive hashes with those tested locally and repeat
installed-artifact checks.

Tagged production release
-------------------------

Preserve the existing GitHub history when updating the repository to the
curated Python tree. Review the final diff and ensure that institutional
datasets, license-review drafts and research outputs are absent. After CI,
scientific review and license readiness pass, tag the reviewed commit with a
version matching pyproject.toml, swmtrack/__init__.py and CITATION.cff.

Configure a PyPI Trusted Publisher with the exact GitHub owner, repository,
workflow filename and protected environment. Restrict publication to reviewed
version tags and require the configured environment's approval. A future
publisher job should have id-token: write, consume the verified build artifacts
without rebuilding them, and use PyPA's publishing action. Keep publishing
permissions separate from pull-request jobs. The manual TestPyPI workflow is
separate from this production-release setup; no production publisher is enabled
here while release ownership/license and account configuration are unresolved.

References
----------

* https://docs.github.com/en/actions/tutorials/build-and-test-code/python
* https://packaging.python.org/en/latest/guides/using-testpypi/
* https://docs.pypi.org/trusted-publishers/using-a-publisher/
