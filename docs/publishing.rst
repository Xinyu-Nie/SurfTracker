Release rehearsal and publication
=================================

Current status
--------------

The PyPI distribution and repository/project name are SurfTracker. The Python
import package and console command remain swmtrack. This preparation does not
reserve a PyPI name, create accounts,
configure repository environments, push commits or upload anything. The legal
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
read-only repository permissions and has no publishing step. Its intended
matrix is Linux/Python 3.10--3.13 plus macOS/Windows on Python 3.12. These are
planned CI targets, not completed platform validation: only locally executed
results should be claimed until the workflow actually passes on GitHub.

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
permissions separate from pull-request jobs. No publisher workflow is enabled
here while release ownership/license and account configuration are unresolved.

References
----------

* https://docs.github.com/en/actions/tutorials/build-and-test-code/python
* https://packaging.python.org/en/latest/guides/using-testpypi/
* https://docs.pypi.org/trusted-publishers/using-a-publisher/
