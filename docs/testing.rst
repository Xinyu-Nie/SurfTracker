Portable tests and artifact verification
========================================

The default tests/ suite contains runtime, geometry, SH projection, segmentation,
seed-policy, registration-interface, checkpoint, endpoint, I/O and helper tests.
Fixtures are constructed locally from synthetic meshes/images/annotations.
Two legacy-payload checks skip when the separately supplied optional payload is
absent. Default Python operation does not require that payload.

Historical algorithm experiments and their tests remain in the original
research workspace. They are not silently skipped by an import guard or shipped
as broken tests: they are explicitly outside the clean release tree. The
recommended-profile tests use analytic FEM and Gaussian-likelihood expectations
and check CLI overrides directly without importing experimental drivers.

From a clean checkout::

    python -m pip install '.[test,release]'
    OPENBLAS_NUM_THREADS=1 python -m pytest -q
    python -m build
    python scripts/check_distributions.py dist/*.whl dist/*.tar.gz
    python -m twine check dist/*

The source archive explicitly contains its portable test helpers and scripts.
Extract it into a new directory and run pytest there without adding the original
research workspace to PYTHONPATH. Also rebuild a wheel from that extracted
source; archive checks should report matching package source and version.

To test the installed wheel independently, install the runtime dependencies and
pytest in a fresh environment, install that wheel with --no-deps, and invoke
pytest against the extracted archive's tests/ from an unrelated empty directory.
Use --import-mode=prepend (pytest's default) so synthetic test-fixture modules
can import one another. Verify that swmtrack.__file__ points into the installed
environment, not the original checkout. Do not add the original scripts/ or
research directory to PYTHONPATH. The generic script tests invoke the archive's
helper files in temporary working directories using the package under test.

No PyPI/GitHub upload is part of these checks. Licensing, README rendering and
supported-platform CI execution are separate final release requirements. The content
checker accepts the current license-pending development archives, reports
project_license_ready=false, and rejects them under --require-license.

Automated offline artifact check
--------------------------------

With runtime/test/build dependencies already installed::

    python scripts/verify_artifacts.py \
        --wheel dist/surftracker-0.3.0.dev16-py3-none-any.whl \
        --sdist dist/surftracker-0.3.0.dev16.tar.gz \
        --output-dir release-audit/artifact-tests

This extracts the archive safely to a fresh temporary directory, runs its tests,
installs the wheel into a separate target, verifies the imported package path,
and repeats the tests from an unrelated working directory. It also verifies SurfTracker distribution metadata, exercises the registered
swmtrack console entry point and module CLI, and rebuilds the wheel from the
extracted source. No original
workspace or research-script path is added. Dependency resolution and platform
coverage require separate fresh-environment CI; this check reuses the current
interpreter's numerical/test/build dependencies and makes no network requests.

Prepared continuous integration
--------------------------------

The GitHub Actions workflow in .github/workflows/ci.yml runs the above artifact
checks and strict Twine checks on each matrix environment, plus the synthetic
tracking example. It does not upload to PyPI/TestPyPI. The intended matrix
requires actual remote execution before claiming additional tested platforms.
See publishing.rst for the publishing rehearsal and account setup.
