Binary downloads
================

Current version: 0.3.0.dev18. Install the supplied wheel in a fresh CPython 3.12
Linux x86_64 environment. The wheel is tagged manylinux2014 / glibc 2.17; the
executed scientific validation used a glibc 2.31 host. Other Python versions,
architectures and operating systems require their own tested binaries.

The binaries/0.3.0.dev18/ directory contains the wheel, complete binary kit,
validated dependency constraints, manifest.json and SHA256SUMS. The kit supplies
the same wheel plus examples, documentation and notices. After downloading,
verify checksums with sha256sum -c SHA256SUMS from the matching directory.
The wheel hash is aad58cb904e7bdb253dd0443391c1020459f84cdb57c4820f9fdafad0e47d3da.

Installation
------------

Unzip the complete kit, activate a fresh Python 3.12 environment and run::

    python -m pip install -c requirements-validation.txt *.whl
    swmtrack doctor
    swmtrack run --help
    python scripts/synthetic_demo.py --output-dir output/synthetic

The constraints reproduce the NumPy/SciPy/Numba/NiBabel/meshio versions used
for the scientific comparisons. Dependencies are installed from PyPI. No
compiler or Cython runtime is required. Numba remains a dependency for thread
configuration compatibility; the core performs no JIT compilation.

The distribution name is SurfTracker, and both the import package and console
command are swmtrack. Use the compiled download for this version. The public
repository contains documentation and supporting helpers rather than a source
package to build with pip install . . Earlier package-index artifacts are
separate from this compiled download; this update makes no PyPI/TestPyPI upload.

Public and private scope
------------------------

The wheel contains ten native scientific modules and readable CLI, I/O,
FreeSurfer import, registration, validation, metadata/API adapters and notices.
Protected scientific Python sources, generated C/Cython files, caches, build
logs and human datasets are excluded. Dataclasses and ordinary Python call
signatures are preserved; Python source inspection and Numba Dispatcher
internals are outside the compiled interface. See python_interface.rst.

Core implementation and scientific tests are maintained privately. This public
main branch contains documentation, examples, general QC/audit helpers and
compiled downloads. Historical public commits and the legacy and dev-c++
branches remain accessible. Changing the current tree does not withdraw earlier
copies, and compilation does not prevent reverse engineering.

Contributor credits and citations are preserved. Project-wide licensing terms
remain unspecified; the dependency notices apply to their respective components
and do not constitute a license grant for SurfTracker itself.

Updates and caches
------------------

Record the installed version, dependency versions, processing parameters and
artifact checksum with a study. The workflow records installed implementation
hashes in new caches and checkpoints. Changed native modules or settings require
a new output directory and recomputation. Older portable NPZ caches without a
build identity remain readable under their existing algorithm-version checks;
recompute them for a new production run. Do not mix the old swmtrack distribution
and SurfTracker in one environment.
