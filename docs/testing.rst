Checking the public downloads
=============================

The public repository supplies installation/QC/audit helpers and a synthetic
example. Scientific implementation tests are maintained privately. To inspect
a Git checkout of this documentation/download tree without importing the core::

    python scripts/check_public_repository.py

This verifies approved files, local Markdown links, download hashes, wheel
payload/RECORD contents and the binary-kit boundary. To inspect only a wheel::

    python scripts/check_binary_wheel.py binaries/0.3.0.dev18/*.whl

In an extracted binary kit, the wheel is in its root directory rather than
binaries/0.3.0.dev18/. Use SHA256SUMS and check_binary_wheel.py for its integrity.

After installing the wheel in a fresh Python 3.12 Linux environment, run::

    python -m pip check
    swmtrack doctor
    python -m swmtrack run --help
    python scripts/synthetic_demo.py --output-dir output/synthetic
    python scripts/validate_installation.py --demo output/synthetic/demo.json

The synthetic case exercises the installed native package, endpoint constraints
and independently checked TCK/TRK/VTK exports. It contains no human data. The
validator confirms native module origins, interfaces/defaults and the expected
98 accepted streamlines/digest. These are download checks rather than cortical
segmentation or anatomical validation. Full compiled/source comparison scope is
reported in validation.rst. No check here builds source distributions or uploads
a package to a package index.
