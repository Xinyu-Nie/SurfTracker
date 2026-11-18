Validation of the compiled build
================================

Version 0.3.0.dev18 was built and tested on Linux x86_64 / CPython 3.12 with
NumPy 2.3.3, SciPy 1.16.2, Numba 0.68.1, llvmlite 0.50.0, NiBabel 5.4.2 and
meshio 5.3.5. The wheel bytes match the privately tested artifact; see
binaries/0.3.0.dev18/manifest.json and SHA256SUMS. Core implementation commit:
96895c083f26b141aae29fee75188e6429052173 (private source provenance).

Installed package and API
-------------------------

* 155 portable tests pass against the installed native wheel outside its source
  checkout. These cover geometry, SH projection, tracking, preprocessing, ROI
  seeds, endpoint rules, I/O, registration contracts and workflow/cache behavior.
* All 107 checked source API signatures/dataclass definitions match, including
  resolved runtime field types and defaults.
* A clean environment without Cython passes the actual swmtrack command,
  python -m swmtrack, imports, pip check and the synthetic example.
* The synthetic example produces 98 accepted streamlines from 256 attempts.
  Its streamline digest is
  8c4b8046e76d972584759a17a05cdef0f5d929214c2a538b783c4253103266bb.
* Wheel/RECORD/native-library audits and README rendering pass. Eleven
  invalid source/cache/data/tamper payload variants are rejected. The private
  build rejects source-archive creation.

Scientific implementation agreement
------------------------------------

Five thousand SplitMix64 draws match exactly, including unsigned overflow cases.
Thirty-five prepared cases cover the four-subject central, temporal, merged
whole-parietal and whole-hemisphere cohort plus existing native MRtrix3 raw and
normalized FOD conditions. Each uses 512 attempts with one and four workers and
different batch sizes: 70 runs and 35,840 attempts per implementation. Seed
arrays, accepted IDs, stops, endpoints, lengths and coordinates match exactly.

Four fresh FreeSurfer/FOD cases each use 512 attempts:

===========================  ===============  ===========================
Mode                         Accepted tracks  Maximum observed difference
===========================  ===============  ===========================
Precentral/postcentral LH    131              0 for features/coordinates
Superior/middle temporal RH  76               0 for features/coordinates
Merged whole-parietal LH     163              0 for features/coordinates
Whole hemisphere RH          195              0 for features/coordinates
===========================  ===============  ===========================

Initial/cleaned/final labels, skeletons, masks, geometry, projection and seed
arrays match exactly. Both implementations pass verified workflow resume and
tractogram export checks. These comparisons used the recommended settings with
the attempt count fixed for validation. Scientific defaults are unchanged.

Limits
------

These checks establish compiled/source agreement on the tested inputs, not
anatomical correctness. Prior cohort/sensitivity work showed subject-dependent
temporal coverage and responses to inward offset and FOD threshold. Review
segmentation, seed coverage, endpoints and anatomical alignment in each study.

Only the stated Linux/Python/dependency combination was executed for the
compiled build. Native dependency auditing permits manylinux2014 / glibc 2.17
tags; execution used glibc 2.31. Other operating systems, architectures, Python
ABIs and dependency combinations need validation. Optional isolated ANTs worker
contracts remain covered; new ANTs fits were not part of these comparisons.
No human input datasets or cohort tractograms are included in the downloads.

Public download checks
-----------------------

Public CI audits the documentation/download tree, checks artifact hashes and
native payloads, installs the wheel and runs the standalone synthetic example
and API/default checks. It does not rebuild the private core or claim to rerun
the full scientific suite. See testing.rst and the repository Actions results.
