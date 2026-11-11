Validation scope
================

The current numerical baseline was evaluated with the dev10/dev11 artifacts
before this release-tree cleanup. Scientific outputs and dataset files remain
in the original research workspace and are not redistributed here. Subsequent
release-preparation revisions preserve the numerical algorithms and defaults.
Distribution metadata, version strings and installation messages may differ.

Fresh FreeSurfer + FOD targeted tests
-------------------------------------

Four subjects (100307, 100408, 101309, 102008), both hemispheres, and central,
superior/middle temporal and merged whole-parietal tracking were evaluated in
24 independent public-CLI runs using recommended defaults. No segmentation or
projection cache was supplied. All 24 runs were nonempty and passed numerical
export/endpoint checks. Accepted totals were 185,445 central, 50,692 temporal
and 802,044 whole-parietal trajectories. These runs reproduced the frozen
recommended-configuration trajectory coordinates exactly.

Whole-hemisphere tests
----------------------

The same four subjects and both hemispheres were evaluated in eight fresh
public-CLI runs. Sulcal cortical vertices supplied seeds; all gyral cortical
regions formed one binary termination target. All eight outputs were nonempty,
with 3,860,835 accepted trajectories in total. Every saved curvature, initial/
cleaned/final label, depth, flux and skeleton array exactly matched the
corresponding preceding targeted run.

Endpoint-classified subsets contained 186,772 central, 60,599 superior/middle
temporal and 829,985 whole-parietal trajectories. The full tractograms retained
1,507,858 same-parcel, 2,350,170 cross-parcel and 2,807 ambiguous-parcel paths.
Ambiguous endpoints were excluded from named subsets. Named subsets are based
on endpoints and may leave/re-enter their selected ROI; they do not reproduce
the restricted-domain targeted experiment exactly.

Checks and limits
-----------------

Checks include independent TCK/TRK readback, finite/in-image coordinates for all
points, physical-surface membership sampling, every endpoint against its recorded
triangle and required target, integrated length limits, and exact VTK point/
polyline/cell-scalar readback. QC label surfaces retain original vertex identities
and point arrays. The default tracking range is 8--150 mm.

Existing HCP scanner/FOD alignment was explicitly declared; these cohorts did
not fit an additional ANTs transform. Optional registration has separate
synthetic affine/landmark and workflow tests. Coordinate bounds and successful
optimization do not establish anatomical registration.

Both eight-hemisphere contact sheets were inspected, with nonempty central,
temporal and whole-parietal subsets. Projected views sampled at most 1,000 paths
per selection; full VTKs retain all accepted paths. Tract counts, spatial overlap
proxies and numerical checks do not establish anatomical correctness, complete
sulcal coverage, biological sensitivity or accuracy of individual U-fibers.
Expert review and appropriately scoped claims remain necessary for interpretation.

Portable packaging tests
-------------------------

The release tree's tests generate synthetic inputs. Research experiment tests
remain in the original workspace. Source, extracted-sdist and installed-wheel
tests run independently without putting that workspace on PYTHONPATH. The
release-audit directory in the prepared local tree records the exact test counts,
archive/source hashes and artifact checks; it is excluded from source archives
and future Git tracking. See testing.rst for the reproduction procedure.

Small tracking sensitivity screen
---------------------------------

A dev14 screen used subjects 100307 and 100408, left hemisphere only, for
central and superior/middle temporal between-gyri tracking. Each case used
2,000 attempts, not the default 20 attempts/mm2. Saved recommended segmentation
and ROI/seed/endpoint masks were fixed. Seven conditions per selection tested
RNG 42/43/44, inward offset 0.4/0.5/0.6 mm, and minimum FOD 0.008/0.01/0.012,
changing one setting at a time. Nonbaseline offsets recomputed FOD projection
on the shifted surface. All 28 conditions were nonempty and passed physical
length, endpoint-target, image-bound and independent TCK/TRK readback checks.

Accepted counts (each condition has 2,000 attempts)::

    Subject  Selection  Baseline  RNG 43/44  Offset .4/.6  FOD .008/.012
    100307   central         523  562/554    528/560        575/488
    100307   temporal        258  242/245    211/256        321/202
    100408   central         581  565/567    576/608        642/537
    100408   temporal        252  278/272    263/267        272/243

For 100307 temporal tracking, lowering/raising the FOD threshold changed
accepted counts by approximately +24%/-22%, while reducing the inward offset
to 0.4 mm changed counts by about -18%. RNG changes alone produced smaller
changes there; 100408 showed a different response. This supports documenting
threshold/offset sensitivity and reviewing temporal coverage. It does not
establish parameter robustness or justify retuning defaults to maximize counts.

The recommended defaults remain 0.5-mm offset, minimum FOD 0.01 and RNG 42.
This fixed-count screen supplements the preceding full-density cohort runs; it
does not replace all-hemisphere review, test registration sensitivity, or prove
anatomical correctness. Full scientific outputs stay outside the package. HCP
data acknowledgment and data-use conditions are in example_data.rst.
