QC outputs
==========

The pipeline writes the full preprocessed white VTK and the inward tracking
surface, plus tracks.vtk, tracks.tck, tracks.trk and numerical validation reports.
The portable source helper scripts/export_qc.py exports label-only white and
optional native inflated surfaces for every requested completed hemisphere.
The helper never reruns or modifies segmentation, registration or tracking.

Run the helper from a source checkout or extracted source archive with an
installed swmtrack package::

    python scripts/export_qc.py --run-dir output/whole \
        --output-dir qc/whole --inflated

A new QC directory contains HEMI/labels.white.vtk, HEMI/labels.inflated.vtk when
requested, and qc.json with source/output hashes and exact roundtrip checks.
Without --inflated, FreeSurfer geometry is not read again. With --inflated, the
recorded subject_dir locates surf/HEMI.inflated; --subject-dir overrides it if
inputs have moved. Vertex/face identity must match the white surface.

Point fields
------------

* candidate_segmentation: 0 gyrus, 1 sulcus, 2 outside tracking domain.
* seed: eligible seed vertices; roi: full propagation domain.
* endpoint_region: 0 outside target, 1 A or merged target, 2 B.
* anatomical_roi: original annotation IDs, retained alongside the original
  annotation array, source_vertex_id and numerical preprocessing fields.
* RGB arrays for segmentation, anatomical parcels and seed eligibility.

Whole-hemisphere tract cell fields
----------------------------------

endpoint_a_id and endpoint_b_id correspond to the first/last trajectory points.
pair_low_id and pair_high_id identify an unordered pair. endpoint_class is
1 same parcel, 2 cross parcel, 3 ambiguous. ID 0 means the gyral vertices of the
terminal face disagree on parcel ID; no arbitrary nearest parcel is substituted.
length_mm is integrated physical length; accepted_trial_id links to diagnostics.
CSV/NPZ tables have the same order as the full tractogram. TCK/TRK carry geometry
and use the accompanying tables for parcel metadata.

Tighter length/endpoint selection can be applied to retained tracks through
these cell scalars or tables. A wider range cannot recover tracks excluded by
the original propagation/length limit. Named subsets can leave/re-enter their
selected ROI because the underlying domain was the whole hemisphere.

Viewing
-------

Tracks use the saved 0.5-mm inward surface. Overlay physical tractograms with
projection/tracking_surface.world.vtk or a transparent white surface; an opaque
white surface can hide inward tracks. Native inflated coordinates are only for
label QC and do not align with the physical tractogram. Inspect segmentation,
seeds, parcels and endpoints, together with actual T1/diffusion overlays.
Numerical checks and a nonempty bundle do not establish anatomical accuracy,
complete sulcal coverage or correct registration.
