Python interface
================

Install the compiled wheel first. Existing swmtrack module paths and ordinary
call signatures are retained. For example::

    from swmtrack.tracking import (
        TrackingConfig, Prepared, prepare, track, save_prepared, load_prepared)
    from swmtrack.preprocessing import PreprocessConfig, preprocess
    from swmtrack.workflow import run_pipeline

End-to-end targeted tracking
----------------------------

Use an explicit anatomical alignment choice and SH convention::

    report = run_pipeline(
        subject_dir="/data/subject/freesurfer",
        fod="/data/subject/wmfod.nii.gz",
        output_dir="output/central",
        hemispheres=("lh",),
        sh_basis="legacy-real-even",
        assume_scanner_aligned=True,
        gyrus_a=["precentral"],
        gyrus_b=["postcentral"],
    )

Within-region tracking accepts a list of original ROI names or packed annotation
IDs. The pipeline constructs the union; a new annotation label is not needed::

    report = run_pipeline(
        subject_dir="/data/subject/freesurfer",
        fod="/data/subject/wmfod.nii.gz",
        output_dir="output/parietal",
        hemispheres=("lh",),
        sh_basis="legacy-real-even",
        assume_scanner_aligned=True,
        gyrus_a=["superiorparietal", "inferiorparietal", "postcentral",
                 "supramarginal", "precuneus"],
        tracking_config=TrackingConfig(endpoint_policy="within_gyrus"),
    )

Omit gyrus_b in within-region mode. Integer ID lists and their string forms are
accepted in place of names; use original annotation color-table IDs rather than
volume LUT IDs. Whole-hemisphere tracking uses whole_hemisphere=True and omits
ROI selection. See end_to_end.rst for endpoint semantics and classification.

Prepared surfaces
-----------------

A custom VTK can use explicit seed, roi, gyrus_a/gyrus_b or one binary gyrus point
mask. Its points must already lie in FOD-world RAS millimeters::

    prepared = prepare("surface.vtk", "wmfod.nii.gz")
    save_prepared("prepared.npz", prepared)
    prepared = load_prepared("prepared.npz")
    streamlines, tracking_report, diagnostics = track(prepared, TrackingConfig())

For one merged endpoint target, supply the binary gyrus point array and select
TrackingConfig(endpoint_policy="within_gyrus"). Seed/target arrays must follow
the desired sulcal/gyral selection; prepare does not infer missing anatomy.

Configurations and compatibility
--------------------------------

PreprocessConfig.recommended() supplies the validated curvature/segmentation/HJ
settings. No-config preprocess and run_pipeline use that profile. The bare
PreprocessConfig() retains the earlier explicit cotan baseline for compatibility.
TrackingConfig() retains area-weighted seeding, RNG 42 and the documented physical
integration/length defaults. See end_to_end.rst for all recommended settings.

The configuration and data containers remain ordinary dataclasses, including
runtime field types, defaults and array attributes. Core functions load from
native modules. Their source bodies and Numba Dispatcher-specific attributes
are not distributed. Projected caches and resume checkpoints include the
installed implementation identity; see downloads.rst before reusing old caches.
