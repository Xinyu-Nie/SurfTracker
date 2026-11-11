Workflow examples
=================

Install this source tree with ``python -m pip install .``. After publication,
install the PyPI distribution ``SurfTracker``. The Python import and command
remain ``swmtrack``. Run ``swmtrack run --help`` for every input/override.
Replace /data paths with your
own FreeSurfer and supported FOD inputs. Explicit alignment and SH convention
are required; no institutional dataset or cohort script is needed.

Central::

    swmtrack run --subject-dir /data/subject/freesurfer --fod /data/subject/fod.nii.gz \
        --output-dir output/central --assume-scanner-aligned --sh-basis legacy-real-even \
        --gyrus-a precentral --gyrus-b postcentral

Superior--middle temporal::

    swmtrack run --subject-dir /data/subject/freesurfer --fod /data/subject/fod.nii.gz \
        --output-dir output/temporal --assume-scanner-aligned --sh-basis legacy-real-even \
        --gyrus-a superiortemporal --gyrus-b middletemporal

Merged whole-parietal region::

    swmtrack run --subject-dir /data/subject/freesurfer --fod /data/subject/fod.nii.gz \
        --output-dir output/parietal --assume-scanner-aligned --sh-basis legacy-real-even \
        --endpoint-policy within_gyrus \
        --gyrus-a superiorparietal inferiorparietal postcentral supramarginal precuneus

A two-parcel intraparietal-region union uses only superiorparietal inferiorparietal
or the IDs 9221140 14433500. It does not isolate a named sulcus anatomically.

Whole hemisphere with endpoint classification::

    swmtrack run --subject-dir /data/subject/freesurfer --fod /data/subject/fod.nii.gz \
        --output-dir output/whole --assume-scanner-aligned --sh-basis legacy-real-even \
        --whole-hemisphere

Optional rigid registration replaces --assume-scanner-aligned with::

    --register-ants --b0 /data/subject/mean_b0.nii.gz

Or use corrected 4D DWI and b-values::

    --register-ants --dwi /data/subject/dwi.nii.gz --bvals /data/subject/bvals

Add --registration-python /path/to/ants-env/bin/python for an isolated ANTsPy
runtime. DWI/b0 must match the FOD spatial grid. See docs/end_to_end.rst.

Full-resolution label QC from the completed run::

    python scripts/export_qc.py --run-dir output/whole --output-dir qc/whole --inflated

Add --resume to an identical swmtrack run command to reuse verified stages.
Do not reuse a nonempty output directory after changing code, input or settings.

Synthetic tracking example
--------------------------

With the package installed, from the source archive run::

    python scripts/synthetic_demo.py --output-dir output/synthetic

This creates all inputs analytically, then runs preparation and tracking with
validated TCK/TRK/VTK exports. No HCP or FreeSurfer download is needed. This flat
mesh example does not run cortical segmentation. See docs/example_data.rst for
HCP S900 example-data terms and the recommended separate distribution route.
