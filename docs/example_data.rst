Example data and HCP S900 redistribution
========================================

The software package contains no human-subject images, surfaces or tractograms.
The copyright holder of SurfTracker remains unspecified; the proposed software
license does not replace the terms governing example data.

Synthetic example (no dataset access needed)
--------------------------------------------

After installing the package, run this helper from the source archive::

    python scripts/synthetic_demo.py --output-dir output/synthetic

It generates a flat triangulated surface with explicit seed/endpoint masks and
a constant l=0 spherical-harmonic FOD, then runs the real preparation, tracking
and TCK/TRK/VTK export code. It checks both endpoints, physical length limits,
input-surface readback and tractogram readback. Inspect surface.vtk,
tracking_surface.vtk and tracks.vtk. The tracking surface has the default
0.5-mm inward offset. The example uses 256 attempts and one thread for speed;
the scientific tracking defaults are otherwise retained. This is an I/O and
tracking demonstration, not a cortical model or segmentation validation.

HCP S900: redistribution is permitted with conditions
-----------------------------------------------------

The WU-Minn HCP Open Access Data Use Terms dated April 26, 2013, explicitly
permit redistribution of original Open Access data and derived data under the
same terms (section 4). The S900 reference manual directs users to those terms.
This includes HCP-derived surfaces, FODs and tractograms; they must not be
presented as unrestricted data under the proposed BSD software license.

Every recipient working with shared HCP Open Access data must register with
ConnectomeDB and agree to the terms. The current HCP site links to ConnectomeDB
powered by BALSA for registration/acceptance. Recipients must also comply with
their institutional human-subject requirements and must not attempt participant
identification or contact. Restricted demographic, family, clinical or other
Restricted Data are not part of the proposed example.

Primary sources, checked November 10, 2026:

* https://www.humanconnectome.org/study/hcp-young-adult/data-use-terms
* https://www.humanconnectome.org/study/hcp-young-adult/document/wu-minn-hcp-consortium-open-access-data-use-terms
* https://www.humanconnectome.org/storage/app/media/data_use_terms/DataUseTerms-HCP-Open-Access-26Apr2013.pdf
* https://www.humanconnectome.org/storage/app/media/documentation/s900/HCP_S900_Release_Reference_Manual.pdf

Recommended distribution route
------------------------------

Keep human example data outside the wheel, source archive and ordinary Git
history. Use a separately versioned data deposit with a terms-acceptance access
step (for example, an appropriate HCP/BALSA deposit), or have users obtain the
original inputs from HCP and follow examples/README.rst. An unrestricted GitHub
asset download alone does not establish that recipients accepted the HCP terms.
These distribution choices are our recommendation, not a blanket HCP ban on
GitHub hosting. No public dataset deposit or upload is performed by this tree.

A useful first example would be one hemisphere of a confirmed S900 Open Access
subject, covering central and superior/middle temporal between-gyri tracking,
a merged parietal region, and whole-hemisphere endpoint classification. Use
original S900 inputs or record every subsequent preprocessing step explicitly.
The local validation directory named HCP2024 is not by itself evidence of S900
provenance; the example release must establish its actual source release first.

Include a minimal manifest listing relative filenames, SHA256 hashes, sizes,
subject/release, source package, processing software/version, transforms,
SH basis/order/axes, tracking settings, and whether each file is original or
derived. Do not publish institutional paths, credentials or unrestricted copies
of Restricted Data. Supply only required files: matching FreeSurfer surfaces,
annotation and orig.mgz, compatible FOD, and optional b0/brain.mgz if registration
is demonstrated. There is no need to distribute an entire subject directory.
Record the FOD estimator and any additional applicable terms.

Supply the unchanged HCP Open Access terms with the separate data archive and
direct recipients to register and accept them before using it. Preserve the
source subject identifier and processing provenance without adding Restricted
attributes. State that these are example outputs, not anatomical ground truth.

HCP acknowledgment
------------------

The HCP terms request the following acknowledgment for presentations and
publications benefiting from the data; include it in the example-data README:

    Data were provided [in part] by the Human Connectome Project, WU-Minn
    Consortium (Principal Investigators: David Van Essen and Kamil Ugurbil;
    1U54MH091657) funded by the 16 NIH Institutes and Centers that support the
    NIH Blueprint for Neuroscience Research; and by the McDonnell Center for
    Systems Neuroscience at Washington University.

Also cite the acquisition/preprocessing publications relevant to the selected
data, following https://www.humanconnectome.org/study/hcp-young-adult/document/hcp-citations .
The two SurfTracker method citations remain separate scientific credits.
