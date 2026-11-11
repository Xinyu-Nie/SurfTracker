Input and coordinate conventions
================================

FreeSurfer surface RAS, scanner RAS and FOD-world RAS are distinct. The importer
uses orig.mgz to convert surface RAS to scanner RAS and then applies the declared
scanner-to-FOD point transform. Supply exactly one of verified existing scanner
alignment, an explicit forward RAS 4x4 affine, or optional ANTs rigid fitting.
All surface points, tracking lengths, steps and offsets use physical millimeters.
Source vertex identities and original packed annotation IDs are retained.

Supported SH channels are even l, with m=-l..l in ascending order. The real basis
uses the Condon--Shortley phase; negative m is sqrt(2) Im(Y_l,abs(m)), positive m
is sqrt(2) Re(Y_l,m). Declare this through --sh-basis legacy-real-even. Spatial
sampling uses the full NIfTI affine. SH direction axes are a separate convention:
use --sh-to-world for an orthogonal 3x3 direction transform when required.
Channel count alone cannot identify a compatible basis. Convert incompatible
bases before using this pipeline.

Native FreeSurfer annotation names or original packed color-table IDs select
parcels. These IDs are not volume LUT numbers. Within-region lists are explicitly
unioned into a separate binary target while the original annotation survives.
Whole-hemisphere mode selects the cortical annotation domain automatically.

Optional registration requires a 3D mean b0 or corrected 4D DWI plus matching
b-values in the supplied FOD grid. FOD estimation and DWI correction are outside
this pipeline. The rigid T1-to-b0 transform is inverted/converted from image-pull
LPS to forward point RAS; the FOD/SH image is not resampled. See end_to_end.rst.
