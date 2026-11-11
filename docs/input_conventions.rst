Input and coordinate conventions
================================

FreeSurfer surface RAS, scanner RAS and FOD-world RAS are distinct. The importer
uses orig.mgz to convert surface RAS to scanner RAS and then applies the declared
scanner-to-FOD point transform. Supply exactly one of verified existing scanner
alignment, an explicit forward RAS 4x4 affine, or optional ANTs rigid fitting.
All surface points, tracking lengths, steps and offsets use physical millimeters.
Source vertex identities and original packed annotation IDs are retained.

No particular voxel storage order is required. RAS/LAS/LPS and oblique NIfTI
grids are sampled through their full voxel-to-world affine. The FOD is not
rewritten or resampled into a canonical RAS voxel array. RAS+ describes the
physical coordinate frame; it does not imply that anatomy and FOD are registered.
For prepare/track with an externally supplied VTK, its points must already be
in the FOD's physical world RAS mm; that route does not run the FreeSurfer importer.

Supported SH channels are even l, with m=-l..l in ascending order. The real basis
uses the Condon--Shortley phase; negative m is sqrt(2) Im(Y_l,abs(m)), positive m
is sqrt(2) Re(Y_l,m). Declare this through --sh-basis legacy-real-even. Spatial
sampling uses the full NIfTI affine. SH direction axes are a separate convention:
use --sh-to-world for an orthogonal 3x3 direction transform when required.
Channel count alone cannot identify a compatible basis. Convert incompatible
bases before using this pipeline.

Standard MRtrix3 WM FODs (dwi2fod csd or msmt_csd) use this supported orthonormal
basis and channel order. The option name legacy-real-even does not refer to the
older MRtrix 0.2 normalization. Export .mif to a 4D NIfTI with all SH channels::

    mrconvert wmfod.mif wmfod.nii.gz -datatype float32

.mif is not read directly by SurfTracker. Use the WM FOD, not GM/CSF maps.
An existing MRtrix3 NIfTI FOD can be used directly. Format conversion does not
change the basis or perform registration. See the
`MRtrix3 SH specification <https://mrtrix.readthedocs.io/en/latest/concepts/spherical_harmonics.html>`_.

Standard MRtrix3 FOD directions are expressed in scanner/world RAS. Leave
--sh-to-world at identity even when voxel axes are not RAS-ordered; the spatial
NIfTI affine is not an additional SH direction rotation. For a known different
frame, this orthogonal matrix maps a coefficient-frame direction to FOD-world
RAS. It does not convert between SH bases or align the anatomy. See
`MRtrix3 coordinates <https://mrtrix.readthedocs.io/en/latest/getting_started/image_data.html#coordinate-system>`_.

If an FOD is separately rotated or warped, the lobes must also be reoriented
with an FOD-aware tool. Changing headers or scalar resampling alone is not
sufficient; --sh-to-world is one global orthogonal transform, not a spatially
varying warp correction. See
`mrtransform FOD handling <https://mrtrix.readthedocs.io/en/latest/reference/commands/mrtransform.html>`_.

Native FreeSurfer annotation names or original packed color-table IDs select
parcels. These IDs are not volume LUT numbers. Within-region lists are explicitly
unioned into a separate binary target while the original annotation survives.
Whole-hemisphere mode selects the cortical annotation domain automatically.

Optional registration requires a 3D mean b0 or corrected 4D DWI plus matching
b-values in the supplied FOD grid. FOD estimation and DWI correction are outside
this pipeline. The rigid T1-to-b0 transform is inverted/converted from image-pull
LPS to forward point RAS; the FOD/SH image is not resampled. See end_to_end.rst.
