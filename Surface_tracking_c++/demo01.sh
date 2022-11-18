# environment: linux

# fod is the nii.gz data of the precomputed FOD volume
# mesh is the triangular mesh for half brain; the file type is obj
# roi is the ROI for tracking; it is binary data representing the vertices of the mesh
# sulcus is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the seed points
# gyrus1 is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary
# gyrus2 is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary
# output_trk is the output track file
# fodshphere is the output preprocess result file

# minFOD is the minimum FOD value we accept in reject sampling
# shrinkdis is the deformed distance of the WM mesh
# trackdegree is the upper bound degree for tracking 
# iteration is the repeating tracks for each seed triangle
# gyrusflag is the flag for Gyruslabel, Gyrusflag=0, for Gyruslabel1 not equals Gyruslabel2; Gyrusflag=1, for Gyruslabel1 equals Gyruslabel2
# preprocess is option if you want do preprocess and output fodsphere file, 1 - 'YES', 0 - 'NO'
# track is option if you want do track and output tracks file, 1 - 'YES', 0 - 'NO'
dir=$(dirname $0)
$dir/main \
    -fod $dir/../101107/101107_FOD_Order16_WholeVol.nii.gz \
    -mesh $dir/../101107/101107_white_lh.obj \
    -roi $dir/../101107/101107_white_lh_PrePostcentralLabel.raw \
    -sulcus $dir/../101107/101107_white_lh_PrePostcentralCommon2mm.raw \
    -gyrus1 $dir/../101107/101107_white_lh_PrecentralGyralSkeletons_dilated.raw \
    -gyrus2 $dir/../101107/101107_white_lh_PostcentralGyralSkeletons_dilated.raw \
    -output_trk $dir/../101107/101107_UFiber_Pre2Post.trk \
    -fodsphere $dir/../101107/101107_FOD2D_ROI.raw \
    -minFOD 0.01 \
    -shrinkdis 0.5 \
    -trackdegree 7 \
    -iteration 6 \
    -gyrusflag 0 \
    -preprocess 1 \
    -track 1
