dir=../
./main -fod $dir/101107/101107_FOD_Order16_WholeVol.nii \
    -mesh $dir/101107/101107_white_lh.obj \
    -roi $dir/101107/101107_white_lh_PrePostcentralLabel.raw \
    -sulcus $dir/101107/101107_white_lh_PrePostcentralCommon2mm.raw \
    -gyrus1 $dir/101107/101107_white_lh_PrecentralGyralSkeletons_dilated.raw \
    -gyrus2 $dir/101107/101107_white_lh_PostcentralGyralSkeletons_dilated.raw \
    -output_trk $dir/101107/101107_UFiber_Pre2Post.trk \
    -fodsphere $dir/101107/101107_FOD2D_ROI.raw \
    -minFOD 0.01 \
    -shrinkdis 0.5 \
    -trackdegree 7 \
    -iteration 1 \
    -gyrusflag 0 \
    -preprocess 1 \
    -track 1
