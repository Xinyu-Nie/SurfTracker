# export PATH=/ifshome/ibagari/gcc-7.2.0/:$PATH
# export PATH=/ifshome/ibagari/gcc-7.2.0/bin:$PATH
# export LD_LIBRARY_PATH=/ifshome/ibagari/gcc-7.2.0/lib64/:$LD_LIBRARY_PATH
# export CFLAGS=-fPIC
# export CC=/ifshome/ibagari/gcc-7.2.0/bin/gcc
# export CXX=/ifshome/ibagari/gcc-7.2.0/bin/g++
# alias gcc=${CC}
# alias "g++"=${CXX}
# make -B
subj_id=$1
# cnt=$2
# echo "cnt=$2";
echo "subject id is : $subj_id";
/ifs/loni/faculty/shi/spectrum/Student_2020/jruan/SurfaceTrack_paper/main \
       -fod /ifs/loni/faculty/shi/spectrum/yshi/ConnectomeAnalysis/FOD_Order16_ThreeCompModel/$subj_id\_FOD_Order16_WholeVol.nii.gz \
       -mesh /ifs/loni/faculty/shi/spectrum/yshi/CorticalPatternAnalysis/HCP/$subj_id/DTISpace/lh_white_DTISpace_rm.obj \
       -roi /ifs/loni/faculty/shi/spectrum/Student_2020/yuanli/SWM_UFiber/HCP/for_yao/$subj_id/DTISpace/UFiber/lh_MotorSensory_SulcalPatch_ROI_Label.raw \
       -sulcus /ifs/loni/faculty/shi/spectrum/Student_2020/yuanli/SWM_UFiber/HCP/for_yao/$subj_id/DTISpace/UFiber/lh_MotorSensory_SulcalPatch_SulcalLabel.raw \
       -gyrus1 /ifs/loni/faculty/shi/spectrum/Student_2020/yuanli/SWM_UFiber/HCP/for_yao/$subj_id/DTISpace/UFiber/lh_MotorSensory_SulcalPatch_GyrusLabel1.raw \
       -gyrus2 /ifs/loni/faculty/shi/spectrum/Student_2020/yuanli/SWM_UFiber/HCP/for_yao/$subj_id/DTISpace/UFiber/lh_MotorSensory_SulcalPatch_GyrusLabel2.raw \
       -output_trk /ifs/loni/faculty/shi/spectrum/Student_2020/jruan/SurfaceTrack_paper/$subj_id/test.trk \
       -fodsphere /ifs/loni/faculty/shi/spectrum/Student_2020/jruan/SurfaceTrack_paper/$subj_id/fodsphere.raw \
       -minFOD 0.1 \
       -shrinkdis 0.55 \
       -trackdegree 10 \
       -iteration 1 \
       -gyrusflag 0 \
       -preprocess 1 \
       -track 1