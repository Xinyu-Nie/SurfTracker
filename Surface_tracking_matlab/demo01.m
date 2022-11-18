clear all
close all

%FOD is the nii.gz data of the precomputed FOD volume.
%ROI is the ROI for tracking; it is binary data representing the vertices of the mesh.
%Mesh is the triangular mesh for half brain; the file type is obj.
%Sulcuslabel is a binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the seed points
%Gyruslabel1 is a binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary
%Gyruslabel2 is a binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary
%Num of iteration is the repeating tracks for each seed triangle
%FOD_THD is the minimum FOD value we accept in reject sampling
%Trackdegree is the upper bound degree for tracking 
%Gyrusflag is the flag for Gyruslabel, Gyrusflag=0, for Gyruslabel1 not equals Gyruslabel2; Gyrusflag=1, for Gyruslabel1 equals Gyruslabel2
%ShrinkDis is the deformed distance of the WM mesh

%output is the .trk file for tractography; if you prefer a .tck file, please transfer using the function TrackVis2MRTrixTracks.

NumSamplePerSeed = '6';
FOD_THD = '0.01';
Angle_THD  = '7';
GyrusFlag = '0';
ShrinkDis = '0.5';

FOD = '..\101107\101107_FOD_Order16_WholeVol.nii.gz';
Mesh= '..\101107\101107_white_lh.obj';
SulcusLabel = '..\101107\101107_white_lh_PrePostcentralCommon2mm.raw';
GyrusLabel1 = '..\101107\101107_white_lh_PrecentralGyralSkeletons_dilated.raw';
GyrusLabel2 = '..\101107\101107_white_lh_PostcentralGyralSkeletons_dilated.raw';
ROIonMesh = '..\101107\101107_white_lh_PrePostcentralLabel.raw';
FODProjOntoROI='..\101107\101107_FOD2D_ROI.mat';
OutputTract_Imagespace= '..\101107\101107_UFiber_Pre2Post.trk';
OutputTract_Physicalspace= '..\101107\101107_UFiber_Pre2Post.tck';


RunUFiberTrackingnew(Mesh,FOD,ROIonMesh,FODProjOntoROI, ShrinkDis, SulcusLabel,GyrusLabel1,GyrusLabel2,OutputTract_Imagespace,NumSamplePerSeed,FOD_THD, Angle_THD, GyrusFlag );

TrackVis2MRTrixTracks(OutputTract_Imagespace, FOD,OutputTract_Physicalspace );