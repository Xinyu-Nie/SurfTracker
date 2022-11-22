% Contribution
%  Author: Jialiang Ruan
%  Created: 2022/11/10
%  Copyright: The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
and Informatics Institute of USC Laboratory of NeuroImaging USC Stevens Neuroimaging and Informatics Institute
%  email: jruan@usc.edu

Please see the demo01.sh for the usage of an example from HCP data.

The main function of the tractography is Ufibertracking.cpp.
Inputs:
fod is the nii.gz data of the precomputed FOD volume
mesh is the triangular mesh for half brain; the file type is obj
roi is the ROI for tracking; it is binary data representing the vertices of the mesh
sulcus is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the seed points
gyrus1 is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary
gyrus2 is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary
output_trk is the output track file
fodshphere is the output preprocess result file

minFOD is the minimum FOD value we accept in reject sampling
shrinkdis is the deformed distance of the WM mesh
trackdegree is the upper bound degree for tracking 
iteration is the repeating tracks for each seed triangle
gyrusflag is the flag for Gyruslabel, Gyrusflag=0, for Gyruslabel1 not equals Gyruslabel2; Gyrusflag=1, for Gyruslabel1 equals Gyruslabel2
preprocess is option if you want do preprocess and output fodsphere file, 1 - 'YES', 0 - 'NO'
track is option if you want do track and output tracks file, 1 - 'YES', 0 - 'NO'


if you prefer a .tck file, please transfer using the function TrackVis2MRTrixTracks in matlab version.


Notice: The code will automatically use the parallel computation.
