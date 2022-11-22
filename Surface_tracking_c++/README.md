Contribution <br />
Author: Jialiang Ruan <br />
Created: 2022/11/10  <br />
Copyright: The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
and Informatics Institute of USC Laboratory of NeuroImaging USC Stevens Neuroimaging and Informatics Institute <br />
email: jruan@usc.edu

Please see the demo01.sh for the usage of an example from HCP data.

The main function of the tractography is Ufibertracking.cpp. <br />
Inputs: <br />
fod is the nii.gz data of the precomputed FOD volume <br />
mesh is the triangular mesh for half brain; the file type is obj <br />
roi is the ROI for tracking; it is binary data representing the vertices of the mesh <br />
sulcus is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the seed points <br />
gyrus1 is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary <br />
gyrus2 is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary <br />
output_trk is the output track file <br />
fodshphere is the output preprocess result file <br />

minFOD is the minimum FOD value we accept in reject sampling <br />
shrinkdis is the deformed distance of the WM mesh <br />
trackdegree is the upper bound degree for tracking  <br />
iteration is the repeating tracks for each seed triangle <br />
gyrusflag is the flag for Gyruslabel, Gyrusflag=0, for Gyruslabel1 not equals Gyruslabel2; Gyrusflag=1, for Gyruslabel1 equals Gyruslabel2 <br />
preprocess is option if you want do preprocess and output fodsphere file, 1 - 'YES', 0 - 'NO' <br />
track is option if you want do track and output tracks file, 1 - 'YES', 0 - 'NO' <br />


if you prefer a .tck file, please transfer using the function TrackVis2MRTrixTracks in matlab version.


Notice: The code will automatically use the parallel computation.
