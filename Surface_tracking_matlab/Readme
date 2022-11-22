Contribution
Author: Xinyu Nie
Created: 2022/11/10
Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
and Informatics Institute of USC Laboratory of NeuroImaging USC Stevens Neuroimaging and Informatics Institute
email: xnie@usc.edu

Please see the demo01.m for the usage of an example from HCP data.

The main function of the tractography is RunUFiberTrackingnew.
Inputs:
FOD is the nii.gz data of the precomputed FOD volume.
ROI is the ROI for tracking; it is binary data representing the vertices of the mesh.
Mesh is the triangular mesh for half brain; the file type is obj.
Sulcuslabel is a binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the seed points
Gyruslabel1 is a binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary
Gyruslabel2 is a binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary
Num of iteration is the repeating tracks for each seed triangle; if you need more streamlines, increase it.
FOD_THD is the minimum FOD value we accept in reject sampling
Trackdegree is the upper bound degree for tracking 
Gyrusflag is the flag for Gyruslabel, set Gyrusflag=0, for Gyruslabel1 not equals Gyruslabel2; set Gyrusflag=1, for Gyruslabel1 equals Gyruslabel2
ShrinkDis is the deformed distance of the WM mesh

output is the .trk file for tractography; if you prefer a .tck file, please transfer using the function TrackVis2MRTrixTracks.


Notice: The code will automatically use the parallel computation in Matlab.
