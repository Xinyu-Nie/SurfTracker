# SurfTracker
***A Probabilistic Tractography in Superficial White Matter:*** This algorithm enables the reconstruction of the fiber tracts in Superficial White Matter (SWM), which has a sheet structure with high curvature. This tracking algorithm works on the triangular mesh representation of the SWM and overcomes  difficulties in conventional volume-based tractography.
<br />
The algorithm was developed in both C++ and Matlab versions. Please check the ***Readme*** and ***demo*** in the corresponding folder for usage details.  We provide scripts and libraries to help compile the C++ code for Linux, Mac, and Windows. <br />
One subject data from the Human Connectome projects (HCP) is uploaded by ifs for illustration, and you can also use the Google Drive link below to download the dataset.<br />
https://drive.google.com/file/d/1cjGqADdFo9s_4v5IIwN057iBGrawxx4B/view?usp=share_link

***Necessary  Data for Tractography:*** 
Fiber Orientation Distribution (FOD) data, A preregistered triangular mesh of the white matter (from T1 space to diffusion space), A binary label for ROI on the mesh, A binary label for seeds of tractography, one or two binary labels for stopping criteria. 

***Publication:*** X. Nie, J. Ruan, M. C. G. Otaduy, L. T. Grinberg, J. Ringman, and Y. Shi, "Surface-based Probabilistic Fiber Tracking in Superficial White Matter," in IEEE Transactions on Medical Imaging, doi: 10.1109/TMI.2023.3329451.


Acknowledgment: The Matlab codes that load the NIfTI document are from Jimmy Shen (2022). Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).
