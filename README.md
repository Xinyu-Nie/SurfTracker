# SurfTracker
***A Probabilistic Tractography in Superficial White Matter:*** This method enables the reconstruction of the fiber tracts in Superficial White Matter, which has a sheer structure with high curvature. This tracking algorithm works on the triangular mesh of the white matter and overcomes the difficulties of conventional tractography.
<br />
The algorithm was developed in both C++ and Matlab versions. Please check the ***Readme*** and ***demo01*** in the corresponding folder for details of usage, and one subject data from the Human connectome projects (HCP) is uploaded by ifs for illustration; you can also use the google drive link blow to download the dataset.<br />
Google drive link: https://drive.google.com/file/d/1cjGqADdFo9s_4v5IIwN057iBGrawxx4B/view?usp=share_link

***Necessary  Data for Tractography:*** 
Fiber Orientation Distribution (FOD) data, A preregistered triangular mesh of the white matter (from T1 space to diffusion space), A binary label for ROI on the mesh, A binary label for seeds of tractography, one or two binary labels for stopping criteria. 




Acknowledgment: The Matlab codes that load the NIfTI document are from: Jimmy Shen (2022). Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), MATLAB Central File Exchange. Retrieved November 22, 2022.
