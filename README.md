# SurfTracker
A probabilistic tractography in Superficial White Matter. This method enables the reconstruction of the fiber tracts in Superficial White Matter, which has a sheer structure with high curvature. This tracking algorithm works on the triangular mesh of the white matter and overcomes the difficulties of conventional tractography.
<br />
The algorithm was developed in both C++ and Matlab versions. Please check the ***Readme*** and ***demo01*** in the corresponding folder for details of usage, and one subject data from the Human connectome projects (HCP) is used for illustration.

***Necessary  Data for Tractography:*** 
Fiber Orientation Distribution (FOD) data, A preregistered triangular mesh of the white matter (from T1 space to diffusion space), A binary label for ROI on the mesh, A binary label for seeds of tractography, one or two binary labels for stopping criteria. 

https://drive.google.com/file/d/1cjGqADdFo9s_4v5IIwN057iBGrawxx4B/view?usp=share_link
