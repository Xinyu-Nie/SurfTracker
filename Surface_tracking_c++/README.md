# SurfTracker C/C++ v1.0 
## Contribution <br />
Author: Jialiang Ruan <br />
Corresponding author: Xinyu Nie <br />
Created: 2023/5/10  <br />
Copyright: The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
and Informatics Institute of USC Laboratory of NeuroImaging USC Stevens Neuroimaging and Informatics Institute <br />
Email: jruan@usc.edu <br />
Corresponding email: xnie@usc.edu <br />

## Enviroments
### git lfs
### zlib for linux&macos, no zlib for win yet
### Visual Studio or MSVC for win
### boost lib for win&macos
## Compile options
### Windows<br />
./build_windows.bat
### MacOS<br/>
bash build_mac.sh
### LINUX
sh build_linux.sh<br/>
<!-- The main function of the tractography is Ufibertracking.cpp. <br /> -->

## Library and Compilation
The **Boost Library** is needed for compliance; please download it and put it on the right path. Here is the link to this library:<br />
https://www.boost.org/ <br />

Please also check if the CMake version in the script is the same as yours. <br />

# Demo Usage
### run ./demo.sh (./demo.bat for win)<br/>
### set dir=(your data directory)
### start with main (main.exe for win) </br>
### follow with -arg value
## Arguments <br />
**-fod**  fod is the nii.gz data of the precomputed FOD volume <br />
**-mesh** mesh is the triangular mesh for half brain; the file type is obj <br />
**-roi** roi is the ROI for tracking; it is binary data representing the vertices of the mesh <br />
**-sulcus** sulcus is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary  data to generate the seed points <br />
**-gyrus1** gyrus1 is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary <br />
**-gyrus2** gyrus2 is the binary data to label the sulcus region for the mesh; it will intersect with the ROI binary data to generate the stopping boundary <br />
**-output_trk** output_trk is the output track file <br />
**-fodsphere** fodshphere is the output preprocess result file <br />

**-minFOD** minFOD is the minimum FOD value we accept in reject sampling <br />
**-shrinkdis** shrinkdis is the deformed distance of the WM mesh <br />
**-trackdegree** trackdegree is the upper bound degree for tracking  <br />
**-iteration** iteration is the repeating tracks for each seed triangle <br />
**-gyrusflag** gyrusflag is the flag for Gyruslabel, Gyrusflag=0, for Gyruslabel1 not equals Gyruslabel2; Gyrusflag=1, for Gyruslabel1 equals Gyruslabel2 <br />
**-preprocess** preprocess is option if you want do preprocess and output fodsphere file, 1 - 'YES', 0 - 'NO' <br />
track is option if you want do track and output tracks file, 1 - 'YES', 0 - 'NO' <br />

## Tip
if you prefer a .tck file, please transfer using the function TrackVis2MRTrixTracks in matlab version.


## Notice: The code will automatically use the parallel computation.
