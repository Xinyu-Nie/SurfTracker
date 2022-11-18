function x = InterpFODontoTRI(FOD, Tri)
% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

nii = load_untouch_nii(FOD);
N = size(nii.img,4);
coord=Tri;
x = zeros(size(coord,1),N);
for i = 1:N
    x(:,i) = interp3(permute(nii.img(:,:,:,i),[2 1 3]),coord(:,1),coord(:,2),coord(:,3));
end