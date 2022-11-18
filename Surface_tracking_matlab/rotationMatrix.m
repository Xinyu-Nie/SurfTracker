function R = rotationMatrix(direction,theta)
% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

ux=direction(1);
uy=direction(2);
uz=direction(3);
cth=cos(theta);
sth=sin(theta);
R=[cth+(ux^2)*(1-cth) ux*uy*(1-cth)-uz*sth ux*uz*(1-cth)+uy*sth; uy*ux*(1-cth)+uz*sth cth+(uy^2)*(1-cth) uy*uz*(1-cth)-ux*sth; uz*ux*(1-cth)-uy*sth uz*uy*(1-cth)+ux*sth cth+(uz^2)*(1-cth)];
end