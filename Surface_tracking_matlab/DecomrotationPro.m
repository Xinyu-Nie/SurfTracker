function [alfa,beta,gamma] = DecomrotationPro(R)
% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

xr=R(:,1);
%yr=R(:,2);
zr=R(:,3);
xq=sqrt(sum(xr.^2));
xr=xr/xq;
%yq=sqrt(sum(yr.^2));
%yr=yr/yq;
zq=sqrt(sum(zr.^2));
zr=zr/zq;
x=[1;0;0];
y=[0;1;0];
z=[0;0;1];
n=cross(z,zr);
qn=sqrt(sum(n.^2));
if  qn>0
    n=n/qn;
end
  alfa=acos(sum(n.*y));
  if n(1)>0
      alfa =2*pi - alfa;
  end
  Rz=[cos(alfa) -sin(alfa) 0; sin(alfa) cos(alfa) 0; 0 0 1];
  x1=Rz*x;
  %y1=n;
  beta=acos(sum(z.*zr));
  Rn=rotationMatrix(n,beta);
  zr1=Rn*z;
  pz1=abs(sum(zr1.*zr));
  Rnr=rotationMatrix(n,2*pi-beta);
  zr11=Rnr*z;
  pz2=abs(sum(zr11.*zr));
  if pz1 < pz2
      beta = 2*pi-beta;
  end
  x2=rotationMatrix(n,beta)*x1;
  gamma=acos(sum(x2.*xr));
  Rzr=rotationMatrix(zr,gamma);
  xr1=Rzr*x2;
  px1=abs(sum(xr1.*xr));
  Rzrr=rotationMatrix(zr,2*pi-gamma);
  xr11=Rzrr*x2;
  px2=abs(sum(xr11.*xr));
  if px1< px2
      gamma=2*pi -gamma;
  end

end