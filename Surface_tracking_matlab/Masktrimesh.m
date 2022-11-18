function [vertices,faces,faceOri] = Masktrimesh(faces0, vertices0,mask)
% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

Position=find(mask>0);
Lv=length(Position);
L2=size(vertices0,2);
vertices=zeros(Lv,L2);
for h=1:Lv
    vert=vertices0(Position(h),:);
    vertices(h,:)=vert;
end

[Lf1,Lf2]=size(faces0);
faces=zeros(Lf1,Lf2);
faceOri=zeros(Lf1,1);
t=1;
for k=1:Lf1
   v=faces0(k,:);
    if mask(v(1))+mask(v(2))+mask(v(3))>2.5
        a1=find(Position==v(1));
        a2=find(Position==v(2));
        a3=find(Position==v(3));
        faces(t,:)=[a1 a2 a3];
        faceOri(t)=k;
        t=t+1;
    end 
end
faces=faces(1:(t-1),:);
faceOri=faceOri(1:(t-1));
end