function SurfaceFODPronew(FOD,Mesh,ROI,ouput,shrinkdis)
%FOD is the FOD data
%ROI is the ROI to track, it is a binary data for the vertices of the Mesh
%Mesh is the triangular mesh for half brain, the file type is obj.
%output is the .mat computed FOD2D data on the triangular mesh
%shrinkdis is the distace to shrink the cortical surface into the WM

% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

fid = fopen(ROI);
MTL = fread(fid,'int32');
fclose(fid);
nii = load_untouch_nii(FOD);
N = size(nii.img,4);
res = nii.hdr.dime.pixdim(2:4);
[coord,nS,trg] = ReadObjShape(Mesh);
coord = coord - shrinkdis*nS;
for i = 1:3
    coord(:,i) = coord(:,i)./res(i) + 1;
end
[coordS,trgS] = Masktrimesh(trg,coord,MTL);
L1=size(trgS,1);
fai=pi/90;
theta=pi/90;
FAI=(pi/180):fai:pi;
THETA=(pi/180):theta:(2*pi);
k1=length(FAI);
k2=length(THETA);
weight=sin(FAI);
K=k1*k2;
TRIcen=zeros(size(trgS));
for t=1:L1
 TRIcen(t,:)=(coordS(trgS(t,1),:)+coordS(trgS(t,2),:)+coordS(trgS(t,3),:))/3;
end
x = InterpFODontoTRI(FOD, TRIcen);
SPHARM_Order = (sqrt(8*N+1)-3)/2;
FODProjected=zeros(k2,L1);

FODcoordinate=zeros(3,K);
for h2=1:k2
        R2=rotationMatrix([0 0 1]',THETA(h2));
        dr2=R2*[1 0 0]';
    for h1=1:k1
            Nvv=cross([0 0 1]',dr2);
            q2=sqrt(sum(Nvv.^2));
            Nvv=Nvv/q2;
          R1=rotationMatrix(Nvv,FAI(h1));  
          dr1=R1*[0 0 1]';
          qq=sqrt(sum(dr1.^2));
          dr1=dr1/qq;
          FODcoordinate(:,(h2-1)*k1+h1)=dr1;
     end
end
 BS = matrixB(FODcoordinate',SPHARM_Order);

%figure
%scatter3(FODcoordinate(1,:)',FODcoordinate(2,:)',FODcoordinate(3,:)')

R2C = Real2complex(SPHARM_Order);
C2R = Complex2real(SPHARM_Order);
R2C1=transpose(R2C);
C2R1=transpose(C2R);
thetay=pi/2;
UY90p= SpYrotationnew (SPHARM_Order,thetay);
thetay=-pi/2;
UY90n= SpYrotationnew (SPHARM_Order,thetay);
 
parfor t=1:L1
    v1=(coordS(trgS(t,2),:)-coordS(trgS(t,1),:))';
    v2=(coordS(trgS(t,3),:)-coordS(trgS(t,1),:))';
    Nv=cross(v1,v2);
    q=sqrt(sum(Nv.^2));
    if q>0
    Nv=Nv/q;
    q1=sqrt(sum(v1.^2));
    v1=v1/q1;
    vy=cross(Nv,v1);
    qy=sqrt(sum(vy.^2));
    if qy>0
       vy=vy/qy;
    end
    Rt=[v1 vy Nv];
    [alfa,beta,gamma] = DecomrotationPro(Rt);
    Ylr= x(t,:)';
    Ylc= C2R1*Ylr;
    Ualfa = SpZrotation(SPHARM_Order,alfa+pi/2);
    Ylc1=Ualfa*Ylc;
    UZbeta=SpZrotation(SPHARM_Order,beta);
    Ylc2=UY90n*UZbeta*UY90p*Ylc1;
    Ugamma = SpZrotation(SPHARM_Order,gamma-pi/2);
    Ylc3=Ugamma*Ylc2;
    Ylr3=R2C1*Ylc3;
    Ylr3= real(Ylr3);
    FODPROJ=BS*Ylr3;
    FODPROJ=max(FODPROJ,0);
      for h=1:k2
        fodsphere=FODPROJ((h-1)*k1+1:h*k1);
        Fp=sum((weight').*fodsphere);
        FODProjected(h,t)= Fp*fai;
     end
    end
end

 [Vert2FaceMap,TrgNbrMap] = BuildMeshNbrhood(coordS,trgS);
 trinerbor=TrgNbrMap;
 trinerbor(:,1)=TrgNbrMap(:,3);
 trinerbor(:,3)=TrgNbrMap(:,1);
 dim=size(nii.img);
 save(ouput,'coordS','trgS','TRIcen','THETA','FODProjected','trinerbor','res','dim');
end