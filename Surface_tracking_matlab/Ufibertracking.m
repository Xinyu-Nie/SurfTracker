function Ufibertracking(Mesh,ROI,Sulcuslabel,Gyruslabel1,Gyruslabel2,output,iteration,coordS,trgS,TRIcen,THETA,FODProjected,trinerbor,dim,res,minFOD,trackdegree,Gyrusflag)
%ROI is the ROI to track, it is a binary data for the vertices of the Mesh
%Mesh is the triangular mesh for half brain, the file is obj.
%Sulcuslabel is a binary data to label the sulcus region for the Mesh, it
%will intersect with the ROI binary data to generate the seed points
%Gyruslabel1 is a binary data to label the sulcus region for the Mesh, it
%will intersect with the ROI binary data to generate the stopping points
%Gyruslabel2 is a binary data to label the sulcus region for the Mesh, it
%will intersect with the ROI binary data to generate the stopping points
%output is the .trk file for tractogram
%iteration is the repeating times for each seed triangle
%FODmin is the minmium fod value we accept.
%Trackdegree is the upper bound degree for tracking 
%Gyrusflag is the flag for Gyruslabel, Gyrusflag=0, for Gyruslabel1 ~= Gyruslabel2; Gyrusflag=1, for Gyruslabel1 == Gyruslabel2
%Other input generated from the SurfaceFODpreprossing function

% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

rng shuffle
numseed=10000; %An upper bound for the number of seed points
fid = fopen(ROI);
MTL = fread(fid,'int32');
fclose(fid);
fid = fopen(Gyruslabel1);
Gyrus1 = fread(fid,'int32');
fclose(fid);
fid = fopen(Gyruslabel2);
Gyrus2 = fread(fid,'int32');
fclose(fid);
fid = fopen(Sulcuslabel);
Suc = fread(fid,'int32');
fclose(fid);
if length(MTL) ~= length(Gyrus1) || length(MTL) ~= length(Suc)
    error('The labels do not match')
end
[coord,nS,trg] = ReadObjShape(Mesh);
[SucdeepS,f1,Faceori]=Masktrimesh(trg,Suc,MTL);
GY1=Masktrimesh(trg,Gyrus1,MTL);
GY2=Masktrimesh(trg,Gyrus2,MTL);
Ltri=size(trgS,1);
GyruLabel1=zeros(Ltri,1);
GyruLabel2=GyruLabel1;
SeedLabel=zeros(Ltri,1);
for k=1:Ltri
 triK=trgS(k,:);
 for m=1:3
     if GY1(triK(m))>0
      GyruLabel1(k)=1;
     end
     if GY2(triK(m))>0
        GyruLabel2(k)=2;
     end
     if SucdeepS(triK(m))>0
         SeedLabel(k)=1;
     end
 end
end

SeedP=find(SeedLabel>0);
numSuc=numel(SeedP);
if numSuc<numseed
    seedflag=SeedP;
else
  flag=randperm(numSuc,numseed);
  seedflag=sort(SeedP(flag));
end
 Biflag=zeros(length(seedflag)*iteration,2);
 cuvature=cell(length(seedflag)*iteration,1);
 PtV=cuvature;
 Tricoord=PtV;
 Ltract=zeros(length(seedflag)*iteration,1);
 Nseed=length(seedflag);
 flagend=Biflag;
 for m=1:iteration
     bm=(m-1)*Nseed;
   parfor t=1:Nseed
     [tractline,flagendinge,flagendingp,biterminflag,TriBaryC]=trackingonsufaceparallel(coordS,trgS,TRIcen,THETA,FODProjected,seedflag(t),trinerbor,GyruLabel2,GyruLabel1,minFOD,trackdegree);
     Ltract(bm+t)=size(tractline,1);
     Biflag(bm+t,:)=biterminflag;
     PtV{bm+t}=tractline;
     flagend(bm+t,:)=[flagendinge flagendingp];
   end
 end
if Gyrusflag==1
 flag=0;
 att=1;
 for k=1:length( Ltract)
     lu=Ltract(k);
     if lu>20
         if Biflag(k,1)==4 || Biflag(k,1)==4.5 
             if Biflag(k,2)==4  || Biflag(k,2)==4.5     
 if flag>0 
         LineInd(end+1)=LineInd(end)+lu;
         PtVec(end+1:end+lu,:)=PtV{k};
         else
         LineInd=[lu];
         PtVec=PtV{k}; 
         flag=1;
 end
             end
         end
     end
  att=att+1;
 end
elseif Gyrusflag==0
    flag=0;
 att=1;
 for k=1:length( Ltract)
     lu=Ltract(k);
     if lu>20
         if Biflag(k,1)==4 || Biflag(k,1)==4.5 
             if Biflag(k,2)==4  || Biflag(k,2)==4.5 
              if Biflag(k,1)~=Biflag(k,2)
 if flag>0 
         LineInd(end+1)=LineInd(end)+lu;
         PtVec(end+1:end+lu,:)=PtV{k};    
         else
         LineInd=[lu];
         PtVec=PtV{k}; 
         flag=1;
 end
               end
              end
         end
     end
  att=att+1;
 end
end

 PtVec1=zeros(size(PtVec));
 PtVec1(:,end)=PtVec(:,end);
 for i = 1:3
    PtVec1(:,i) = (PtVec(:,i)-1).*res(i);
 end
 LineInd1=LineInd;
WriteTrackVisFile2(LineInd1',PtVec1,dim(1:3),res,output);
end