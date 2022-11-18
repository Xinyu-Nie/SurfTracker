function RunUFiberTrackingnew(Mesh,FOD,ROIonMesh,FODProjOntoROI, ShrinkDis, SulcusLabel,GyrusLabel1,GyrusLabel2,OutputTract,NumSamplePerSeed,FOD_THD, Angle_THD, GyrusFlag)
% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

ShrinkDis = str2num(ShrinkDis);

SurfaceFODPronew(FOD,Mesh,ROIonMesh,FODProjOntoROI,ShrinkDis);

load(FODProjOntoROI);
NumSamplePerSeed = str2num(NumSamplePerSeed);
FOD_THD = str2num(FOD_THD);
Angle_THD = str2num(Angle_THD);
GyrusFlag = str2num(GyrusFlag);
Ufibertracking(Mesh,ROIonMesh,SulcusLabel,GyrusLabel1,GyrusLabel2,OutputTract,NumSamplePerSeed,coordS,trgS,TRIcen,THETA,FODProjected,trinerbor,dim,res,FOD_THD,Angle_THD,GyrusFlag);

end