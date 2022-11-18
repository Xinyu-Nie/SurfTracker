function [tract,flagendinge,flagendingp,biterminflag, TriBaryC] = trackingonsufaceparallel(coordS,trgS,tricenter,THETA,FOD,seedflag,trinerbor,Triendflag,Tribeginflag,minFOD,trackdegree)
% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu
seedpoint=tricenter(seedflag,:);
seedFOD=FOD(:,seedflag);
tractline=seedpoint;
tractbidirection=seedpoint;
alfa=rejectsampling(seedFOD,2,minFOD);
    v1=coordS(trgS(seedflag,1),:)';
    v2=coordS(trgS(seedflag,2),:)';
    v3=coordS(trgS(seedflag,3),:)';
    v0=seedpoint';
    v12=v2-v1;
    v13=v3-v1;
    v23=v3-v2;
    Nv=cross(v12,v13);
    q=sqrt(sum(Nv.^2));
if q>0
    Nv=Nv/q;
    q1=sqrt(sum(v12.^2));
    v12=v12/q1;
    R=rotationMatrix(Nv,THETA(alfa));
    vd=R*v12;
    a1=sum(cross(vd,v12).^2);
    a2=sum(cross(vd,v13).^2);
    a3=sum(cross(vd,v23).^2);
    if a1>0
        A1=[vd (v1-v2)];
        b1=v1-v0;
        [x1,flagx]=lsqr(A1,b1);
        if x1(2)>= 0 && x1(2)<= 1
            if x1(1)>0
                tractline(end+1,:)=(x1(2)*v2+(1-x1(2))*v1)';
                cpedge=1;
                BaryP=[seedflag (1-x1(2)) x1(2) 0];
            else
                tractbidirection(end+1,:)=(x1(2)*v2+(1-x1(2))*v1)';
                cnedge=1;
                BaryN=[seedflag (1-x1(2)) x1(2) 0];
            end
        end
    end
     if a2>0
        A1=[vd (v1-v3)];
        b1=v1-v0;
        [x1,flagx]=lsqr(A1,b1);
        if x1(2)>= 0 && x1(2)<= 1
            if x1(1)>0
                tractline(end+1,:)=(x1(2)*v3+(1-x1(2))*v1)';
                cpedge=2;
                BaryP=[seedflag (1-x1(2)) 0 x1(2)];
            else
                tractbidirection(end+1,:)=(x1(2)*v3+(1-x1(2))*v1)';
                cnedge=2;
                BaryN=[seedflag (1-x1(2)) 0 x1(2)];
            end
        end
    end
      if a3>0
        A1=[vd (v2-v3)];
        b1=v2-v0;
        [x1,flagx]=lsqr(A1,b1);
        if x1(2)>= 0 && x1(2)<= 1
            if x1(1)>0
                tractline(end+1,:)=(x1(2)*v3+(1-x1(2))*v2)';
                cpedge=3;
                BaryP=[seedflag 0 (1-x1(2)) x1(2)];
            else
                tractbidirection(end+1,:)=(x1(2)*v3+(1-x1(2))*v2)';
                cnedge=3;
                BaryN=[seedflag 0 (1-x1(2)) x1(2)];
            end
        end
      end
flagpd=0;
flagnd=0;
currentpositivetrig=seedflag;
currentnegativetrig=seedflag;
Ptri=seedflag;
Ntri=seedflag;
biterminflag=zeros(1,2);
flap=1;
flagendingp=0;
terminflag=0;
BaryP0=[];
BaryN0=[];
if exist('cpedge')
    BaryP0=BaryP;
while (flagpd==0 &&  cpedge~=0 && flap<1000 && terminflag==0)
    Tric=trinerbor(currentpositivetrig,:);
    if Tric(cpedge)>0
    if Triendflag(Tric(cpedge))>0
      flagendingp=1;
      biterminflag(1)=4;
      break
    end
    if  Tribeginflag(Tric(cpedge))>0
      flagendingp=2;
      biterminflag(1)=4.5;
      break
    end
    [outpoint,outedgef,terminflag,outbaryc] = traveltrigparallel(coordS,trgS,tractline(end,:),tractline(end-1,:),THETA,FOD,Tric(cpedge),currentpositivetrig,cpedge,minFOD,trackdegree);
    currentpositivetrig=Tric(cpedge);
    cpedge=outedgef;
    if outedgef>0
    tractline(end+1,:)=outpoint;
    Ptri(end+1)=currentpositivetrig;
    BaryP0(end+1,:)=outbaryc;
    flap=flap+1;
    end
    else
    flagpd=1; 
    biterminflag(1)=5;
    end
end
if biterminflag(1)==0
if terminflag>0
biterminflag(1)=terminflag;
elseif flap==1000
       biterminflag(1)=6;
end
end
else
    biterminflag(1)=7;
end
flae=1;
flagendinge=0;
terminflag=0;
if exist('cnedge')
    BaryN0=BaryN;
while (flagnd==0 &&  cnedge~=0 && flae<1000 && terminflag==0)
   Tric=trinerbor(currentnegativetrig,:);
    if Tric(cnedge)>0
    if Triendflag(Tric(cnedge))>0
      flagendinge=1;
      biterminflag(2)=4;
      break
    end
    if Tribeginflag(Tric(cnedge))>0
      flagendinge=2;
      biterminflag(2)=4.5;
      break
    end
    [outpoint,outedgef,terminflag,outbaryc] = traveltrigparallel(coordS,trgS,tractbidirection(end,:),tractbidirection(end-1,:),THETA,FOD,Tric(cnedge),currentnegativetrig,cnedge,minFOD,trackdegree);
    currentnegativetrig=Tric(cnedge);
    cnedge=outedgef;
    if outedgef>0
    tractbidirection(end+1,:)=outpoint;
    Ntri(end+1)=currentnegativetrig;
    BaryN0(end+1,:)=outbaryc;
    flae=flae+1;
    end
   else
     flagnd=1; 
     biterminflag(2)=5;
    end
end
if biterminflag(2)==0
if terminflag>0
biterminflag(2)=terminflag;
elseif flae==1000
       biterminflag(2)=6;
end
end
else
    biterminflag(2)=7;
end

tractbidirection=flip(tractbidirection,1);
tract=tractbidirection;
L=size(tractline,1);
if L>1
tract(end+1:end+L-1,:)=tractline(2:end,:);
end
Ntri=flip(Ntri);
TriCount=Ntri;
Lt=length(Ptri);
Ln=size(BaryN0,1);
Lp=size(BaryP0,1);
TriBaryC=[];
if Ln>0 && Lp >0
TriBaryC=flip(BaryN0,1);
TriBaryC(end+1:end+Lp,:)=BaryP0;
elseif  Lp >0
  TriBaryC=BaryP0;
end

if Lt>0
 TriCount(end+1:end+Lt)=Ptri(1:end);
end
  
end
