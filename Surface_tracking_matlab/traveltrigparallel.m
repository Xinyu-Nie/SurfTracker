function [outpoint,outedgef,terminflag,outbaryc] = traveltrigparallel(coordS,trgS,inpoint,prepoint,THETA,FOD,intri,pretri,inedge,minFOD,trackdegree)
% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu
threshold=cos(trackdegree*pi/180); % threshold for theata of streamline
toler=10^(-4);
terminflag=3;
outbaryc=[];
vin1=coordS(trgS(intri,1),:)';
vin2=coordS(trgS(intri,2),:)';
vin3=coordS(trgS(intri,3),:)';
v12=vin2-vin1;
v13=vin3-vin1;
v23=vin3-vin2;
q12=sqrt(sum(v12.^2));
v12=v12/q12;
q13=sqrt(sum(v13.^2));
v13=v13/q13;
q23=sqrt(sum(v23.^2));
v23=v23/q23;
e=[trgS(intri,1) trgS(intri,2);trgS(intri,3) trgS(intri,1);trgS(intri,2) trgS(intri,3)];
Nvin=cross(v12,v13);
qin=sqrt(sum(Nvin.^2));
Nvin=Nvin/qin;
v0=inpoint';
vpre1=coordS(trgS(pretri,1),:)';
vpre2=coordS(trgS(pretri,2),:)';
vpre3=coordS(trgS(pretri,3),:)';
v12pre=vpre2-vpre1;
v13pre=vpre3-vpre1;
v23pre=vpre3-vpre2;
q12=sqrt(sum(v12pre.^2));
v12pre=v12pre/q12;
q13=sqrt(sum(v13pre.^2));
v13pre=v13pre/q13;
q23=sqrt(sum(v23pre.^2));
v23pre=v23pre/q23;
Nvpre=cross(v12pre,v13pre);
qpre=sqrt(sum(Nvpre.^2));
Nvpre=Nvpre/qpre;
if inedge == 1
        edge=[trgS(pretri,1) trgS(pretri,2)];
    elseif inedge ==2
        edge=[trgS(pretri,3) trgS(pretri,1)];
    elseif inedge ==3
        edge=[trgS(pretri,2) trgS(pretri,3)];
end
    flage=zeros(3,1);
    for m=1:3
        em=e(m,:);
        u=intersect(em,edge);
        if length(u)>1.5
           flage(m)=1;
           break
        end
    end
    if sum(flage)==0
        error('nerbor error');
        terminflag=2;
    end
if em(1)==edge(1) && em(2)==edge(2)
    Nvpre=-Nvpre;
end
%Rotation matrix
direction1=cross(Nvpre,Nvin);
Len=sqrt(sum(direction1.^2));
direction=direction1/Len;
ux=direction(1);
uy=direction(2);
uz=direction(3);
cth=sum(Nvpre.*Nvin)/sqrt(sum(Nvin.^2)*sum(Nvpre.^2));
sth=sqrt(sum(direction1.^2))/sqrt(sum(Nvin.^2)*sum(Nvpre.^2));
R=[cth+(ux^2)*(1-cth) ux*uy*(1-cth)-uz*sth ux*uz*(1-cth)+uy*sth; uy*ux*(1-cth)+uz*sth cth+(uy^2)*(1-cth) uy*uz*(1-cth)-ux*sth; uz*ux*(1-cth)-uy*sth uz*uy*(1-cth)+ux*sth cth+(uz^2)*(1-cth)];
%
vfix=coordS(edge(1),:)';
rvector=prepoint'-vfix;
rvector=R*rvector;
pret=vfix+rvector;
invector=inpoint'-pret;
qinv=sqrt(sum(invector.^2));
if qinv > (10^(-10))
invector=invector/qinv;
end
if sum(invector.*Nvin)> toler 
     error('invector rotaion failed');
     terminflag=2;
end

N=1000;
tn=0;
while(terminflag==3)
   tn=tn+1;
if tn>N
    terminflag=1;
    break
end
cFOD=FOD(:,intri);
%alfa=Inversesampling(cFOD);
alfa=rejectsampling(cFOD,2,minFOD);
R12=rotationMatrix(Nvin,THETA(alfa));
vd=R12*v12;
qd=sqrt(sum(vd.^2));
if qd>0
vd=vd/qd;
end
if sum(vd.*Nvin)> toler
     error('vector rotaion failed');
end
    a1=sum(cross(vd,v12).^2);
    a2=sum(cross(vd,v13).^2);
    a3=sum(cross(vd,v23).^2);
     if a1>toler && m~=1
        A1=[vd (vin1-vin2)];
        b1=vin1-v0;
        [x1,flagx]=lsqr(A1,b1);
        if x1(2)>= 0 && x1(2)<= 1 && x1(1)>0
            outpoint=(x1(2)*vin2+(1-x1(2))*vin1)';
            outedgef=1;
            outbaryc=[intri (1-x1(2)) x1(2) 0];
        end
    end
     if a2>toler && m~=2
        A1=[vd (vin1-vin3)];
        b1=vin1-v0;
        [x1,flagx]=lsqr(A1,b1);
        if x1(2)>= 0 && x1(2)<= 1 && x1(1)>0
          outpoint=(x1(2)*vin3+(1-x1(2))*vin1)';
          outedgef=2;
          outbaryc=[intri (1-x1(2)) 0 x1(2)];
        end
    end
      if a3>toler && m~=3
        A1=[vd (vin2-vin3)];
        b1=vin2-v0;
        [x1,flagx]=lsqr(A1,b1);
        if x1(2)>= 0 && x1(2)<= 1 && x1(1)>0
             outpoint=(x1(2)*vin3+(1-x1(2))*vin2)';
             outedgef=3;
             outbaryc=[intri 0 (1-x1(2)) x1(2)];
        end
      end
      if exist('outpoint')
         vout=(outpoint-inpoint)';
         qout=sqrt(sum(vout.^2));
         vout=vout/qout;
          if sum(vout.*invector)> threshold
              terminflag=0;
              break
          end
      end
end

if ~ exist('outedgef')
    outedgef=0;
end

if ~ exist('outpoint')
    outpoint=inpoint;
end
