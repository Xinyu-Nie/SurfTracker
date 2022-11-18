function direction = rejectsampling(FOD,threshold,minfod)
%  Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

Lf=length(FOD);
fodmax=max(FOD);
if fodmax >0
FODnorm=FOD/sum(FOD); 
fmax=max(FODnorm);
fmin=minfod;
M=threshold*fmax;
Y=unidrnd(Lf);
U = rand(1);
flag=FODnorm(Y)/M;
test=1;
dirmax=find(FODnorm==fmax);
direction=dirmax(1);
while(U>flag || FOD(Y)<fmin)
    test=test+1;
    if test>100
      dirmax=find(FODnorm==fmax);
      direction=dirmax(1);
      break
    end
Y=unidrnd(Lf);
U = rand(1);
flag=FODnorm(Y)/M;
direction=Y; 
end
else
    direction=unidrnd(Lf);
end

end