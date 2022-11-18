function U = SpYrotationnew (maxOrder,beta)
% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

N=(maxOrder+1)*(maxOrder+2)/2;
U=zeros(N,N);
U(1,1)=1;
t=1;
for order = 2:2:maxOrder
    len=2*order+1;
    D = WignerDmatrixnew(order,beta);
    U(t+1:t+len,t+1:t+len)=D;
    t=t+len;
end

end