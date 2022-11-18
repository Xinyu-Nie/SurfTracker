function Pval = PolyJacobi(n,alfa,beta,x)
% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

Pval=0;
for s= 0:1:n
    Pval = Pval + nchoosek(n+alfa,n-s)*nchoosek(n+beta,s)*(((x-1)/2)^s)*(((x+1)/2)^(n-s));
end
end