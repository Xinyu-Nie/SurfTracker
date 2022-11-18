function D = WignerDmatrixnew(n,beta)
% Contribution
%  Author: Xinyu Nie
%  Created: 2022/11/10
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu
len=2*n+1;
D=zeros(len,len);
for m=0:1:n
    for k=m:n
      Pbeta=PolyJacobi(n-k,k-m,k+m,cos(beta));
      D(n+1+m,n+1+k)=(1/(2^k))*sqrt(factorial(n-k)*factorial(n+k)/(factorial(n-m)*factorial(n+m)))*(sin(beta)^(k-m))*((1+cos(beta))^m)*Pbeta;
       if k>m
        D(n+1+k,n+1+m)=D(n+1+m,n+1+k)*((-1)^(k-m));
        D(n+1-k,n+1-m)=D(n+1+m,n+1+k);
       end
      D(n+1-m,n+1-k)=D(n+1+m,n+1+k)*((-1)^(k-m));
     if m>0
     Pbetar=jacobiP(n-k,k-m,k+m,cos(pi-beta));
      D(n+1+m,n+1-k)=((-1)^(n-m))*(1/(2^k))*sqrt(factorial(n-k)*factorial(n+k)/(factorial(n-m)*factorial(n+m)))*(sin(pi-beta)^(k-m))*((1+cos(pi-beta))^m)*Pbetar;
      if k>m
       D(n+1+k,n+1-m)=D(n+1+m,n+1-k);
       D(n+1-k,n+1+m)=((-1)^(-k-m))*D(n+1+m,n+1-k);
      end
      D(n+1-m,n+1+k)=((-1)^(-k-m))*D(n+1+m,n+1-k);
     end
    end
end
end