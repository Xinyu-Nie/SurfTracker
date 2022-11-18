function B = matrixB(A,maxOrder)
%% Input: maxOrder = 0,2,4,6,8
%         A = matrix of size n x 3  
%% P = P^m_l(cos(theta))
% P = (P_1(cos(theta_1)).....P_R(cos(theta_1))
%      P_1(cos(theta_2)).....P_R(cos(theta_2))
%
%      .....................
%      P_1(cos(theta_n)).....P_R(cos(theta_n)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = (...sin(2*phi_1) sin(phi_1) 1 sin(phi_1) sin(2*phi_1)...
%  ...sin(2*phi_2) sin(phi_2) 1 sin(phi_2) sin(2*phi_2)...
%      ..........................................
%  ...sin(2*phi_n) sin(phi_n) 1 sin(phi_n) sin(2*phi_n)...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rows,cols] = size(A);

thetavec = zeros(rows,1);
thetavec = atan2(sqrt(A(:,1).^2+A(:,2).^2),A(:,3));

phivec = zeros(rows,1);
phivec = atan2(A(:,2),A(:,1));

order = (0:2:maxOrder)';
R = zeros(length(order),1);
R = (order+1).*(order+2)/2;

P = zeros(rows,R(end));
T = zeros(rows,R(end));
scale = zeros(1,R(end));

P(:,1) = legendre(0,cos(thetavec));
T(:,1) = ones(rows,1);
scale(1) = sqrt(1/(4.0*pi));

indx =0;
for i = 2 : length(order)
  temp = legendre(order(i),cos(thetavec)); temp = temp';
  P(:,(R(i-1)+1):R(i)) = [temp(:,order(i)+1:-1:2) temp];
  T(:,(R(i-1)+1):R(i)) = [sin(repmat(phivec,1,order(i)).*repmat(order(i):-1:1,rows,1)) ones(rows,1) cos(repmat(phivec,1,order(i)).*repmat(1:1:order(i),rows,1))];
  indx =0;
  for m = 1:order(i)
      indx = indx+1;
      tempv(indx) = sqrt(((2*order(i)+1)/(2*pi))/prod((order(i)-m+1):1:(order(i)+m)));
  end
  scale(R(i-1)+1:R(i)) = [tempv(end:-1:1) sqrt((2*order(i)+1)/(4.0*pi)) tempv(1:1:end)];
end
B = repmat(scale,rows,1).*P.*T;

    




