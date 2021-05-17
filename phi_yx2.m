%
%     function to calculate the  cross correlation between output y(t) 
%     and GWN input x(t) :
%
%      phi_yx(s) = INT y(t) x(t-s) dt
%
%   call: FI = phi_yx(y, x, M)
% 
%   used for the Lee-Schetzen method
%    of non-linear systems identification.
%    
%     note that s>=0 (causality)
%
function FI = phi_yx2(y,x, M)

if nargin<3 
    M=300;
end

N = length(x);
FI=zeros(1, M+1);

%
%  make a delay matrix of x
%

xdelay = zeros(M+1, N);
xdelay(1,:) = x;
for s=2:M+1  
   xshift=circshift(x,[0 s]); 
   for k=1:s
       xshift(k) = 0;
   end
   xdelay(s,:) = xshift;
end

figure(5)
clf
imagesc(xdelay);
axis xy
title('Delayed input signal')

    for t1=1:M+1   
       FI(t1) = sum(y .* xdelay(t1,:)); 
    end

   % integral = (sum(abs(FI)));
   % FI = FI/integral;   % normalize to 1.0
    
return