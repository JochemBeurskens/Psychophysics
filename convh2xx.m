%
%    2e order Wiener or Volterra functional
%
%  G2(t) = int int h2(s1, s2) x(t-s1) x(t-s2) ds1 ds2 - P * int h2(s,s) ds 
%   or
%  G2(t) = int int h2(s1, s2) x(t-s1) x(t-s2) ds1 ds2
%
%  with h2(s1, s2) the second-order homogeneous Volterra kernel
%
%           call:   G2 = convh2xx(h2, x, 1) for Volterra
%                   G2 = convh2xx(h2, x)    for Wiener
%
%     You use this function to make model predictions with the 2nd order kernel 
%
%     See as an example, also wiener.m
%
%     Note 1: you have to ensure that the normalisation is correct. This is not implemented in this function. 
%     Try to calculate first your own simple known model (with wiener.m) and reconstruct until  you make it 
%     as good as possible....
%
%     Note 2: if you want to use the 2nd order convolution for Volterra
%     prediction (not with GWN as an input, call the function with the 4th
%     input set to 1
%
%      check wiener.m for a simple example that you can use as inspiration
%
function G2 = convh2xx(h2, x, P, V_flag)

      if nargin<4 
          V_flag = 0;  % Wiener functional, not Volterra
      end
      
      N=length(x);
      NH = size(h2, 1);
      G2 = zeros(1,N);       
      M=size(h2,1);  % this is the memory length of the filter
     
%
%  make a delay matrix of x
%

xdelay = zeros(M, N);
xdelay(1,:) = x;
for s=2:M  
   xshift=circshift(x,[0 s]); 
   for k=1:s
       xshift(k) = 0;
   end
   xdelay(s,:) = xshift;  
end
 
offset = 0;
if V_flag == 0       % Wiener's offset term
   for k = 1:M
       offset = offset + P*h2(k, k);
   end
end

for t=1:N
  for t2=1:M
    for t1=1:M
       G2(t) = G2(t) + h2(t1, t2) * xdelay(t2,t) * xdelay(t1,t); 
    end
  end    
end
  
G2 = G2 - offset;  

return
