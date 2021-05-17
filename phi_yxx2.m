%
%     function toe calculate the third-order cross-correlation function between output y(t) 
%     and GWN input x(t):
%
%        phi_yxx(s1, s2) = INT y(t) x(t-s1) x(t-s2) dt
%
%        call:      A = phi_yxx(y, x, M)
%        y and x are arrays length N, the cross-correlation function has size MxM
%        (based on the  system's 'memory' (default: M=300 time steps)
%         A is a (M+1 bij M+1) matrix
% 
%    this function is needed for the Lee-Schetzen cross-correlation technique 
%     for non-linear systems identification. It ensures that only s1,s2>=0
%     are given (causality)
%
%    
function FI = phi_yxx2(y,x, M)

if nargin<3
     M=300;
end

N = length(x);
FI=zeros(M+1,M+1);
h=waitbar(0,'Please, be patient...');

% first make a matrix of the stimuli: (M+1)xN

xdelay = zeros(M+1, N);
xdelay(1,:) = x;
for s=2:M+1  
   xshift=circshift(x,[0 s]); 
   for k=1:s
       xshift(k) = 0;
   end
   xdelay(s,:) = xshift;   
end

      
  for t2=1:M+1
    for t1=1:M+1
          FI(t1, t2) = sum(y .* xdelay(t2,:) .* xdelay(t1,:));       % no normalisation needed; we assume that dt = 1       
          perc = 100*t1*t2/(M+1)^2; 
          if abs(perc - floor(perc))<0.01
             waitbar(perc/100, h);
          end
        
    end  % t1  rows
  end    % t2  columns

  % integral = sum(sum(abs(FI)))
  % FI = FI/integral;
   
return

