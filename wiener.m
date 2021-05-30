function [h0, h1, h2] = wiener(a, b)
%
%      [h0, h1, h2] = wiener(a, b);
%
%       simulates a simple static non-linearity: au+bu^2 that follows a linear filter
%        Note that when b=0 the model is linear
%
%       the linear filter is: g(t) = a * exp(-t/T1) * sin(2*pi*f0*t) with T1=50 ms and f0 = 0.03 Hz/ms
%
%       E.g., try  [h0,h1,h2] = wiener(0,3); (no linearity)
%                  [h0,h1,h2] = wiener(1,0); (no nonlinearity)
%                  [h0,h1,h2] = wiener(1,3);
%
%      A couple of instances will require some tweaking with the time base
%      (to prevent too much contribution from random noise) 
%      and amplitude (to prevent spilling over of information to the next functional)
%
%      This is part of the exercise!
%
close all
if nargin<2 
    b=2;
end
if nargin<1
   a= 1;
end
  
gwnt = 0.5 * randn(1,4000);    % gaussian noise with power (close to) 0.5
tf = [0:1:500];          % time for the filter  unit: ms
time=[0:1:3999];         % time for the signal  unit: ms

T1=50; T2=8;              % time constants in ms
ht = exp(-tf/T1) .* sin(2 * pi * 0.03 * tf);   % impulse response of the linear system
s=sum((ht)); %similarly to h1 this is not yet normalised by this s
ht = ht/s;      % normalise the impulse response to make its integral = 1

P = std(gwnt)          % the power of the noise

figure(3)
plot(tf, ht, 'k-');
hold on

ut = conv(gwnt, ht);    % determine u(t): the output of the linear process
ut = ut(1:length(gwnt)); % only the first gwnt samples

figure(4)
plot(time, ut, 'r-');
hold on
yt = a*ut + b*ut.^2;    % the static non-linearity

plot(time, yt, 'k-');


%
%    now determine the zero-order Wiener functional
%
h0=mean(yt);
G0 = h0;       % the zero-order Wiener prediction
y0 = yt - G0;      % the Lee-Schetzen method: subtract h0 from the response

%  first order estimate: crosscorrelate
%
fi1 = phi_yx2(y0, gwnt, length(tf)-1);   %first-order kernel
h1=fi1/P;       % only positive tau! 
h1(200:end)=0; %listened to the hint, this reduced the noise a lot indeed
%although I think a closer approximation is obtained for 200:end being 0,
%as there is some information between 150 and 200 (in the shape of a peak)
s=sum((h1)); %get a sum over h1 of one if we take away the absolute value, but get worse results
h1 = h1/s; % normalise the integral to one (this is a "tweaking point.." )
dt=1;
disp(sum(h1)*dt);

m1=max(h1);
m2=max(ht);
m3=m1/m2;     
 
 %    there may be a scaling issue. You then need to compensate for this, e.g. by aligning
 %    the peaks... 
 %
 %    you may also want to truncate the h1 filter to zero after a certain
 %    time (in this case, e.g. after 150 ms), to get rid of the extra
 %    noise.
 %    You could play with this in your own simulations. Report about this!
 %    
 %    The peaks in figure 4 do seem to be alligned after
 %    truncating the h1 filter at 150 ms, the G1 is roughly equal to the ut
 %    signal, which is the output of the linear process (so that seems
 %    correct)

figure(3)
plot(tf, h1/m3, 'r-');      % the predicted (scaled) filter is close to the real filter
legend('g(tau)','h_1(tau)');
title('Predicted first-order kernel (h1) and the model linear filter (g)');

G1 = conv(gwnt, h1);      % 1e order Wiener prediction with the scaled kernel 
G1 = G1(1:length(gwnt));
figure(4)
plot(G1); %/8.5) %when ht is not normalised then we need to divide h1 by ~8.5 for them to have the same order
legend('u(t)','y(t)','G1');

y1 = yt - G0 - G1;    % residue, corrected for G0 and G1

%
%    second order kernel: 3rd-order cross-correlation
%
fi2=phi_yxx2(y1, gwnt, length(tf)-1);              
h2=fi2/(2*P^2);         % the 2nd order Wiener kernel: only positive tau1, tau2 
%h2(1:200,1:200)=0;
s = sum(sum((h2)));  % total integral of the kernel
h2 = h2/s;              % normalize the kernel to one   
disp(sum(sum(h2))*dt);
%%%%%%%%%%%%%%%%%%%
%  Note: also this scaling may have to be tweaked a little bit to optimize
%  your predictions! Play with this.
%%%%%%%%%%%%%%%%%%
G2 = convh2xx(h2,gwnt,P);      % 1e order Wiener prediction with the scaled kernel 
figure(10)
plot(G2);
hold on
plot(yt);
hold on
plot(G0+G1+G2);
legend('G2','y(t)','total')

disp(corrcoef((G0+G1+G2),yt));
disp(corrcoef((G2),yt));
disp(corrcoef((G1),yt));
figure(2)
imagesc(h2(1:100,1:100));    % Look parallel to the main diagonal!
axis xy
colorbar;
title('Second-order kernel')

figure(5)
Diagh2 = zeros(1,length(tf));
for n=1:length(tf)
    Diagh2(n) = h2(n,n);
end
 mx = max(Diagh2(1:150));
mx2 = max(h1(1:150));
plot(Diagh2(1:150)*mx2/mx,'r-','linewidth',2);
hold on
plot(h1(1:150), 'k-');
title('Diagonal of h2 and h1');

%
%  now you can also calculate G2 with convh2xx() using h2 and the input noise.
%
%  From here on you can use the two functionals to predict the system's
%  output. Try to tweak it such that you seem to get a decent result..
%  you may have to scale some factors for G0, G1 and G2. This is because
%  the noise is not 'perfect'
%
%  Write about this initial stage of your simulations in your report!
%

