function [h0, h1, h2] = wiener(a, b,trunc_start, scaling_factor, mean_compensation)
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
if nargin<5
   mean_compensation= 0; %this factor is added for the case where the predictions seem to be off by a constant, which is the case for a large non-linear part
else
   mean_compensation= 1; 
end

gwnt = 0.5 * randn(1,4000);    % gaussian noise with power (close to) 0.5
%now check the gwnt signal for its properties
close all;
[f,xi] = ksdensity(gwnt,'Function','pdf'); 
sigma=std(gwnt);
mu=mean(gwnt);
x=linspace(-5,5,100);
g=(1/(sqrt(2*pi)*sigma))*exp(-((x-mu).^2)/(2*sigma^2));
%g is the theoretical gaussian, given the mean and std of gwnt
figure(100)
plot(xi,f);
hold on
plot(x,g);
title('Amplitude distribution')
legend('Measured','Theoretical')
xlabel('x')
ylabel('P(x)')

fs = 1000; % sample frequency (Hz), is 1000 Hz as dt is 1ms, as resolution of signals is 1ms

figure(101)
pspectrum(gwnt,fs);
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
title('Power spectrum')
%seems flat, especially when compared to the spectrum of eg. a sinusoid

deviation=3999;
a_corr_gwn_xx=phi_yx2(gwnt,gwnt,deviation);
figure(102)
xlim([0 deviation])
plot(a_corr_gwn_xx)
title('Full autocorrelation')
figure(103)
xlim([0 200])
plot(a_corr_gwn_xx(1:200))
title('First 200 time lags of the autocorrelation')
ylabel('\phi_{xx}')
xlabel('\tau')
%end of property check

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
title('The linear and non-linear parts of the signal and the second order prediction')
hold on
yt = a*ut + b*ut.^2;    % the static non-linearity
plot(time, yt, 'k-');


%
%    now determine the zero-order Wiener functional
%
h0=mean(yt);
G0 = h0;       % the zero-order Wiener prediction
if mean_compensation==0
    y0 = yt - G0;      % the Lee-Schetzen method: subtract h0 from the response
else 
    y0 = yt;
end
%  first order estimate: crosscorrelate
%
fi1 = phi_yx2(y0, gwnt, length(tf)-1);   %first-order kernel
h1=fi1/P;       % only positive tau! 
h1(trunc_start:end)=0; %listened to the hint, this reduced the noise a lot indeed
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
if mean_compensation==0
    y1 = yt - G0 - G1;    % residue, corrected for G0 and G1
else 
    y1 = yt - G1;    % residue, corrected for G0 and G1
end

%
%    second order kernel: 3rd-order cross-correlation
%
display(length(tf)-1)
fi2=phi_yxx2(y1, gwnt, 1000); %originally length(tf)-1 as third parameter              
h2=fi2/(2*P^2);         % the 2nd order Wiener kernel: only positive tau1, tau2 
%h2(150:end,150:end)=0; %this clearly improves performance here
s = sum(sum((h2)));  % total integral of the kernel
h2 = h2/s;              % normalize the kernel to one   
disp(sum(sum(h2))*dt);
%%%%%%%%%%%%%%%%%%%
%  Note: also this scaling may have to be tweaked a little bit to optimize
%  your predictions! Play with this.
%%%%%%%%%%%%%%%%%%
G2 = convh2xx(h2,gwnt,P,1);      % 1e order Wiener prediction with the scaled kernel 
figure(4)
plot(G0+G1+G2); %/8.5) %when ht is not normalised then we need to divide h1 by ~8.5 for them to have the same order
legend('u(t)','y(t)','Total prediction (always includes G0)');

figure(29)
plot(G1);
title('The linear part of the signal, y0, and the prediction from the first order kernel, G1')
hold on
plot(y0);
legend('G1','y0')

figure(10)
plot(G2);
title('The non-linear part of the signal, y1, and the prediction from the second order kernel, G2')
hold on
plot(y1);
legend('G2','y1')
figure(11)
plot(yt);
title('The signal, and the second order prediction')
hold on
if mean_compensation==0
    plot((G0+G1+G2)/scaling_factor);
else 
    plot((G1+G2)/scaling_factor);
end
legend('y(t)','second order')
% disp("corr second order and yt");
% disp(corrcoef((G0+G1+G2),yt));
% disp("corr first order and yt");
% disp(corrcoef((G0+G1),yt));
if mean_compensation==0
    disp("corr second order and yt");
    disp(corrcoef((G0+G1+G2),yt));
    disp("corr first order and yt");
    disp(corrcoef((G0+G1),yt));
else 
    disp("corr second order (without G0) and yt");
    disp(corrcoef((G1+G2),yt));
    disp("corr first order (without G0) and yt");
    disp(corrcoef((G1),yt));
end
disp("corr G2 and y1");
disp(corrcoef((G2),y1));
disp("corr G1 and y0");
disp(corrcoef((G1),y0));
figure(2)
imagesc(h2(1:100,1:100));    % Look parallel to the main diagonal!
axis xy
colorbar;
title('Second-order kernel')

disp("Mean and standard deviation of the difference between the signal and the second order estimation");
disp(mean(abs(yt-((G0+G1+G2)/scaling_factor))));
disp(std(abs(yt-((G0+G1+G2)/scaling_factor))));
disp("Mean and standard deviation of the difference between the signal and the first order estimation");
disp(mean(abs(yt-((G0+G1)/scaling_factor))));
disp(std(abs(yt-((G0+G1)/scaling_factor))));
disp("Mean and standard deviation of the difference between the signal and G1+G2");
disp(mean(abs(yt-((G1+G2)/scaling_factor))));
disp(std(abs(yt-((G1+G2)/scaling_factor))));
disp("Mean and standard deviation of the difference between the signal and G1");
disp(mean(abs(yt-((G1)/scaling_factor))));
disp(std(abs(yt-((G1)/scaling_factor))));

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

