% 
% Demonstrates relative performance of Kalman filter 
% and Rauch-Tung-Striebel smoother on random walk estimation 
% 
clear all; 
close all; 
N = 500;  % Number of samples of process used in simulations 
Q = .0001;  % Variance of random walk increments 
R = .1;    % Variance of sampling noise 
sigw  = sqrt(Q); % Standard deviations 
sigv  = sqrt(R); 
% 
% 
Pcorr      = 100;           % Covariance of initial uncertainty
Ppred = 100;
xuncor(1)    = 0; % Initial value of true process 
xpred(1)   = 0;             % Initial (predicted) value of true process 
sawtooth   = sqrt(Pcorr); 
ts         = 0; 
% 
% Forward pass: filter 
t = (0:(N-1));
U = sin(0.1 * t);
dU = 1*0.1*cos(0.1 * t);

for k=2:N
     
   % uncorrected signal
  xuncor(k)  = xuncor(k-1) + dU(k) + sigw; % Random walk 

  % a priori
  xpred(k) = xcorr(k-1) + dU(k) + sigw*randn; 
  Ppred(k) = Pcorr(k-1) + Q; 
   
   z(k)     = U(k) + sigv*randn;            % Noisy sample 
   K(k)        = Ppred(k)/(Ppred(k)+R); 
   Pcorr(k) = Ppred(k) - K(k)*Ppred(k); 
   
   % corrected state estimate
   xcorr(k) = xpred(k) + K(k)*(z(k) - xpred(k));  % Kalman filter estimate 

end
% 
% Backward pass: smooth 
% 
xsmooth = xcorr; 
for k=N-1:-1:1
   A(k)          = Pcorr(k)/Ppred(k+1); 
%    xsmooth(k) = xsmooth(k) + A(k)*(xsmooth(k+1) - xpred(k+1)); 
   
   xsmooth(k) = xcorr(k) + A(k)*(xsmooth(k+1) - xpred(k+1)); 
end

%%
figure(1)
plot(t, U, t,xuncor,t,z,t,xpred,t,xcorr,t,xsmooth); 
legend('True','Flow','Measured','Predicted','Corrected','Smoothed'); 
title('DEMO #7: Kalman Filter versus Rauch-Tung-Striebel Smoother'); 
xlabel('Discrete Time'); 
ylabel('Random Walk'); 

return

%%
if ishandle(2), close(2); end
figure(2)
plot(U, xpred,'.'); hold on
plot(U, xuncor,'.')
plot(U, xsmooth,'.')

ylim([-2 2])
%%
figure; 
semilogy(ts,sawtooth); 
xlabel('Discrete Time'); 
ylabel('RMS Estimation Uncertainty'); 
title('''Sawtooth Curve'''); 