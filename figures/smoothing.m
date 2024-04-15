clear all; close all; clc

cd('C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\3011\analyzed\mat')
filename = 'pas_005_01_analyzed_Q=001';

load([filename,'.mat']);

N = length(Fdat.Region.FL);

%% recreate a priori estimates
clc
x_apriori(1) = Fdat.Region.Fascicle.alpha{1};
x_apost(1)=  Fdat.Region.Fascicle.alpha{1};
alpha_prev = Fdat.Region.Fascicle.alpha{1};

Qval = 0.01;

P_apost = 0;
P_apriori = P_apost;

for i = 1:N
    X(i) = Fdat.Region.Fascicle.alpha{i};
end

for i = 2:N
    
    w = Fdat.Region.warp(:,:,i-1);

    fas_prev = [Fdat.Region.Fascicle.fas_x{i-1}' Fdat.Region.Fascicle.fas_y{i-1}'];
    
    fas_new = transformPointsForward(w, fas_prev);
    
    alpha_prev = X(i-1);
    
    % Estimate the change in fascicle angle from the change in points
    dalpha = abs(atan2d(diff(fas_new(:,2)), diff(fas_new(:,1)))) - abs(atan2d(diff(fas_prev(:,2)), diff(fas_prev(:,1))));
    alpha_new = alpha_prev + dalpha;

    x_apriori(i) = alpha_new;

    % A priori state estimate
    x_minus = alpha_new;

    dx = abs(dalpha);
    dx(dx<0.005) = 0;
    f.Q = getQ(Qval, dx);
    f.R = Fdat.R(1);

    % apriori estimate from optical flow
    f.x_minus = x_minus;

    % measurement from Hough transform
    f.y = Fdat.geofeatures(i).alpha;

    % previous state covariance
    f.P_prev = P_apost(i-1);

    % run kalman filter
    F = run_kalman_filter(f);
    
    P_apriori(i) = F.P_minus;
    
    % ascribe
    alpha = F.x_plus;
    Kgain(i) = F.K;
    P_apost(i) = F.P_plus;
    
    x_apost(i) = alpha;
    y(i) = Fdat.geofeatures(i).alpha;
end

%%
close all
figure(1)
plot(x_apriori); hold on
plot(x_apost)
plot(X,'--')
plot(y)

legend('apriori','aposteriori','X','y','location','best')
legend boxoff

%%
figure(2)
plot(Kgain); hold on
plot(Fdat.Region.Fascicle.K,'--')
%% smoothing

xsmooth = nan(1,N);
xsmooth(N) = Fdat.Region.Fascicle.alpha{N};
xnew(N) = Fdat.Region.Fascicle.alpha{N};

for i = (N-1):-1:1
    xprev = Fdat.Region.Fascicle.alpha{i+1};
    xnew(i) = Fdat.Region.Fascicle.alpha{i};
    
    Pprev = Fdat.Region.Fascicle.fas_p{i+1}(2);
    Pnew = Fdat.Region.Fascicle.fas_p{i}(2);
    
    C(i) = Pnew / Pprev;
    
    xsmooth(i) = xnew(i) + C(i) * (xsmooth(i+1) - xprev);
end

figure(1)
plot(xnew); hold on
plot(xsmooth,'--')


function [K] = run_kalman_filter(k)
% this assumes we already have the aposteriori state estimate (k.x_minus),
% the measurement (k.y) and the process- and measurement noise covariances (k.R and k.Qvalue)

% a posteriori variance estimate
K.P_minus = k.P_prev + k.Q;

% kalman xshiftcor
K.K = K.P_minus / (K.P_minus + k.R);

if (K.K < 0) || (K.K > 1)
    disp('warning: kalman gain outside 0-1 interval');

    K.K(K.K<0) = 0;
    K.K(K.K>1) = 1;
end

% estimated state
K.x_plus = k.x_minus + K.K * (k.y - k.x_minus);

% estimated variance
K.P_plus = (1-K.K) * K.P_minus;
end

function[Q] = getQ(Qval, dx)

% Optical flow error is proportional to flow
Q = Qval  * dx^2;
end
