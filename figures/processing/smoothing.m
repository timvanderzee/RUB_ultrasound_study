clear all; close all; clc

p = 2;
Qs = [nan, 0, 10.^(-4:0), 1000, inf];
i = 4;
Qval = Qs(i);
k = 3;

foldernames = {'3011', '0812', '1312','1612','1601','1701','1901a','1901b'};
mainfolder = 'C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
filenames = {'*pas_005*.mp4','*pas_30*.mp4','*pas_120*.mp4'};

% get names
foldername = foldernames{p};
cd([mainfolder foldername]);
files = dir(filenames{k});
vidname = files.name(1:end-4);

filename = [vidname,'_analyzed_Q=',strrep(num2str(Qval),'.','')];

cd([mainfolder foldername,'\analyzed\mat']);

if exist([filename,'.mat'],'file')
    load([filename,'.mat']);
else
    disp('Does not exist')
end

N = length(Fdat.Region.FL);

%% recreate a priori estimates
clc
x_apriori(1) = Fdat.Region.Fascicle.alpha{1};
x_apost(1)=  Fdat.Region.Fascicle.alpha{1};
alpha_prev = Fdat.Region.Fascicle.alpha{1};
xflow(1) = Fdat.Region.Fascicle.alpha{1};
y(1) = Fdat.Region.Fascicle.alpha{1};

xkal(1) = Fdat.Region.Fascicle.alpha{1};

P_apost = 0;
P_apriori = P_apost;

for i = 1:N
    X(i) = Fdat.Region.Fascicle.alpha{i};
end

for i = 2:N
    
    fas_prev = [Fdat.Region.Fascicle.fas_x{i-1}' Fdat.Region.Fascicle.fas_y{i-1}'];
    alpha_prev = Fdat.Region.Fascicle.alpha{i-1};
    
    % transform
    w = Fdat.Region.warp(:,:,i-1);
    fas_new = transformPointsForward(w, fas_prev);
    
    % Estimate the change in fascicle angle from the change in points
    dalpha = abs(atan2d(diff(fas_new(:,2)), diff(fas_new(:,1)))) - abs(atan2d(diff(fas_prev(:,2)), diff(fas_prev(:,1))));
    
    % a priori state estimate
    x_minus = alpha_prev + dalpha;
    
    % estimate covariances
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
    
    xflow(i) = xflow(i-1) + dalpha;
    x_apriori(i) = x_minus;
    
    x_apost(i) = alpha;
    y(i) = Fdat.geofeatures(i).alpha;
end


%% smoothing
xsmooth = nan(1,N);
xsmooth(N) = Fdat.Region.Fascicle.alpha{N};
xnew(N) = Fdat.Region.Fascicle.alpha{N};

for i = (N-1):-1:1
%     xprev = x_apost(i);
    xprev = x_apriori(i+1);
    xnew(i) = Fdat.Region.Fascicle.alpha{i};
    
    Pprev = Fdat.Region.Fascicle.fas_p{i+1}(2);
    Pnew = Fdat.Region.Fascicle.fas_p{i}(2);
    
    Pprev = P_apriori(i+1);
    Pnew = P_apost(i);
    
    C(i) = Pnew / Pprev;
%     C(C>1) = 1;
%     C(i) = 1;
    xsmooth(i) = x_apost(i) + C(i) * (xsmooth(i+1) - xprev);
end

%%
close all
figure(1)

subplot(2,4,1:4);
plot(x_apriori); hold on
plot(x_apost)
plot(X,'--')
plot(y)
plot(xsmooth,'--')
plot(xflow,':')

legend('apriori','aposteriori','X','y','xsmooth','xflow','location','best')
legend boxoff

% subplot(312)
% plot(C)


%% plot vs angle
% passive trials

Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('RampTarget.mat','tnew','rampTarget')
load('MVC_EMG.mat');

%% time series
force_conditions = {'pas_005', 'pas_30','pas_120'};

fs = 100;
Wn = 10 / (.5*fs);
[b,a] = butter(2, Wn);

% n = 8201;
n = N;


% angle_rs = nan(n, 8, 3);

    
% load
load([force_conditions{k},'_summary.mat'])
tus = 0:.03:max(tnew);

% down-sample
angle_rs = interp1(tnew, angle(p,:)', tus);



%%


N = min([length(X) length(angle_rs)]);

pen = [X(1:N); x_apost(1:N); xsmooth(1:N); xflow(1:N)];
titles = {'Original','a posteriori', 'smooth', 'flow'};

% figure(p+10)

for m = 1:size(pen,1)

    subplot(2,4,m+4)
   plot(angle_rs(1:N), pen(m,1:N),'.'); hold on
   title(titles{m})
%        ylim([15 30])

end




%%
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
