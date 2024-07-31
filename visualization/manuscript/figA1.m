clear all; close all; clc

dt = 1/30;

Qs = [(.0001:.0001:.001) .002 .003 .004 .005 .01 .1 1 10 100 1000];

N = 601;
t = 0:dt:((N-1)*dt);
pixpermm =  13.6;


%% manual tracking
obs = {'TZ', 'PT', 'BR'};
ls = {'-','--',':'};

Fangle = nan(101,3);
Flength = nan(101,3);
Aangle = nan(101,3);

for i = 1:length(obs)
    load(['Manual_Tracking_',obs{i},'.mat'])
    
    [ts, id] = sort(t(FasData.digitizedFrames));
    
    FasData.FAngle(FasData.FAngle<0) = FasData.FAngle(FasData.FAngle<0) + 180;
    
    Fangle(1:length(id),i) = FasData.FAngle(id)';
    Aangle(1:length(id),i) = ApoData.Angle(id)';
    Flength(1:length(id),i) = FasData.FLength(id)'/pixpermm;
end

Pangle = Fangle - Aangle;

for i = 1:length(Qs)
    if Qs(i) < 1
        f = ['%.',num2str(-floor(log10(Qs(i)))),'f'];
    else
        f = '%.0f';
    end
    filename = ['exampleVideoTA_Verheul_2022_tracked_Q=',strrep(num2str(Qs(i), f),'.',''), '.mat'];
    
    load(filename)
    % load('exampleVideoTA_Verheul_2022.mat')


    ANG = Fdat.Region.ANG;
    PEN = Fdat.Region.PEN;
    FL = Fdat.Region.FL;


    RMSE(i,:) = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]
end


%% HybridTrack
%% Hybrid Track
% cd('C:\Users\u0167448\Documents\hybrid-muscle-tracking-main\TrackingResults')
load('exampleVideoTA_trackingResults.mat')

N = 601;
clear FL PEN
for i = 1:601
    FL(i) = Hybrid(i).length;
    PEN(i) = Hybrid(i).angle;
end

t = 0:dt:((N-1)*dt);


HRMSE = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]

%%
if ishandle(1), close(1); end
figure(1)
subplot(121);
semilogx(Qs, RMSE(:,1),'linewidth',2); hold on
plot([min(Qs) max(Qs)], [HRMSE(1) HRMSE(1)],'k-')   
ylabel('RMSD (mm)')
ylim([0 6])
title('Fascicle length')

subplot(122);
semilogx(Qs, RMSE(:,2),'linewidth',2); hold on
plot([min(Qs) max(Qs)], [HRMSE(2) HRMSE(2)],'k-')
box off
ylim([0 1])
ylabel('RMSD (deg)')
title('Fascicle angle')

for i = 1:2
    subplot(1,2,i)
    xlabel('Process noise covariance parameter c')
    box off
    xlim([min(Qs) max(Qs)])
end