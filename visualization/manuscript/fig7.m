clear all; close all; clc

dt = 1/30;

load('exampleVideoTA_Verheul_2022_tracked_Q=0001.mat')
% load('exampleVideoTA_Verheul_2022.mat')

N = length(Fdat.Region.PEN);
t = 0:dt:((N-1)*dt);
pixpermm =  13.6;

ANG = Fdat.Region.ANG;
PEN = Fdat.Region.PEN;
FL = Fdat.Region.FL;

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

j = 1:2:20;

for i = 1:3
subplot(3,2,j(i))
plot(ts, mean(Flength,2),'color',[.5 .5 .5],'linewidth',2); hold on
plot(ts, max(Flength,[],2),'color',[.5 .5 .5],'linewidth',1)
plot(ts, min(Flength,[],2),'color',[.5 .5 .5],'linewidth',1)

subplot(3,2,j(i)+1)
plot(ts, mean(Pangle,2),'color',[.5 .5 .5],'linewidth',2); hold on
plot(ts, max(Pangle,[],2),'color',[.5 .5 .5],'linewidth',1)
plot(ts, min(Pangle,[],2),'color',[.5 .5 .5],'linewidth',1)
end

color = get(gca,'colororder')

%% UTT
figure(1)
subplot(321)
plot(t, FL,'linewidth',2,'color',color(4,:)); hold on

subplot(322)
plot(t, PEN,'linewidth',2,'color',color(4,:)); hold on

RMSE = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]


%% Hybrid Track
% cd('C:\Users\u0167448\Documents\hybrid-muscle-tracking-main\TrackingResults')
load('exampleVideoTA_trackingResults.mat')

N = 601;
clear FL PEN
for i = 1:601
    FL(i) = Hybrid(i).length;
    PEN(i) = Hybrid(i).angle_aac;
end

t = 0:dt:((N-1)*dt);

subplot(324)
plot(t, PEN,'linewidth',2,'color',color(6,:))

subplot(323)
plot(t, FL,'linewidth',2,'color',color(6,:))

RMSE = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]


%% DL_track
Ldata = readmatrix('exampleVideoTA_Verheul_2022_reduced.xlsx', 'Sheet','Fasc_length');
Pdata = readmatrix('exampleVideoTA_Verheul_2022_reduced.xlsx', 'Sheet','Pennation');
Ldata = Ldata(2:end,2:end);
Pdata = Pdata(2:end,2:end);

N = size(Ldata,1);
t = 0:dt:((N-1)*dt);

PEN = nan(1,size(Ldata,1));
FL = nan(1,size(Ldata,1));

for i = 1:size(Ldata,1)
    FL(i) = median(Ldata(i,:),'omitnan') / pixpermm;
    PEN(i) = median(Pdata(i,:),'omitnan');
end

FasData.digitizedFrames(end) = 600;

subplot(326)
plot(t, PEN,'color',color(5,:))

subplot(325)
plot(t, FL,'color',color(5,:))

RMSE = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]

%% make nice
for i = 1:6
    subplot(3,2,i)
    xlim([0 max(t)]); 
    hold on
    box off

    axis tight
    
    if mod(i,2)
        ylabel('Fascicle length (mm)')
    else
        ylabel('Fascicle angle (deg)')
    end
    
    if i < 5
        set(gca,'XTicklabel',[])
        xlabel('')
    else
        xlabel('Time (s)'); 
    end
end



%% 
figure(1)
set(gcf,'units','normalized','position',[.1 .1 .35 .6])


%%
cd('C:\Users\u0167448\OneDrive\8. Ultrasound comparison - TBD\figures\svg\output from MATLAB')
fig = gcf;
print(fig,'MyPlot.svg','-dsvg');

