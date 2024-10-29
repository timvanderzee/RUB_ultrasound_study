dt = 1/30;

load('exampleVideoTA_Verheul_2022_tracked_Q=0001.mat')

N = length(Fdat.Region.PEN);
t = 0:dt:((N-1)*dt);
pixpermm =  13.6;


%% Hybrid Track: settings
load('exampleVideoTA_trackingResults.mat','settings')

%% manual tracking
obs = {'TZ', 'PT', 'BR'};
ls = {'-','--',':'};

Fangle = nan(101,3);
Flength = nan(101,3);
Aangle = nan(101,3);

for i = 1:length(obs)
cd('C:\Users\u0167448\Documents\GitHub\RUB_ultrasound_study\data\ultrasound\TA_video\manual')
    load(['Manual_Tracking_',obs{i},'.mat'],'FasData','ApoData')
    
    [ts, id] = sort(t(FasData.digitizedFrames));
    
%     FasData.FAngle(FasData.FAngle<0) = FasData.FAngle(FasData.FAngle<0) + 180;
    
%     Fangle(1:length(id),i) = FasData.FAngle(id)';
%     Aangle(1:length(id),i) = ApoData.Angle(id)';
%     Flength(1:length(id),i) = FasData.FLength(id)'/pixpermm;

    % recalc stuff (because of different horizontal and vertical pixtomm)
    for j = 1:length(id)
        fas_dX = (FasData.pts{id(j)}(2,1) - FasData.pts{id(j)}(1,1)) / settings.horzmm;
        fas_dY = (FasData.pts{id(j)}(2,2) - FasData.pts{id(j)}(1,2)) / settings.vertmm;

        apo_dX = (ApoData.pts{id(j)}(2,1) - ApoData.pts{id(j)}(1,1)) / settings.horzmm;
        apo_dY = (ApoData.pts{id(j)}(2,2) - ApoData.pts{id(j)}(1,2)) / settings.vertmm;

        Flength(j,i) = sqrt(fas_dX^2 + fas_dY^2);
        Fangle(j,i) = atan2d(fas_dY, fas_dX);
        Aangle(j,i) = atan2d(apo_dY, apo_dX);
    end

    Fangle(Fangle<0) = Fangle(Fangle<0) + 180;
%     fasA(fasA<0) = fasA(fasA<0) + 180;
end

Pangle = Fangle - Aangle;

j = 1:2:20;

figure(8)
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
load('exampleVideoTA_Verheul_2022_tracked_Q=0001.mat')
ANG = Fdat.Region.ANG;
PEN = Fdat.Region.PEN;
FL = Fdat.Region.FL;

subplot(321)
plot(t, FL,'linewidth',2,'color',color(4,:)); hold on
title('UltraTimTrack: fascicle length')

subplot(322)
plot(t, PEN,'linewidth',2,'color',color(4,:)); hold on
title('UltraTimTrack: fascicle angle')

RMSE(1,:) = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]

%% UTT v2
load('exampleVideoTA_tracked_Q=0001_shifted.mat')
ANG = Fdat.Region.ANG;
PEN = Fdat.Region.PEN;
FL = Fdat.Region.FL;

subplot(321)
plot(t, FL,':','linewidth',1,'color',color(4,:)); hold on
title('UltraTimTrack: fascicle length')

subplot(322)
plot(t, PEN,':','linewidth',1,'color',color(4,:)); hold on
title('UltraTimTrack: fascicle angle')

RMSE(2,:) = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]


%% Hybrid Track
load('exampleVideoTA_trackingResults.mat', 'Hybrid')

N = 601;
clear FL PEN
for i = 1:601
    FL(i) = Hybrid(i).length;
    PEN(i) = Hybrid(i).angle_aac;
end

t = 0:dt:((N-1)*dt);

subplot(324)
plot(t, PEN,'linewidth',2,'color',color(6,:))
title('HybridTrack: fascicle angle')

subplot(323)
plot(t, FL,'linewidth',2,'color',color(6,:))
title('HybridTrack: fascicle length')

RMSE(3,:) = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]


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
title('DL Track: fascicle angle')

subplot(325)
plot(t, FL,'color',color(5,:))
title('DL Track: fascicle length')

RMSE(4,:) = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]

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
figure(8)
set(gcf,'units','normalized','position',[.1 .1 .35 .6])

