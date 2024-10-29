dt = 1/81.6060;

load('calf_raise_tracked_Q=0001.mat')

N = length(Fdat.Region.PEN);
t = 0:dt:((N-1)*dt);
pixpermm = 508 / 65;

ANG = Fdat.Region.ANG;
PEN_UTT = Fdat.Region.PEN;
FL_UTT = Fdat.Region.FL;

%% manual tracking
obs = {'TZ','PT','BR'};
ls = {'-','--',':'};

Fangle = nan(84,3);
Flength = nan(84,3);
Aangle = nan(84,3);

for i = 1:length(obs)
cd('C:\Users\u0167448\Documents\GitHub\RUB_ultrasound_study\data\ultrasound\calf_raise\manual')
    load(['Manual_Tracking_',obs{i},'.mat'],'ApoData','FasData')
    
    [ts, id] = sort(t(FasData.digitizedFrames));
    
    FasData.FAngle(FasData.FAngle<0) = FasData.FAngle(FasData.FAngle<0) + 180;
    
    Fangle(1:length(id),i) = FasData.FAngle(id)';
    Aangle(1:length(id),i) = ApoData.Angle(id)';
    Flength(1:length(id),i) = FasData.FLength(id)'/pixpermm;
end

Pangle = Fangle - Aangle;

j = 1:2:20;

figure(9)
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

color = get(gca,'colororder');


%% UTT
subplot(321)
plot(t, FL_UTT,'linewidth',2,'color',color(4,:)); hold on
title('UltraTimTrack: fascicle length')

subplot(322)
plot(t, PEN_UTT,'linewidth',2,'color',color(4,:)); hold on
title('UltraTimTrack: fascicle angle')

RMSE(1,:) = [sqrt(mean((FL_UTT(FasData.digitizedFrames)' - mean(Flength,2,'omitnan')).^2,'omitnan')) sqrt(mean((PEN_UTT(FasData.digitizedFrames)' - mean(Pangle,2,'omitnan')).^2,'omitnan'))]


% legend('Tim','Paolo','Brent','location','best')
%% Hybrid Track
load('calf_raise_crop_trackingResults.mat', 'Hybrid')

N = length(Hybrid);
clear FL PEN
for i = 1:N
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

RMSE(2,:) = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]


%% DL_track
Ldata = readmatrix('calf_raise.xlsx', 'Sheet','Fasc_length');
Pdata = readmatrix('calf_raise.xlsx', 'Sheet','Pennation');
Ldata = Ldata(2:end,2:end);
Pdata = Pdata(2:end,2:end);

N = size(Ldata,1);
t = 0:dt:((N-1)*dt);

PEN_DLT = nan(1,size(Ldata,1));
FL_DLT = nan(1,size(Ldata,1));

for i = 1:size(Ldata,1)
    FL_DLT(i) = median(Ldata(i,:),'omitnan') / pixpermm;
    PEN_DLT(i) = median(Pdata(i,:),'omitnan');
end


subplot(326)
plot(t, PEN_DLT,'color',color(5,:))
title('DL Track: fascicle angle')

subplot(325)
plot(t, FL_DLT,'color',color(5,:))
title('DL Track: fascicle length')

RMSE(3,:) = [sqrt(mean((FL_DLT(FasData.digitizedFrames)' - mean(Flength,2,'omitnan')).^2,'omitnan')) sqrt(mean((PEN_DLT(FasData.digitizedFrames)' - mean(Pangle,2,'omitnan')).^2,'omitnan'))]

%% make nice
for i = 1:6
    subplot(3,2,i)
    xlim([0 max(t)]); 
    hold on
    box off

    axis tight
    
    if mod(i,2)
        ylabel('Fascicle length (mm)')
        ylim([25 60])
    else
        ylabel('Fascicle angle (deg)')
        ylim([13 32])
    end
    
    if i < 5
        set(gca,'XTicklabel',[])
        xlabel('')
    else
        xlabel('Time (s)'); 
    end
end



%% 
figure(9)
set(gcf,'units','normalized','position',[.1 .1 .35 .6])

