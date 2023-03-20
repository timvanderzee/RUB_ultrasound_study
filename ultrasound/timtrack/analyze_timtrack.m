
clear all; close all; clc
% cd('C:\Users\timvd\Documents\RUB_ultrasound_study')
% Tmax    = readmatrix('max_torques.txt');
% Trest   = readmatrix('rest_torques.txt'); 

%% ramp conditions
force_conditions = {'slow','medium','fast','asym'};
image_qualities = {'low', 'high'};

for i = 1:2
for j = 1:4
    
    
% load
load([force_conditions{j},'_',image_qualities{i},'_timtrack.mat'])

faslen = faslen / 564 * 5; % pixels to cm
thickness = thickness / 564 * 5; % pixels to cm

figure(i)
color = get(gca,'colororder');

subplot(311);
plot(t, mean(phi,'omitnan'),'color',color(j,:),'linewidth',2); hold on
plot(t, mean(phi,'omitnan')-std(phi,1,'omitnan'),'-','color',color(j,:), 'linewidth',.5)
plot(t, mean(phi,'omitnan')+std(phi,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
xlabel('Time (s)'); ylabel('Angle (deg)'); title('Pennation angle'); box off; xlim([0 82])

subplot(312);
plot(t, mean(thickness,'omitnan'),'color',color(j,:),'linewidth',2); hold on
plot(t, mean(thickness,'omitnan')-std(thickness,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
plot(t, mean(thickness,'omitnan')+std(thickness,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
% yline(0,'k--'); yline(50,'k--')
xlabel('Time (s)'); ylabel('Thickness (cm)'); title('Muscle thickness'); box off; xlim([0 82])

subplot(313);
plot(t, mean(faslen,'omitnan'),'color',color(j,:),'linewidth',2); hold on
plot(t, mean(faslen,'omitnan')-std(faslen,1,'omitnan'),'-','color',color(j,:), 'linewidth',.5)
plot(t, mean(faslen,'omitnan')+std(faslen,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
xlabel('Time (s)'); ylabel('Length (cm)'); title('Fascicle length'); box off; xlim([0 82])

if i == 2
% drift

for k = 1:5
    coef_phi       = polyfit(t,  phi(k,:),1);
    drift_phi(k,j) = coef_phi(1);

    coef_len       = polyfit(t,  faslen(k,:),1);
    drift_len(k,j) = coef_len(1);

    % noise
    noise_phi(k, j) = mean(abs(diff(phi(k,:))));
    noise_len(k, j) = mean(abs(diff(faslen(k,:))));
end

end

end
end

%%
sine_conditions = {'sine_020','sine_1020'};

for j = 1:2
    
    figure(3)
    color = get(gca,'colororder');
    
    % load
    load([sine_conditions{j},'_timtrack.mat'])
    
    faslen = faslen / 564 * 5; % pixels to cm
thickness = thickness / 564 * 5; % pixels to cm
    
    subplot(311);
    plot(t, mean(phi,'omitnan'),'color',color(j,:),'linewidth',2); hold on
    plot(t, mean(phi,'omitnan')-std(phi,1,'omitnan'),'-','color',color(j,:), 'linewidth',.5)
    plot(t, mean(phi,'omitnan')+std(phi,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
    xlabel('Time (s)'); ylabel('Angle (deg)'); title('Pennation angle'); box off; xlim([0 82])

    subplot(312);
    plot(t, mean(thickness,'omitnan'),'color',color(j,:),'linewidth',2); hold on
    plot(t, mean(thickness,'omitnan')-std(thickness,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
    plot(t, mean(thickness,'omitnan')+std(thickness,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
    % yline(0,'k--'); yline(50,'k--')
    xlabel('Time (s)'); ylabel('Thickness (pixels)'); title('Muscle thickness'); box off; xlim([0 82])

    subplot(313);
    plot(t, mean(faslen,'omitnan'),'color',color(j,:),'linewidth',2); hold on
    plot(t, mean(faslen,'omitnan')-std(faslen,1,'omitnan'),'-','color',color(j,:), 'linewidth',.5)
    plot(t, mean(faslen,'omitnan')+std(faslen,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
    xlabel('Time (s)'); ylabel('Length (pixels)'); title('Fascicle length'); box off; xlim([0 82])
    
    
for k = 1:5
    coef_phi       = polyfit(t,  phi(k,:),1);
    drift_phi(k,j+4) = coef_phi(1);

    coef_len       = polyfit(t,  faslen(k,:),1);
    drift_len(k,j+4) = coef_len(1);

    % noise
    noise_phi(k, j+4) = mean(abs(diff(phi(k,:))));
    noise_len(k, j+4) = mean(abs(diff(faslen(k,:))));
end
end