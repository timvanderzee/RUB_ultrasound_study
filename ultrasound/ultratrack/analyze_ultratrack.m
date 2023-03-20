clear all; close all; clc
% cd('C:\Users\timvd\Documents\RUB_ultrasound_study')
% Tmax    = readmatrix('max_torques.txt');
% Trest   = readmatrix('rest_torques.txt'); 

%% ramp conditions
force_conditions = {'slow','medium','fast','asym'};
image_qualities = 'high';

cd('C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\ultratrack\all_subjects')

for j = 1:4
    
    
% load
load([force_conditions{j},'_',image_qualities,'_ultratrack.mat'],'faslen','phi')

% faslen = Fdat.Region.FL / 10; % pixels to cm
% phi = Fdat.Region.PEN * 180/pi; % pixels to cm

figure(1)
color = get(gca,'colororder');

t = linspace(0,82,2667);

subplot(211);
plot(t, mean(phi,'omitnan'),'color',color(j,:),'linewidth',2); hold on
plot(t, mean(phi,'omitnan')-std(phi,1,'omitnan'),'-','color',color(j,:), 'linewidth',.5)
plot(t, mean(phi,'omitnan')+std(phi,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
xlabel('Time (s)'); ylabel('Angle (deg)'); title('Pennation angle'); box off; xlim([0 82])

subplot(212);
plot(t, mean(faslen,'omitnan'),'color',color(j,:),'linewidth',2); hold on
plot(t, mean(faslen,'omitnan')-std(faslen,1,'omitnan'),'-','color',color(j,:), 'linewidth',.5)
plot(t, mean(faslen,'omitnan')+std(faslen,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
xlabel('Time (s)'); ylabel('Length (cm)'); title('Fascicle length'); box off; xlim([0 82])

% drift

% for k = 1
% coef_phi       = polyfit(t,  phi(k,:),1);
% drift_phi(k,j) = coef_phi(1);
% 
% coef_len       = polyfit(t,  faslen(k,:),1);
% drift_len(k,j) = coef_len(1);
% 
% % noise
% noise_phi(k, j) = mean(abs(diff(phi(k,:))));
% noise_len(k, j) = mean(abs(diff(faslen(k,:))));
% end


end

%% sine conditions
sine_conditions = {'sine_020','sine_1020'};
cd('C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\ultratrack\all_subjects')

for j = 1:2
    
    
% load
load([sine_conditions{j},'_ultratrack.mat'],'faslen','phi')

figure(2)
color = get(gca,'colororder');

t = linspace(0,82,2667);

subplot(211);
plot(t, mean(phi,'omitnan'),'color',color(j,:),'linewidth',2); hold on
plot(t, mean(phi,'omitnan')-std(phi,1,'omitnan'),'-','color',color(j,:), 'linewidth',.5)
plot(t, mean(phi,'omitnan')+std(phi,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
xlabel('Time (s)'); ylabel('Angle (deg)'); title('Pennation angle'); box off; xlim([0 82])

subplot(212);
plot(t, mean(faslen,'omitnan'),'color',color(j,:),'linewidth',2); hold on
plot(t, mean(faslen,'omitnan')-std(faslen,1,'omitnan'),'-','color',color(j,:), 'linewidth',.5)
plot(t, mean(faslen,'omitnan')+std(faslen,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
xlabel('Time (s)'); ylabel('Length (cm)'); title('Fascicle length'); box off; xlim([0 82])

% drift

% for k = 1
% coef_phi       = polyfit(t,  phi(k,:),1);
% drift_phi(k,j+4) = coef_phi(1);
% 
% coef_len       = polyfit(t,  faslen(k,:),1);
% drift_len(k,j+4) = coef_len(1);
% 
% % noise
% noise_phi(k, j+4) = mean(abs(diff(phi(k,:))));
% noise_len(k, j+4) = mean(abs(diff(faslen(k,:))));
% end


end


