clear all; close all; clc
% cd('C:\Users\timvd\Documents\RUB_ultrasound_study')
Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('RampTarget.mat','tnew','rampTarget')

%% ramp conditions
force_conditions = {'slow','medium','fast','asym'};
image_qualities = {'low', 'high'};

mTrel = nan(4,2);
sTrel = nan(4,2);

for i = 1:2
for j = 1:4
    
% load
load([force_conditions{j},'_',image_qualities{i},'_summary.mat'])

% subtract rest torque, divide by MVC torque, multiply by 100%
Trel = (torque - Trest(:)) ./ Tmax(:) * 100;

mTrel(j,i) = mean(mean(Trel,2,'omitnan'));
sTrel(j,i) = std(mean(Trel,2,'omitnan'));

figure(i)
color = get(gca,'colororder');

subplot(211);
plot(tnew, mean(torque,'omitnan'),'color',color(j,:),'linewidth',2); hold on
plot(tnew, mean(torque,'omitnan')-std(torque,1,'omitnan'),'-','color',color(j,:), 'linewidth',.5)
plot(tnew, mean(torque,'omitnan')+std(torque,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
xlabel('Time (s)'); ylabel('Torque (N-m)'); title('Absolute torque'); box off; xlim([0 82])

subplot(212);
plot(tnew, rampTarget(j,:),'k-'); hold on

plot(tnew, mean(Trel,'omitnan'),'color',color(j,:),'linewidth',2); hold on
plot(tnew, mean(Trel,'omitnan')-std(Trel,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
plot(tnew, mean(Trel,'omitnan')+std(Trel,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
yline(0,'k--'); yline(50,'k--')
xlabel('Time (s)'); ylabel('Torque (%MVC)'); title('Relative torque'); box off; xlim([0 82])

RMSEr(j,i) = rms(mean(Trel(:,2:end-1),'omitnan') - rampTarget(j,2:end-1));



end
end

%% sine conditions
sine_conditions = {'sine_020','sine_1020'};

sineTarget = zeros(2, length(tnew));

O = [10; 15];
A = [10; 5];

for j = 1:2
    sineTarget(j,:) = A(j)*cos(2*pi*1.5*tnew) + O(j);
end

mTrel = nan(1,2);
sTrel = nan(1,2);


for j = 1:2
    
    % load
    load([sine_conditions{j},'_summary.mat'])
    
    % subtract rest torque, divide by MVC torque, multiply by 100%
    Trel = (torque - Trest(:)) ./ Tmax(:) * 100;
    
    mTrel(j) = mean(mean(Trel,2,'omitnan'));
    sTrel(j) = std(mean(Trel,2,'omitnan'));

    figure(3)
    color = get(gca,'colororder');

    subplot(211);
    plot(tnew, mean(torque,'omitnan'),'color',color(j,:),'linewidth',2); hold on
    plot(tnew, mean(torque,'omitnan')-std(torque,1,'omitnan'),'-','color',color(j,:), 'linewidth',.5)
    plot(tnew, mean(torque,'omitnan')+std(torque,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
    xlabel('Time (s)'); ylabel('Torque (N-m)'); title('Absolute torque'); box off; xlim([0 82])

    subplot(212);

    plot(tnew, sineTarget(j,:),'k-'); hold on
    
    plot(tnew, mean(Trel,'omitnan'),'color',color(j,:),'linewidth',2); hold on
    plot(tnew, mean(Trel,'omitnan')-std(Trel,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
    plot(tnew, mean(Trel,'omitnan')+std(Trel,1,'omitnan'),'-','color',color(j,:),'linewidth',.5)
    yline(0,'k--'); yline(10,'k--'); yline(20,'k--')
    xlabel('Time (s)'); ylabel('Torque (%MVC)'); title('Relative torque'); box off; xlim([0 82])
    
    RMSEs(j) = rms(mean(Trel(:,2:end-1),'omitnan') - sineTarget(j,2:end-1));
end

mean([RMSEr; RMSEs], 'all')
std([RMSEr; RMSEs], 1, 'all')

%% save figures
for i = 1:3
    figure(i); 
    saveas(gcf, ['fig',num2str(i)], 'fig')
end