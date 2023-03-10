clear all; close all; clc
% cd('C:\Users\timvd\Documents\RUB_ultrasound_study')
load('MVC_EMG.mat');
% Trest   = readmatrix('rest_EMGs.txt'); 

%% ramp conditions
force_conditions = {'slow','medium','fast','asym'};
image_qualities = {'low', 'high'};

for i = 1:2
for j = 1:4
    
% load
load([force_conditions{j},'_',image_qualities{i},'_EMG.mat'])

% subtract rest torque, divide by MVC torque, multiply by 100%
MGrel = EMG.MG ./ max(MVC.MG,[],2) * 100;
LGrel = EMG.LG ./ max(MVC.LG,[],2) * 100;
TArel = EMG.TA ./ max(MVC.TA,[],2) * 100;
SOrel = EMG.SO ./ max(MVC.SO,[],2) * 100;

figure(i)
color = get(gca,'colororder');

figure(i); plotEMG(tnew, EMG.MG , MGrel, color(j,:))
yline(0,'k--'); yline(50,'k--');
set(gcf,'name','MG')

figure(i+2); plotEMG(tnew, EMG.LG , LGrel, color(j,:))
yline(0,'k--'); yline(50,'k--');
 set(gcf,'name','LG')

figure(i+4); plotEMG(tnew, EMG.SO , SOrel, color(j,:))
yline(0,'k--'); yline(50,'k--');
 set(gcf,'name','SO')

figure(i+6); plotEMG(tnew, EMG.TA , TArel, color(j,:))
yline(0,'k--'); yline(50,'k--');
set(gcf,'name','TA')

end
end


%% sine conditions
sine_conditions = {'sine_020','sine_1020'};

for j = 1:2
    
    figure(10)
    color = get(gca,'colororder');
    
    % load
    load([sine_conditions{j},'_EMG.mat'])
    
    % subtract rest torque, divide by MVC torque, multiply by 100%
    MGrel = EMG.MG ./ max(MVC.MG,[],2) * 100;
    LGrel = EMG.LG ./ max(MVC.LG,[],2) * 100;
    TArel = EMG.TA ./ max(MVC.TA,[],2) * 100;
    SOrel = EMG.SO ./ max(MVC.SO,[],2) * 100;
    
   
    figure(10); plotEMG(tnew, EMG.MG , MGrel, color(j,:))
    yline(0,'k--'); yline(10,'k--'); yline(20,'k--')
    set(gcf,'name','MG')
    
    figure(11); plotEMG(tnew, EMG.LG , LGrel, color(j,:))
    yline(0,'k--'); yline(10,'k--'); yline(20,'k--')
     set(gcf,'name','LG')
     
    figure(12); plotEMG(tnew, EMG.SO , SOrel, color(j,:))
    yline(0,'k--'); yline(10,'k--'); yline(20,'k--')
     set(gcf,'name','SO')
     
    figure(13); plotEMG(tnew, EMG.TA , TArel, color(j,:))
    yline(0,'k--'); yline(10,'k--'); yline(20,'k--')
    set(gcf,'name','TA')
 
end

%% save figures
% for i = 1:3
%     figure(i); 
%     saveas(gcf, ['fig',num2str(i)], 'fig')
% end

function[] = plotEMG(t, EMG, EMGrel, color)
    subplot(211);
    plot(t, mean(EMG,'omitnan'),'color',color,'linewidth',2); hold on
    plot(t, mean(EMG,'omitnan')-std(EMG,1,'omitnan'),'-','color',color, 'linewidth',.5)
    plot(t, mean(EMG,'omitnan')+std(EMG,1,'omitnan'),'-','color',color,'linewidth',.5)
    xlabel('Time (s)'); ylabel('EMG (V)'); title('Absolute EMG'); box off; xlim([0 82])

    subplot(212);
    plot(t, mean(EMGrel,'omitnan'),'color',color,'linewidth',2); hold on
    plot(t, mean(EMGrel,'omitnan')-std(EMGrel,1,'omitnan'),'-','color',color,'linewidth',.5)
    plot(t, mean(EMGrel,'omitnan')+std(EMGrel,1,'omitnan'),'-','color',color,'linewidth',.5)
    xlabel('Time (s)'); ylabel('EMG (%MVC)'); title('Relative EMG'); box off; xlim([0 82])
 
end