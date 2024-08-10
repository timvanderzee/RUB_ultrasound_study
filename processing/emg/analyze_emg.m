clear all; close all; clc
load('MVC_EMG.mat');

%% ramp conditions
force_conditions = {'slow','medium','fast','asym'};
image_qualities = {'low', 'high'};

mErel = nan(4,4,2);
sErel = nan(4,4,2);

for i = 1:2
for j = 1:4
    
% load
load([force_conditions{j},'_',image_qualities{i},'_EMG.mat'])

% subtract rest torque, divide by MVC torque, multiply by 100%
MGrel = EMG.MG.filt ./ max(MVC.MG,[],2) * 100;
LGrel = EMG.LG.filt ./ max(MVC.LG,[],2) * 100;
TArel = EMG.TA.filt ./ max(MVC.TA,[],2) * 100;
SOrel = EMG.SO.filt ./ max(MVC.SO,[],2) * 100;

mErel(j,:,i) = [mean(mean(MGrel,2,'omitnan')) mean(mean(LGrel,2,'omitnan')) mean(mean(SOrel,2,'omitnan')) mean(mean(TArel,2,'omitnan'))];
sErel(j,:,i) = [std(mean(MGrel,2,'omitnan'))  std(mean(LGrel,2,'omitnan'))  std(mean(SOrel,2,'omitnan'))  std(mean(TArel,2,'omitnan'))];

figure(i)
color = get(gca,'colororder');

figure(i); plotEMG(tnew, EMG.MG.filt , MGrel, color(j,:))
yline(0,'k--'); yline(50,'k--');
set(gcf,'name','MG')

figure(i+2); plotEMG(tnew, EMG.LG.filt , LGrel, color(j,:))
yline(0,'k--'); yline(50,'k--');
 set(gcf,'name','LG')

figure(i+4); plotEMG(tnew, EMG.SO.filt , SOrel, color(j,:))
yline(0,'k--'); yline(50,'k--');
 set(gcf,'name','SO')

figure(i+6); plotEMG(tnew, EMG.TA.filt , TArel, color(j,:))
yline(0,'k--'); yline(50,'k--');
set(gcf,'name','TA')

mMGrel(j,i) = mean(MGrel(:),'omitnan');
mLGrel(j,i) = mean(LGrel(:),'omitnan');
mTArel(j,i) = mean(TArel(:),'omitnan');
mSOrel(j,i) = mean(SOrel(:),'omitnan');

end
end

mCalfrel = mean([mMGrel mLGrel mSOrel],'all')
sCalfrel = std([mMGrel mLGrel mSOrel],[],'all')

mDorsrel = mean(mTArel,'all')
sDorsrel = std(mTArel,[],'all')


%% sine conditions
sine_conditions = {'sine_020','sine_1020'};
mErel = nan(2,4);
sErel = nan(2,4);


for j = 1:2
    
    

    figure(10)
    color = get(gca,'colororder');
    
    % load
    load([sine_conditions{j},'_EMG.mat'])
    
    % subtract rest torque, divide by MVC torque, multiply by 100%
    MGrel = EMG.MG.filt ./ max(MVC.MG,[],2) * 100;
    LGrel = EMG.LG.filt ./ max(MVC.LG,[],2) * 100;
    TArel = EMG.TA.filt ./ max(MVC.TA,[],2) * 100;
    SOrel = EMG.SO.filt ./ max(MVC.SO,[],2) * 100;
    
   
    mErel(j,:) = [mean(mean(MGrel,2,'omitnan')) mean(mean(LGrel,2,'omitnan')) mean(mean(SOrel,2,'omitnan')) mean(mean(TArel,2,'omitnan'))];
    sErel(j,:) = [std(mean(MGrel,2,'omitnan'))  std(mean(LGrel,2,'omitnan'))  std(mean(SOrel,2,'omitnan'))  std(mean(TArel,2,'omitnan'))];
    
    figure(10); plotEMG(tnew, EMG.MG.filt , MGrel, color(j,:))
    yline(0,'k--'); yline(10,'k--'); yline(20,'k--')
    set(gcf,'name','MG')
    
    figure(11); plotEMG(tnew, EMG.LG.filt , LGrel, color(j,:))
    yline(0,'k--'); yline(10,'k--'); yline(20,'k--')
     set(gcf,'name','LG')
     
    figure(12); plotEMG(tnew, EMG.SO.filt , SOrel, color(j,:))
    yline(0,'k--'); yline(10,'k--'); yline(20,'k--')
     set(gcf,'name','SO')
     
    figure(13); plotEMG(tnew, EMG.TA.filt , TArel, color(j,:))
    yline(0,'k--'); yline(10,'k--'); yline(20,'k--')
    set(gcf,'name','TA')
 
end


%% passive conditions
pas_conditions = {'pas_005','pas_30','pas_120'};
mErel = nan(3,4);
sErel = nan(3,4);

for j = 1:3
    
   
    figure(10)
    color = get(gca,'colororder');
    
    % load
    load([pas_conditions{j},'_EMG.mat'])
    
    % subtract rest torque, divide by MVC torque, multiply by 100%
    MGrel = EMG.MG.filt ./ max(MVC.MG,[],2) * 100;
    LGrel = EMG.LG.filt ./ max(MVC.LG,[],2) * 100;
    TArel = EMG.TA.filt ./ max(MVC.TA,[],2) * 100;
    SOrel = EMG.SO.filt ./ max(MVC.SO,[],2) * 100;
    
       
    mErel(j,:) = [mean(mean(MGrel,2,'omitnan')) mean(mean(LGrel,2,'omitnan')) mean(mean(SOrel,2,'omitnan')) mean(mean(TArel,2,'omitnan'))];
    sErel(j,:) = [std(mean(MGrel,2,'omitnan'))  std(mean(LGrel,2,'omitnan'))  std(mean(SOrel,2,'omitnan'))  std(mean(TArel,2,'omitnan'))];
    
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