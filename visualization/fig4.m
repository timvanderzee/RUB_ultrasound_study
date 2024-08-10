close all; clc
% typical example participant, torques, angles, EMGs
% for passive trials

Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('MVC_EMG.mat', 'MVC', 'tnew');

%% time series
force_conditions = {'pas_005', 'pas_30','pas_120'};
p = 1;

titles = {'Slow passive', 'Medium passive', 'Fast passive'};

figure(1)
color = get(gca,'colororder');
N = 7;
M = length(force_conditions);

for j = 1:length(force_conditions)
    
% load
load([force_conditions{j},'_summary.mat'], 'angle')

% load
load([force_conditions{j},'_EMG.mat'], 'EMG')

% subtract rest torque, divide by MVC torque, multiply by 100%
MGrel = EMG.MG.raw ./ max(MVC.MG,[],2) * 100;
LGrel = EMG.LG.raw ./ max(MVC.LG,[],2) * 100;
TArel = EMG.TA.raw ./ max(MVC.TA,[],2) * 100;
SOrel = EMG.SO.raw ./ max(MVC.SO,[],2) * 100;

MGrel_filt = EMG.MG.filt ./ max(MVC.MG,[],2) * 100;
LGrel_filt = EMG.LG.filt ./ max(MVC.LG,[],2) * 100;
TArel_filt = EMG.TA.filt ./ max(MVC.TA,[],2) * 100;
SOrel_filt = EMG.SO.filt ./ max(MVC.SO,[],2) * 100;

t = tnew(1:5:end);

% subplot(length(force_conditions)*2,1,j*2-1);
subplot(N, length(force_conditions), j)
plot(t, angle(p,1:5:end),'linewidth',2,'color',color(5,:)); hold on
box off; 
% ylim([-5 120])

title(titles{j})

% subplot(length(force_conditions)*2,1,j*2);
subplot(N, length(force_conditions), j+length(force_conditions))
plot(t, MGrel(p,1:5:end),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, MGrel_filt(p,1:5:end),'linewidth',2,'color',color(5,:)); hold on
box off; 
ylim([-200 200])

subplot(N, length(force_conditions), j+length(force_conditions)*2)
plot(t, LGrel(p,1:5:end),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, LGrel_filt(p,1:5:end),'linewidth',2,'color',color(5,:)); hold on
box off; 
ylim([-200 200])

subplot(N, length(force_conditions), j+length(force_conditions)*3)
plot(t, SOrel(p,1:5:end),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, SOrel_filt(p,1:5:end),'linewidth',2,'color', color(5,:)); hold on
box off; 
ylim([-200 200])

subplot(N, length(force_conditions), j+length(force_conditions)*4)
plot(t, TArel(p,1:5:end),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, TArel_filt(p,1:5:end),'linewidth',2,'color',color(5,:)); hold on
box off; 
ylim([-200 200])

end

%%
for i = 1:(length(force_conditions)*N)
    subplot(N, length(force_conditions), i);
    xlim([0 12])
    
    if i < 4
        ylim([-10 50])
    end
    
    if i > 3 && i < 16
        ylim([-10 10])
    end
end

%% ultrasound
filenames = {'*pas_005*','*pas_30*','*pas_120*'};

dcolor = [color(2,:)+[0 .2 .2]; color(6,:); color(4,:)];

for k = 1:length(filenames)
    cd([mainfolder, 'RUB_ultrasound_study\data\ultrasound\GM\p1'])
    files = dir(filenames{k});

    for i = 1:length(files)
        filename = files(i).name;

        if exist(filename,'file')
            load(filename);

            % recreate t
            t = 0:.03:((length(Fdat.Region.PEN)-1)*.03);

            subplot(N, M, k+M*5)
            plot(t,Fdat.Region.PEN,'color',dcolor(i,:),'linewidth',2); hold on
            ylim([20 40])
            box off

            subplot(N, M, k+M*6)  
            plot(t,Fdat.Region.FL,'color',dcolor(i,:),'linewidth',2); hold on
            ylim([35 75])
            box off
            xlabel('Time (s)')
        end
    end
end

sgtitle('Passive rotation trials')

%% make nice
for i = 1:(length(force_conditions)*N)
    subplot(N, length(force_conditions), i);
    xlim([0 16])
end

figure(1)
set(gcf,'units','normalized','position',[0 0 .4 .99])

subplot(N,3,1); ylabel('Angle (deg)')
subplot(N,3,4); ylabel('MG (%MVC)')
subplot(N,3,7); ylabel('LG (%MVC)')
subplot(N,3,10); ylabel('SOL (%MVC)')
subplot(N,3,13); ylabel('TA (%MVC)')
subplot(N,3,16); ylabel('Pennation (deg)')
subplot(N,3,19); ylabel('Length (mm)')


