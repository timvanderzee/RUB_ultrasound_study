% typical example participant, torques, angles, EMGs
% for all trials except passive

Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('MVC_EMG.mat', 'MVC', 'tnew');

%% contruct target
tpoints = [0 1 3.5 4.5 7 8; 
           0 2 3.5 4.5 6 8;
           0 3 3.5 4.5 5 8;
           0 1 3.5 4.5 4.5000001 8];
       
Target = [0 0 .5 .5 0 0];

%% time series
force_conditions = {'slow','medium','fast','asym'};
image_qualities = {'low', 'high'};

p = 1;
i = 2;

titles = {'Slow ramp', 'Medium ramp', 'Fast ramp', 'Asymmetric ramp'};

figure(2)
color = get(gca,'colororder');
N = 7;
M = length(force_conditions);

for j = 1:M
    
% load
load([force_conditions{j},'_',image_qualities{i},'_summary.mat'], 'torque')

% load
load([force_conditions{j},'_',image_qualities{i},'_EMG.mat'], 'EMG')

% subtract rest torque, divide by MVC torque, multiply by 100%
MGrel = EMG.MG.raw ./ max(MVC.MG,[],2) * 100;
LGrel = EMG.LG.raw ./ max(MVC.LG,[],2) * 100;
TArel = EMG.TA.raw ./ max(MVC.TA,[],2) * 100;
SOrel = EMG.SO.raw ./ max(MVC.SO,[],2) * 100;

MGrel_filt = EMG.MG.filt ./ max(MVC.MG,[],2) * 100;
LGrel_filt = EMG.LG.filt ./ max(MVC.LG,[],2) * 100;
TArel_filt = EMG.TA.filt ./ max(MVC.TA,[],2) * 100;
SOrel_filt = EMG.SO.filt ./ max(MVC.SO,[],2) * 100;

% subtract rest torque, divide by MVC torque, multiply by 100%
Trel = (torque(p,:) - Trest(p)) ./ Tmax(p) * 100;

t = tnew - 1;

% subplot(length(force_conditions)*2,1,j*2-1);
subplot(N, length(force_conditions), j)
plot(t, interp1(tpoints(j,:), Target * Tmax(p), t), '--','color', [.5 .5 .5],'linewidth',2); hold on
plot(t, torque(p,:)-Trest(p),'linewidth',2,'color',color(5,:)); hold on

box off; 
ylim([-5 120])

title(titles{j})

% subplot(length(force_conditions)*2,1,j*2);
subplot(N, M, j+M)
plot(t, MGrel(p,:),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, MGrel_filt(p,:),'linewidth',2,'color',color(5,:)); hold on
box off; 
ylim([-200 200])

subplot(N, M, j+M*2)
plot(t, LGrel(p,:),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, LGrel_filt(p,:),'linewidth',2,'color',color(5,:)); hold on
box off; 
ylim([-200 200])

subplot(N, M, j+M*3)
plot(t, SOrel(p,:),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, SOrel_filt(p,:),'linewidth',2,'color', color(5,:)); hold on
box off; 
ylim([-200 200])

subplot(N, M, j+M*4)
plot(t, TArel(p,:),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, TArel_filt(p,:),'linewidth',2,'color',color(5,:)); hold on
box off; 
ylim([-200 200])

end

%%  ultrasound
filenames = {'*slow_high*','*medium_high*','*fast_high*','*asym_high*'}; 
dcolor = [color(2,:)+[0 .2 .2]; color(6,:); color(4,:)];

for k = 1:length(filenames)
    cd([datafolder, '\fascicle_tracking_estimates\p1'])
    files = dir(filenames{k});

    for i = 1:length(files)
        filename = files(i).name;

        if exist(filename,'file')
            load(filename);

            % recreate t
            t = 0:.03:((2667-1)*.03);

            subplot(N, M, k+M*5)
            plot(t,Fdat.Region.PEN,'color',dcolor(i,:),'linewidth',2); hold on
            ylim([15 40])
            box off

            subplot(N, M, k+M*6)  
            plot(t,Fdat.Region.FL,'color',dcolor(i,:),'linewidth',2); hold on
            ylim([35 75])
            box off
            xlabel('Time (s)')
        end
    end
end

sgtitle('Ramp-and-hold trials')

for i = 1:N*M
    subplot(N, M, i);
    xlim([0 8])
end

%%
% set(gcf,'units','normalized','position',[.2 0 .2 .9])
figure(2)
set(gcf,'units','normalized','position',[0 0 .4 .99])

subplot(N,M,1); ylabel('Torque (N-m)')
subplot(N,M,M+1); ylabel('MG (%MVC)')
subplot(N,M,M*2+1); ylabel('LG (%MVC)')
subplot(N,M,M*3+1); ylabel('SOL (%MVC)')
subplot(N,M,M*4+1); ylabel('TA (%MVC)')
subplot(N,M,M*5+1); ylabel('Pennation (deg)')
subplot(N,M,M*6+1); ylabel('Length (mm)')
