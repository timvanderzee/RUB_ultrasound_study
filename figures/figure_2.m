clear all; close all; clc
% typical example participant, torques, angles, EMGs
% for all trials except passive

Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('RampTarget.mat','tnew','rampTarget')
load('MVC_EMG.mat');

%% time series
force_conditions = {'sine_020','sine_1020'};
image_qualities = {'low', 'high'};

p = 1;
i = 2;

titles = {'0-20 %MVC sine','10-20 %MVC sine'};

close all

figure(p)
color = get(gca,'colororder');
N = 7;
M = length(force_conditions);

for j = 1:M
    
% load
load([force_conditions{j},'_',image_qualities{i},'_summary.mat'])

% load
load([force_conditions{j},'_',image_qualities{i},'_EMG.mat'])

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

% subplot(M*2,1,j*2-1);
subplot(N, M, j)
plot(t, torque(p,:)-Trest(p),'linewidth',2); hold on
box off; 
ylim([-5 120])

title(titles{j})

% subplot(M*2,1,j*2);
subplot(N, M, j+M)
plot(t, MGrel(p,:),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, MGrel_filt(p,:),'linewidth',2,'color',color(1,:)); hold on
box off; 
ylim([-200 200])

subplot(N, M, j+M*2)
plot(t, LGrel(p,:),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, LGrel_filt(p,:),'linewidth',2,'color',color(1,:)); hold on
box off; 
ylim([-200 200])

subplot(N, M, j+M*3)
plot(t, SOrel(p,:),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, SOrel_filt(p,:),'linewidth',2,'color', color(1,:)); hold on
box off; 
ylim([-200 200])

subplot(N, M, j+M*4)
plot(t, TArel(p,:),'linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, TArel_filt(p,:),'linewidth',2,'color',color(1,:)); hold on
box off; 
ylim([-200 200])

end

%% add ultrasound
N = 7;
Qs = [nan, 0, 10.^(-4:0), 1000, inf];
color = get(gca,'colororder');

mainfolder = 'C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
% mainfolder = 'C:\Users\u0167448\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
subfolders = dir(mainfolder);

foldernames = {'3011', '0812', '1312','1612','1601','1701','1901a','1901b'};
filenames = {'*sine_020*.mp4','*sine_1020*.mp4'}; 

participants = foldernames;

j = 1;
i = 1;

foldername = foldernames{j};

dcolor = [5 .5 .5; .8 .8 .8; color(1,:)];

is = [1 length(Qs) 5];
m = 0;
for i = is
    if ismember(i,is)
        m = m+1;
    end
for k = 1:length(filenames)
    cd([mainfolder foldername]);
    files = dir(filenames{k});
    vidname = files.name(1:end-4);

    filename = [vidname,'_analyzed_Q=',strrep(num2str(Qs(i)),'.',''),'_v2'];

    cd([mainfolder foldername,'\analyzed\mat']);

    if exist([filename,'.mat'],'file')
        load([filename,'.mat']);
        
        % recreate t
        t = 0:.03:((2667-1)*.03);

        subplot(N, M, k+M*5)
        plot(t,Fdat.Region.PEN*180/pi,'color',dcolor(m,:),'linewidth',2); hold on
        ylim([15 40])
        box off

        subplot(N, M, k+M*6)  
        plot(t,Fdat.Region.FL,'color',dcolor(m,:),'linewidth',2); hold on
        ylim([35 80])
        box off
        xlabel('Time (s)')
    end


end
end

%%
for i = 1:N*M
    subplot(N, M, i);
    xlim([0 8])
end


% set(gcf,'units','normalized','position',[.2 0 .2 .9])
figure(1)
set(gcf,'units','normalized','position',[0 0 .4 .99])

subplot(N,M,1); ylabel('Torque (N-m)')
subplot(N,M,M+1); ylabel('MG (%MVC)')
subplot(N,M,M*2+1); ylabel('LG (%MVC)')
subplot(N,M,M*3+1); ylabel('SOL (%MVC)')
subplot(N,M,M*4+1); ylabel('TA (%MVC)')
subplot(N,M,M*5+1); ylabel('Pennation (deg)')
subplot(N,M,M*6+1); ylabel('Length (mm)')

