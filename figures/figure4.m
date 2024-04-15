clear all; close all; clc
% cycle- and participant average
% passive trials

Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('RampTarget.mat','tnew','rampTarget')
load('MVC_EMG.mat');

%% time series
force_conditions = {'pas_005', 'pas_30','pas_120'};

titles = {'Slow passive', 'Medium passive', 'Fast passive'};

close all
N = 3;

fs = 100;
Wn = 10 / (.5*fs);
[b,a] = butter(2, Wn);

n = 8201;
tus = 0:.03:((n-1)*.03);

angle_rs = nan(n, 8, 3);

for j = 1:length(force_conditions)
    
    % load
    load([force_conditions{j},'_summary.mat'])

    % down-sample
    for p = 1:8
        angle_rs(:,p,j) = interp1(tnew, angle(p,:)', tus);
    
        figure(p);
        subplot(N,3,j);
        plot(tnew, angle(p,:)); hold on
        plot(tus, angle_rs(:,p,j), '.','markersize',5)
        box off
    end
end

%% ultrasound
Qs = [nan, 0, 10.^(-4:0), 1000, inf];
color = get(gca,'colororder');

mainfolder = 'C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
% mainfolder = 'C:\Users\u0167448\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
subfolders = dir(mainfolder);

foldernames = {'3011', '0812', '1312','1612','1601','1701','1901a','1901b'};
filenames = {'*pas_005*.mp4','*pas_30*.mp4','*pas_120*.mp4'};

participants = foldernames;

j = 1;
i = 1;



ls = {'-','--'};
for z = 1:2
n = 8201;
pen = nan(n, 8, 3, 3);
len = nan(n, 8, 3, 3);

for j = 1
foldername = foldernames{j};

dcolor = [.8 .8 .8; .5 .5 .5; color(1,:)];

is = [length(Qs) 1 5];
m = 0;

for i = is
    m = m+1;
    
for k = 1:length(filenames)
    cd([mainfolder foldername]);
    files = dir(filenames{k});
    vidname = files.name(1:end-4);

    if z == 1
        filename = [vidname,'_analyzed_Q=',strrep(num2str(Qs(i)),'.','')];
        cd([mainfolder foldername,'\analyzed\mat']);
    else
        filename = [vidname,'_tracked_Q=',strrep(num2str(Qs(i)),'.','')];
        cd([mainfolder foldername,'\try_bidirectional\tracked']);
    end
    
    if exist([filename,'.mat'],'file')
        load([filename,'.mat']);
        
        figure(j)
        n = length(Fdat.Region.PEN);
        t = 0:.03:((n-1)*.03);
        subplot(N, length(force_conditions), k+length(force_conditions))
        plot(t,Fdat.Region.PEN*180/pi,'color',dcolor(m,:),'linewidth',2,'linestyle',ls{z}); hold on
        ylim([10 40])
        box off

        subplot(N, length(force_conditions), k+length(force_conditions)*2)  
        plot(t,Fdat.Region.FL,'color',dcolor(m,:),'linewidth',2,'linestyle',ls{z}); hold on
        ylim([50 100])
        box off
        xlabel('Time (s)')
        
        % save
        pen(1:n,j,k,m) = Fdat.Region.PEN(:)*180/pi;
        len(1:n,j,k,m) = Fdat.Region.FL(:);
    end


end
end
end
end
%%
for p = 1:8
    figure(p)
    
    for i = 1:(length(force_conditions)*3)
        subplot(N, length(force_conditions), i);
        xlim([2 12])
    end
end

%% as a function of angle
close all
figure(12)
for p = 1
% figure(p+10)

    for m = 1:3
        for i = 1:(length(force_conditions))

        subplot(2,3,i)
        plot(angle_rs(:,p,i), pen(:,p,i,m),'.','color',color(m,:)); hold on

        subplot(2,3,i+3)
        plot(angle_rs(:,p,i), len(:,p,i,m),'.','color',color(m,:)); hold on
        end
    end
end
