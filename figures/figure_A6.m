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

n = 2667;
tus = 0:.03:((n-1)*.03);

angle_rs = nan(n, 8, 3);
p = 1;

for j = 1:length(force_conditions)
    
    % load
    load([force_conditions{j},'_summary.mat'])

    % down-sample
    angle_rs(:,p,j) = interp1(tnew, angle(p,:)', tus,[],'extrap');

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

p = 1;

nn = 2667;
pen = nan(nn, 8, 3, 3,2);
len = nan(nn, 8, 3, 3,2);


foldername = foldernames{p};


is = [length(Qs) 1 5];
m = 0;

for i = is
    m = m+1;
    
for k = 1:length(filenames)
    cd([mainfolder foldername]);
    files = dir(filenames{k});
    vidname = files.name(1:end-4);
    
    cd([mainfolder foldername,'\analyzed\mat']);
    

        filename = [vidname,'_analyzed_Q=',strrep(num2str(Qs(i)),'.',''),'_v2'];

    if exist([filename,'.mat'],'file')
        load([filename,'.mat']);
        
        
        nn = min(2667, length(Fdat.Region.PEN));

        % save
        pen(1:nn,p,k,m) = Fdat.Region.PEN(1:nn)*180/pi;
        len(1:nn,p,k,m) = Fdat.Region.FL(1:nn);
    else
        disp('Does not exist')
    end

    
end
end


%% as a function of angle
figure(1)

n = 2667;
dcolor = [[229 197 184]/255; .5 .5 .5; color(1,:)];
dcolors = [241 224 217; [.9 .9 .9]*255; 203 235 255]/255;
NN = 100;

clear pens lens

for m = 1:3
    for i = 1:(length(force_conditions))

        [asort, ii] = sort(angle_rs(1:n,p,i)-min(angle_rs(1:n,p,i)));

        mpen = movmean(pen(ii,p,i,m),NN,'omitnan');
        spen = movstd(pen(ii,p,i,m),NN,'omitnan');

        mlen = movmean(len(ii,p,i,m),NN,'omitnan');
        slen = movstd(len(ii,p,i,m),NN,'omitnan');

        subplot(2,3,i)

        id = isfinite(asort);

        coord_combine = [[asort(id) mpen(id)+spen(id)]; flipud([asort((id)) mpen(id)-spen(id)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:),'LineStyle','none','FaceAlpha',.5); hold on


%             plot(asort, pen(ii,p,i,m,z),'.','color',dcolors(m,:)); hold on

        plot(asort(id), mpen(id),'-','color',dcolor(m,:),'linewidth',2);      hold on

        ylim([15 40])

         xlim([0 45])
           box off

        subplot(2,3,i+3)

        coord_combine = [[asort(id) mlen(id)+slen(id)]; flipud([asort((id)) mlen(id)-slen(id)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:),'LineStyle','none','FaceAlpha',.5); hold on


%             plot(asort, len(ii,p,i,m,z),'.','color',dcolors(m,:)); hold on
         plot(asort(id), mlen(id),'-','color',dcolor(m,:),'linewidth',2); hold on  

         xlim([0 45])
        ylim([35 70])
        box off
        xlabel('Angle (deg)')

        % interpolate
        alin = linspace(0,45,101);

        % save
%             angles(:,p,i) = asort;
        pens(:,p,i,m) = interp1(asort(id), mpen(id), alin);
        lens(:,p,i,m) = interp1(asort(id), mlen(id), alin);

    end
end


subplot(231); ylabel('Pennation (deg)')
title('5 deg/s')

subplot(232); title('30 deg/s')
subplot(233); title('120 deg/s')
subplot(234); ylabel('Length (mm)')

set(gcf,'units','normalized','position',[.2 .2 .4 .4])