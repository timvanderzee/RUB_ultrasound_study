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

for j = 1:length(force_conditions)
    
    % load
    load([force_conditions{j},'_summary.mat'])

    % down-sample
    for p = 1:8
        angle_rs(:,p,j) = interp1(tnew, angle(p,:)', tus,[],'extrap');
    
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

p = 1;
i = 1;


nn = 2667;
pen = nan(nn, 8, 3, 3,2);
len = nan(nn, 8, 3, 3,2);

ls = {'-','--'};
z = 2;

for p = 1:8
foldername = foldernames{p};

dcolor = [.8 .8 .8; .5 .5 .5; color(1,:)];

is = [length(Qs) 1 5];
m = 0;

for i = is
    m = m+1;
    
for k = 1:length(filenames)
    cd([mainfolder foldername]);
    files = dir(filenames{k});
    vidname = files.name(1:end-4);
    
    cd([mainfolder foldername,'\analyzed\mat']);
    
    if z == 1
        filename = [vidname,'_analyzed_Q=',strrep(num2str(Qs(i)),'.','')];
    
    else
        filename = [vidname,'_analyzed_Q=',strrep(num2str(Qs(i)),'.',''),'_v2'];
%         cd([mainfolder foldername,'\try_bidirectional\tracked']);
    end
    
    if exist([filename,'.mat'],'file')
        load([filename,'.mat']);
        
        figure(p)
        nn = length(Fdat.Region.PEN);
        t = 0:.03:((nn-1)*.03);
        subplot(N, length(force_conditions), k+length(force_conditions))
        plot(t,Fdat.Region.PEN*180/pi,'color',dcolor(m,:),'linewidth',2,'linestyle',ls{z}); hold on
        ylim([10 40])
        box off

        subplot(N, length(force_conditions), k+length(force_conditions)*2)  
        plot(t,Fdat.Region.FL,'color',dcolor(m,:),'linewidth',2,'linestyle',ls{z}); hold on
        ylim([50 100])
        box off
        xlabel('Time (s)')
        
        nn = min(2667, length(Fdat.Region.PEN));

        % save
        pen(1:nn,p,k,m,z) = Fdat.Region.PEN(1:nn)*180/pi;
        len(1:nn,p,k,m,z) = Fdat.Region.FL(1:nn);
    else
        disp('Does not exist')
    end

    
%     x = angle_rs(3:2500,p,k);
%     y = pen(3:2500,p,k,m,z);
%     
%     xnorm = (x - mean(x,'omitnan'))/ std(x,1,'omitnan');
%     ynorm = (y - mean(y,'omitnan'))/ std(y,1,'omitnan');
% 
%     figure(10)
%     subplot(2, length(force_conditions), k + (z-1)*3)
%     plot(xnorm,'k','linewidth',2); hold on
%     plot(ynorm,'color',dcolor(m,:),'linewidth',2);

%     d = finddelay(xnorm,ynorm);
    
end
end
end


return
%%

for p = 1:8
    figure(p)
    
    for i = 1:(length(force_conditions)*3)
        subplot(N, length(force_conditions), i);
        xlim([2 12])
    end
end

%% determine delay
p=1;
i=1;

x = angle_rs(:,p,i);
y = pen(:,p,i,3);

figure(10)
plot((x - mean(x,'omitnan'))/ std(x,1,'omitnan')); hold on
plot((y - mean(y,'omitnan'))/ std(y,1,'omitnan'));

%% as a function of angle
close all
figure(12)

n = 2667;
dcolor = [.7 .7 .7; .5 .5 .5; color(1,:)];
dcolors = dcolor+.2;
z = 2;
close all

NN = 100;

clear pens lens

for p = 1:8
figure(p)
% figure(100)

    for m = 1:3
        for i = 1:(length(force_conditions))

            [asort, ii] = sort(angle_rs(1:n,p,i)-min(angle_rs(1:n,p,i)));
            
            mpen = movmean(pen(ii,p,i,m,z),NN,'omitnan');
            spen = movstd(pen(ii,p,i,m,z),NN,'omitnan');
            
            mlen = movmean(len(ii,p,i,m,z),NN,'omitnan');
            slen = movstd(len(ii,p,i,m,z),NN,'omitnan');
            
            subplot(2,3,i)
              
            id = isfinite(asort);
            
            coord_combine = [[asort(id) mpen(id)+spen(id)]; flipud([asort((id)) mpen(id)-spen(id)])];
            h = fill(coord_combine(:,1),coord_combine(:,2),'b');
            set(h,'FaceColor',dcolors(m,:),'LineStyle','none','FaceAlpha',.5); hold on

           
%             plot(asort, pen(ii,p,i,m,z),'.','color',dcolors(m,:)); hold on

            plot(asort(id), mpen(id),'-','color',dcolor(m,:),'linewidth',2);      hold on
            
            ylim([10 40])

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
end
%%

if ishandle(100), close(100); end
figure(100)
for m = 1:3
for i = 1:3
    
    mpens = mean(pens(:,:,i,m)-pens(1,:,i,m),2);
    spens = std(pens(:,:,i,m)-pens(1,:,i,m),1,2);
    
    mlens = mean(lens(:,:,i,m)-lens(1,:,i,m),2);
    slens = std(lens(:,:,i,m)-lens(1,:,i,m),1,2);
    
    subplot(2,3,i)
%     plot(alin, pens(:,:,i)-pens(1,:,i)); hold on
    
    coord_combine = [[alin(:) mpens+spens]; flipud([alin(:) mpens-spens])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',dcolors(m,:),'LineStyle','none','FaceAlpha',.5); hold on

    
    subplot(2,3,i+3)
%     plot(alin, lens(:,:,i)-lens(1,:,i)); hold on
    coord_combine = [[alin(:) mlens+slens]; flipud([alin(:) mlens-slens])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',dcolors(m,:),'LineStyle','none','FaceAlpha',.5); hold on


end
end

for m = 1:3
    for i = 1:3
            mpens = mean(pens(:,:,i,m)-pens(1,:,i,m),2);
    spens = std(pens(:,:,i,m)-pens(1,:,i,m),1,2);
    
    mlens = mean(lens(:,:,i,m)-lens(1,:,i,m),2);
    slens = std(lens(:,:,i,m)-lens(1,:,i,m),1,2);
         subplot(2,3,i)
    plot(alin, mpens,'linewidth',2,'color',dcolor(m,:))
                ylim([0 20])

             xlim([0 45])
               box off
    
     subplot(2,3,i+3)
     plot(alin,mlens,'linewidth',2,'color',dcolor(m,:))
                  xlim([0 45])
            ylim([-30 0])
            box off
            xlabel('Angle (deg)')
            
            
    end
end

subplot(231); ylabel('Pennation (deg)')
title('5 deg/s')

subplot(232); title('30 deg/s')
subplot(233); title('120 deg/s')
subplot(234); ylabel('Length (mm)')
%%
i=1;
p=1;
m=3;

shifts = -20:30;

for j = 1:length(shifts)
shift = shifts(j);

figure(100)
subplot(311);
plot(angle_rs(100:(2500-shift),p,i), pen((100+shift):2500,p,i,m,1),'.','color',color(m,:)); 
title(num2str(shift))

subplot(312);
plot(angle_rs(100:(2500-shift),p,i), pen((100+shift):2500,p,i,m,2),'.','color',color(m,:))

subplot(313);
plot(angle_rs(100:(2500-shift),p,i), pen((100+shift):2500,p,i,2,2),'.','color',color(m,:))

pause
end
