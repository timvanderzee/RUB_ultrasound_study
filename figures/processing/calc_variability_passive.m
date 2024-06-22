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

nn = 2667;
pen = nan(nn, 8, 3, 3);
len = nan(nn, 8, 3, 3);

for p = 1:8
foldername = foldernames{p};

is = [1 length(Qs) 5];
m = 0;

for i = is
    m = m+1;
    
for k = 1:length(filenames)
    cd([mainfolder foldername]);
    files = dir(filenames{k});
    vidname = files.name(1:end-4);
    
    filename = [vidname,'_tracked_Q=',strrep(num2str(Qs(i)),'.','')];
    cd([mainfolder foldername,'\Tracked']);
    
%     cd([mainfolder foldername,'\analyzed\mat']);
%         filename = [vidname,'_analyzed_Q=',strrep(num2str(Qs(i)),'.',''),'_v2'];
    
    if exist([filename,'.mat'],'file')
        load([filename,'.mat']);
        
        nn = min(2667, length(Fdat.Region.PEN));

        % save
        pen(1:nn,p,k,m) = Fdat.Region.PEN(1:nn);
        len(1:nn,p,k,m) = Fdat.Region.FL(1:nn);
    else
        disp('Does not exist')
    end

  
end
end
end


%%


% 
% if ishandle(1), close(1); end
% figure(1)

titles = {'Slow','Medium','Fast'};


dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];

mpen = nan(50,8,3,3);
spen = nan(50,8,3,3);
mlen = nan(50,8,3,3);
slen = nan(50,8,3,3);

for p = 1:8
% for a = 1:3
%     
%     for k = 1:3
% 
%         subplot(1,3,k)
%         plot(angle_rs(:,p,k), pen(:,p,k,a),'.'); hold on
%         title(titles{k})
%     end
%     
% end


for a = 1:3
        for k = 1:3

            % fit polynomial


            
%             id = isfinite(angle_rs(:,p,k)) & isfinite(pen(:,p,k,a));

                as = floor(min(angle_rs(:,p,k))):1:(ceil(max(angle_rs(:,p,k)))-1);
            for i = 1:length(as)-1
                id = (angle_rs(:,p,k) > as(i)) & (angle_rs(:,p,k) < as(i+1));
            
                mpen(i,p,k,a) = mean(pen(id,p,k,a));
                spen(i,p,k,a) = std(pen(id,p,k,a));
                
                mlen(i,p,k,a) = mean(len(id,p,k,a));
                slen(i,p,k,a) = std(len(id,p,k,a));
            end

                figure(p)
                subplot(1,3,k)
              errorbar(as(1:end-1), mpen(1:(length(as)-1),p,k,a), spen(1:(length(as)-1),p,k,a), 'color', dcolor(a,:)); hold on
                
%             coef = polyfit(angle_rs(id,p,k), pen(id,p,k,a), 3);
% 
%             subplot(1,3,k)
%             angle_lin = linspace(min(angle_rs(id,p,k)), max(angle_rs(id,p,k)), 100);
%             plot(angle_lin, polyval(coef, angle_lin),'k-')
%             
%             pen_pred = polyval(coef, angle_rs(id,p,k));
            
%             pen_sd(a,p,k) = mean(sqrt((pen_pred - pen(id,p,k,a)).^2));
        
        end
end
end
% legend('TT','UT','UTT')

%% average over angles
for p = 1:8
    for a = 1:3
        for k = 1:3
            mspen(a,p,k) = mean(spen(:,p,k,a),'omitnan');
            mslen(a,p,k) = mean(slen(:,p,k,a),'omitnan');
        end
    end
end

%% plot
figure(100)
for a = 1:3
    errorbar(1:3, squeeze(mean(mspen(a,:,:),2)), squeeze(std(mspen(a,:,:),1,2))); hold on
end

%%
cd('C:\Users\timvd\Documents\RUB_ultrasound_study\figures\data')
% cd('C:\Users\u0167448\Documents\GitHub\RUB_ultrasound_study\figures\data')
save('variability_passive.mat')
%%
% figure(2)
% for k = 1:3
%     subplot(1,3,k)
%     errorbar(1:3, mean(pen_sd(:,:,k),2), std(pen_sd(:,:,k),1,2));
%     ylim([0 3])
% end
