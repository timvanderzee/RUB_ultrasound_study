clear all; close all; clc
% typical example participant, torques, angles, EMGs

Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('RampTarget.mat','tnew','rampTarget')
load('MVC_EMG.mat');

%% time series
force_conditions = {'slow','medium','fast','asym'};
image_qualities = {'low', 'high'};
i = 2;

titles = {'Slow ramp', 'Medium ramp', 'Fast ramp', 'Asymmetric ramp'};

N = 7;
M = length(force_conditions);

for p = 1:8
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
%     Trel = (torque(p,:) - Trest(p));

    % cut the first and last second
    tcut = tnew(101:end-101)-1;
    Tcut = Trel(101:end-101);
    MGrel_cut = MGrel_filt(p,101:end-101);
    LGrel_cut = LGrel_filt(p,101:end-101);
    TArel_cut = TArel_filt(p,101:end-101);
    SOrel_cut = SOrel_filt(p,101:end-101);

    T = 1/0.125;
    Ts = (0:10)*T;
    ts = linspace(0, T, 101);
    
    for n = 1:length(Ts)-1
        id = (tcut >= Ts(n)) & (tcut < Ts(n+1));
        tmod = mod(tcut(id),T);
               
        T_pc(n,:) = interp1(tmod, Tcut(id), ts,[], 'extrap');
        MG_pc(n,:) = interp1(tmod, MGrel_cut(id), ts,[], 'extrap');
        LG_pc(n,:) = interp1(tmod, LGrel_cut(id), ts,[], 'extrap');
        TA_pc(n,:) = interp1(tmod, TArel_cut(id), ts,[], 'extrap');
        SO_pc(n,:) = interp1(tmod, SOrel_cut(id), ts,[], 'extrap');
    end
    
    mu =    [mean(T_pc); mean(MG_pc); mean(LG_pc); mean(SO_pc); mean(TA_pc)];
    sigma = [std(T_pc); std(MG_pc); std(LG_pc); std(SO_pc); std(TA_pc)];
    
    % save the cycle average
    mus(:,:,p,j,i) = mu;
end
end

%% add ultrasound
N = 7;
Qs = [nan, 0, 10.^(-4:0), 1000, inf];
color = get(gca,'colororder');

mainfolder = 'C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
% mainfolder = 'C:\Users\u0167448\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
subfolders = dir(mainfolder);

foldernames = {'3011', '0812', '1312','1612','1601','1701','1901a','1901b'};

% force_conditions = {'slow','medium','fast','asym','sine_020','sine_1020'};

% fignames = {'sine 0-20', 'sine 0-10', 'passive 5 deg/s','passive 30 deg/s','passive 120 deg/s','ramp: asymmetric','ramp: slow','ramp: medium','ramp: fast'};
participants = foldernames;

%% Figure 3: cycle average
% dcolor = [.5 .5 .5; .8 .8 .8; color(1,:)];
dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
ix = [1 length(Qs) 5];
% ix = 1;
% ix = 5;
% ix = 6;

iis = [1:4, 6:length(Qs), 5];

v = 10; % number of contraction
msdphis = nan(v, 8, 6, length(Qs),2);
msdlens = nan(v, 8, 6, length(Qs),2);

mtsus = nan(101, 8, 6, length(Qs),2);
mphis = nan(101, 8, 6, length(Qs),2);
mlens = nan(101, 8, 6, length(Qs),2);

    
for p = 1:8
for j = 1:2
    if j == 1
        filenames = {'*slow_low*.mp4','*medium_low*.mp4','*fast_low*.mp4','*asym_low*.mp4'}; 
        
        if p > 4
            fs = 25;
        else
            fs = 100/3;
        end
        
    else
        filenames = {'*slow_high*.mp4','*medium_high*.mp4','*fast_high*.mp4','*asym_high*.mp4'}; 
        fs = 100/3;
    end
    m = 0;
    foldername = foldernames{p};
    
    disp(foldername)
    
for i = iis
    if ismember(i,ix)
        m = m+1;
    end
    
for k = 1:M
    cd([mainfolder foldername]);
    files = dir(filenames{k});
    vidname = files.name(1:end-4);

    filename = [vidname,'_tracked_Q=',strrep(num2str(Qs(i)),'.','')];
    cd([mainfolder foldername,'\Tracked']);
    
%     filename = [vidname,'_analyzed_Q=',strrep(num2str(Qs(i)),'.',''),'_v2'];
%     cd([mainfolder foldername,'\analyzed\mat']);

    if exist([filename,'.mat'],'file')
        load([filename,'.mat'],'Fdat');
    else
        disp(['Not found: ', filename,'.mat'])
    end

    T = 8;
    Ts = (0:10)*T;
    dt = 1/fs;
    
    if j == 1
        tall = 0:dt:(Ts(end)-dt);
    
    else
        tall = 0:dt:Ts(end);
    end
    
    msd_phi = nan(length(Ts)-1,1);
    msd_faslen = nan(length(Ts)-1,1);
    pen_rs = nan(length(Ts)-1, 101);
    len_rs = nan(length(Ts)-1, 101);
    
    for n = 1:length(Ts)-1
%         id = (tall <= Ts(n));
        id = (tall >= Ts(n)) & (tall <= Ts(n+1));

        % cut
        pen = Fdat.Region.PEN(id);
        len = Fdat.Region.FL(id);
        tnew = mod(tall(id),T);
        
        
        % unique
        [tnew,IA,IC] = unique(tnew);
        upen = pen(IA);
        ulen = len(IA);
        
%         figure(100)
%         plot(tnew, pen); hold on

        % sort
        [tnew2, is] = sort(tnew);
        
        

        % resample
        ts = linspace(0, T, 101);
        pen_rs(n,:) = interp1(tnew2, upen(is), ts,[],'extrap');
        len_rs(n,:) = interp1(tnew2, ulen(is), ts,[],'extrap');    
    
        % mean and standard deviation
        phi = mean(pen_rs,'omitnan');
        sd_phi = std(pen_rs,1,'omitnan');
        faslen = mean(len_rs,'omitnan');
        sd_faslen = std(len_rs,1,'omitnan');

        msd_phi(n,1) = mean(sd_phi,'omitnan');
        msd_faslen(n,1) = mean(sd_faslen,'omitnan');
    end
    
    % downsample for display
%     ti = linspace(min(ts), max(ts), N);
%     phi_ds = interp1(ts, phi, ti);
%     faslen_ds = interp1(ts, faslen, ti);

    x2 = 0:100;
    
    % save
    mtsus(:,p,k,i,j) = ts;
    mphis(:,p,k,i,j) = phi;
    mlens(:,p,k,i,j) = faslen;
    
    msdphis(1:length(Ts)-1,p,k,i,j) = msd_phi;
    msdlens(1:length(Ts)-1,p,k,i,j) = msd_faslen;

end
end
end
end

%% save
cd('C:\Users\timvd\Documents\RUB_ultrasound_study\figures\data')
% cd('C:\Users\u0167448\Documents\GitHub\RUB_ultrasound_study\figures\data')
save('cycle_averages_ramps.mat','msdphis','mtsus','mphis','mlens','msdlens')
