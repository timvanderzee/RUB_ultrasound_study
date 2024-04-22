clear all; close all; clc
% typical example participant, torques, angles, EMGs

Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('RampTarget.mat','tnew','rampTarget')
load('MVC_EMG.mat');

figure(1)
color = get(gca,'colororder');

dcolor = [.5 .5 .5; [229 197 184]/255; color(1,:)];
dcolors = [[.9 .9 .9]*255; 241 224 217; 203 235 255]/255;

%% time series
force_conditions = {'slow','medium','fast','asym'};
image_qualities = {'low', 'high'};
i = 2;

titles = {'Slow ramp', 'Medium ramp', 'Fast ramp', 'Asymmetric ramp'};

N = 7;
M = length(force_conditions);

p = 1;

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
        
    figure(p)
    color = get(gca,'colororder');
    
    for kk = 1:size(mu,1)
    
        x1 = 0:100;

        subplot(N, M, j+M*(kk-1))
        coord_combine = [[x1; mu(kk,:)+sigma(kk,:)], fliplr([x1; mu(kk,:)-sigma(kk,:)])];
        h = fill(coord_combine(1,:),coord_combine(2,:),'b');
        set(h,'FaceColor',dcolors(end,:),'LineStyle','none'); hold on
        plot(x1, mu(kk,:),'linewidth',2,'color',dcolor(end,:)); hold on
        ylim([-2 60])
    end
    
    % save the cycle average
    mus(:,:,p,j) = mu;
end


%% add ultrasound
N = 7;
Qs = [nan, 0, 10.^(-4:0), 1000, inf];
color = get(gca,'colororder');

mainfolder = 'C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
% mainfolder = 'C:\Users\u0167448\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
subfolders = dir(mainfolder);

foldernames = {'3011', '0812', '1312','1612','1601','1701','1901a','1901b'};
filenames = {'*slow*.mp4','*medium*.mp4','*fast*.mp4','*asym*.mp4'}; 

% force_conditions = {'slow','medium','fast','asym','sine_020','sine_1020'};

% fignames = {'sine 0-20', 'sine 0-10', 'passive 5 deg/s','passive 30 deg/s','passive 120 deg/s','ramp: asymmetric','ramp: slow','ramp: medium','ramp: fast'};
participants = foldernames;

%% Figure 3: cycle average

ix = [length(Qs) 1 5];
% ix = 1;
% ix = 5;
% ix = 6;

iis = [1:4, 6:length(Qs), 5];

v = 10; % number of contraction
msdphis = nan(v, 8, 6, length(Qs));
msdlens = nan(v, 8, 6, length(Qs));

mtsus = nan(101, 8, 6, length(Qs));
mphis = nan(101, 8, 6, length(Qs));
mlens = nan(101, 8, 6, length(Qs));


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

    filename = [vidname,'_analyzed_Q=',strrep(num2str(Qs(i)),'.',''),'_v2'];

    cd([mainfolder foldername,'\analyzed\mat']);

    if exist([filename,'.mat'],'file')
        load([filename,'.mat'],'Fdat');
    else
        disp(['Not found: ', filename,'.mat'])
    end

    T = 8;
    Ts = (0:10)*T;
    tall = 0:.03:Ts(end);
    
    msd_phi = nan(length(Ts)-1,1);
    msd_faslen = nan(length(Ts)-1,1);
    pen_rs = nan(length(Ts)-1, 101);
    len_rs = nan(length(Ts)-1, 101);
    
    for n = 1:length(Ts)-1
%         id = (tall <= Ts(n));
        id = (tall >= Ts(n)) & (tall <= Ts(n+1));

        % cut
        pen = Fdat.Region.PEN(id)*180/pi;
        len = Fdat.Region.FL(id);
        tnew = mod(tall(id),T);
        
%         figure(100)
%         plot(tnew, pen); hold on

        % sort
        [tnew2, is] = sort(tnew);

        % resample
        ts = linspace(0, T, 101);
        pen_rs(n,:) = interp1(tnew2, pen(is), ts,[],'extrap');
        len_rs(n,:) = interp1(tnew2, len(is), ts,[],'extrap');    
    
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

    if ismember(i,ix)
        figure(p)
        subplot(N, M, k+M*5)
    %     plot(ts, Fdat.Region.PEN(is)*180/pi,'.','color',[.5 .5 .5 .5]); hold on
        coord_combine = [[x2;phi+sd_phi] fliplr([x2;phi-sd_phi])];

        h = fill(coord_combine(1,:),coord_combine(2,:),'b');
        set(h,'FaceColor',dcolors(m,:),'LineStyle','none'); hold on
        plot(x2, phi, 'linewidth',2,'color',dcolor(m,:)); hold on
        ylim([15 40])

        subplot(N, M, k+M*6)  
    %     plot(ts, Fdat.Region.FL(is),'.','color',[.5 .5 .5 .5]); hold on
        coord_combine = [[x2;faslen+sd_faslen] fliplr([x2;faslen-sd_faslen])];
              h = fill(coord_combine(1,:),coord_combine(2,:),'b');
        set(h,'FaceColor',dcolors(m,:),'LineStyle','none'); hold on
        plot(x2, faslen, '-','linewidth',2,'color',dcolor(m,:)); hold on
        
        ylim([35 70])
        xlabel('Contraction cycle (%)');
    end
    
    % save
    mtsus(:,p,k,i) = ts;
    mphis(:,p,k,i) = phi;
    mlens(:,p,k,i) = faslen;
    
    msdphis(1:length(Ts)-1,p,k,i) = msd_phi;
    msdlens(1:length(Ts)-1,p,k,i) = msd_faslen;

end
end


%%
for i = 1:N*M
    subplot(N, M, i)
    box off; 
    xlim([0 100])
end

% set(gcf,'units','normalized','position',[.2 0 .2 .9])
figure(1)
set(gcf,'units','normalized','position',[0 0 .4 .99])

subplot(N,M,1); ylabel('Torque (%MVC)')
subplot(N,M,M+1); ylabel('MG (%MVC)')
subplot(N,M,M*2+1); ylabel('LG (%MVC)')
subplot(N,M,M*3+1); ylabel('SOL (%MVC)')
subplot(N,M,M*4+1); ylabel('TA (%MVC)')
subplot(N,M,M*5+1); ylabel('Pennation (deg)')
subplot(N,M,M*6+1); ylabel('Length (mm)')
