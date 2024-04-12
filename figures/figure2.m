clear all; close all; clc
% typical example participant, torques, angles, EMGs

Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('RampTarget.mat','tnew','rampTarget')
load('MVC_EMG.mat');

%% time series
force_conditions = {'slow','medium','fast','asym','sine_020','sine_1020'};
image_qualities = {'low', 'high'};
i = 2;

titles = {'Slow ramp', 'Medium ramp', 'Fast ramp', 'Asymmetric ramp', '0-20 %MVC sine','10-20 %MVC sine'};

close all
N = 7;

for p = 1:8
for j = 1:length(force_conditions)
    
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
    % Trel = (torque(p,:) - Trest(p)) ./ Tmax(p) * 100;
    Trel = (torque(p,:) - Trest(p));

    % cut the first and last second
    tcut = tnew(101:end-101)-1;
    Tcut = Trel(101:end-101);
    MGrel_cut = MGrel_filt(p,101:end-101);
    LGrel_cut = LGrel_filt(p,101:end-101);
    TArel_cut = TArel_filt(p,101:end-101);
    SOrel_cut = SOrel_filt(p,101:end-101);

    % cut into 8 sec intervals
    f = 0.125;
    t = mod(tcut,1/f);

    % reshape
    ts = reshape(t,     [length(t)/10, 10]);
    Ts = reshape(Tcut,  [length(t)/10, 10]);
    MGs = reshape(MGrel_cut,  [length(t)/10, 10]);
    LGs = reshape(LGrel_cut,  [length(t)/10, 10]);
    TAs = reshape(TArel_cut,  [length(t)/10, 10]);
    SOs = reshape(SOrel_cut,  [length(t)/10, 10]);

    figure(p)
    color = get(gca,'colororder');

    subplot(N, length(force_conditions), j)
    plot(mean(ts,2), mean(Ts,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(ts,2), mean(Ts,2)+std(Ts,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(ts,2), mean(Ts,2)-std(Ts,1,2),'linewidth',1,'color',color(1,:))
    ylim([-5 120])

    title(titles{j})

    subplot(N, length(force_conditions), j+length(force_conditions)*1)
    plot(mean(ts,2), mean(MGs,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(ts,2), mean(MGs,2)+std(MGs,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(ts,2), mean(MGs,2)-std(MGs,1,2),'linewidth',1,'color',color(1,:))
    ylim([0 60])

    subplot(N, length(force_conditions), j+length(force_conditions)*2)
    plot(mean(ts,2), mean(LGs,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(ts,2), mean(LGs,2)+std(LGs,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(ts,2), mean(LGs,2)-std(LGs,1,2),'linewidth',1,'color',color(1,:))
    ylim([0 60])

    subplot(N, length(force_conditions), j+length(force_conditions)*3)
    plot(mean(ts,2), mean(SOs,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(ts,2), mean(SOs,2)+std(SOs,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(ts,2), mean(SOs,2)-std(SOs,1,2),'linewidth',1,'color',color(1,:))
    ylim([0 60])

    subplot(N, length(force_conditions), j+length(force_conditions)*4)
    plot(mean(ts,2), mean(TAs,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(ts,2), mean(TAs,2)+std(TAs,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(ts,2), mean(TAs,2)-std(TAs,1,2),'linewidth',1,'color',color(1,:))
    ylim([0 60])
    
    % save the cycle average
    mts(:,p,j) = mean(ts,2);
    mTs(:,p,j) = mean(Ts./ Tmax(p) * 100,2);
    mMGs(:,p,j) = mean(MGs,2);
    mLGs(:,p,j) = mean(LGs,2);
    mSOs(:,p,j) = mean(SOs,2);
    mTAs(:,p,j) = mean(TAs,2);

end


%%
% set(gcf,'units','normalized','position',[.2 0 .2 .9])
figure(p)

subplot(N,6,1); ylabel('Torque (N-m)')
subplot(N,6,7); ylabel('MG (%MVC)')
subplot(N,6,13); ylabel('LG (%MVC)')
subplot(N,6,19); ylabel('SOL (%MVC)')
subplot(N,6,25); ylabel('TA (%MVC)')
end

%% add ultrasound
N = 7;
Qs = [nan, 0, 10.^(-4:0), 1000, inf];
color = get(gca,'colororder');

% mainfolder = 'C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
mainfolder = 'C:\Users\u0167448\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
subfolders = dir(mainfolder);

foldernames = {'3011', '0812', '1312','1612','1601','1701','1901a','1901b'};
filenames = {'*slow*.mp4','*medium*.mp4','*fast*.mp4','*asym*.mp4','*sine_020*.mp4','*sine_1020*.mp4'}; 

% force_conditions = {'slow','medium','fast','asym','sine_020','sine_1020'};

fignames = {'sine 0-20', 'sine 0-10', 'passive 5 deg/s','passive 30 deg/s','passive 120 deg/s','ramp: asymmetric','ramp: slow','ramp: medium','ramp: fast'};
participants = foldernames;

j = 1;
i = 1;

dcolor = [.8 .8 .8; .5 .5 .5; color(1,:)];

ix = [length(Qs) 1 5];
ix = 5;
m = 2;

for p = 1:8
    m = 2;
    foldername = foldernames{p};
for i = ix
    m = m+1;
    
for k = 1:length(force_conditions)
    cd([mainfolder foldername]);
    files = dir(filenames{k});
    vidname = files.name(1:end-4);

    filename = [vidname,'_analyzed_Q=',strrep(num2str(Qs(i)),'.','')];

    cd([mainfolder foldername,'\analyzed\mat']);

    if exist([filename,'.mat'],'file')
        load([filename,'.mat'],'Fdat');
    else
        disp(['Not found: ', filename,'.mat'])
    end

    f = 0.125;
    t = 0:.03:((2667-1)*.03);
    
    % cut
    tnew = mod(t,1/f);

    % interpolate
    [ts, is] = sort(tnew);

    % moving average
    phi = movmean(Fdat.Region.PEN(is)*180/pi, 100);
    faslen = movmean(Fdat.Region.FL(is), 100);

    % moving s.d.
    sd_phi = movstd(Fdat.Region.PEN(is)*180/pi, 100);
    sd_faslen = movstd(Fdat.Region.FL(is), 100);

    % downsample for display
%     ti = linspace(min(ts), max(ts), N);
%     phi_ds = interp1(ts, phi, ti);
%     faslen_ds = interp1(ts, faslen, ti);

    figure(p)
    subplot(N, length(force_conditions), k+length(force_conditions)*5)
%     plot(ts, Fdat.Region.PEN(is)*180/pi,'.','color',[.5 .5 .5 .5]); hold on
    plot(ts, phi, 'linewidth',2,'color',dcolor(m,:)); hold on
    plot(ts, phi + sd_phi,'-','color',dcolor(m,:));
    plot(ts, phi - sd_phi,'-','color',dcolor(m,:));
    ylim([10 40])

    subplot(N, length(force_conditions), k+length(force_conditions)*6)  
%     plot(ts, Fdat.Region.FL(is),'.','color',[.5 .5 .5 .5]); hold on
    plot(ts, faslen, '-','linewidth',2,'color',dcolor(m,:)); hold on
    plot(ts, faslen + sd_faslen,'-','color',dcolor(m,:));
    plot(ts, faslen - sd_faslen,'-','color',dcolor(m,:));
    ylim([50 100])
    xlabel('Time (s)');
    
    % save
    mtsus(:,p,k) = ts;
    mphis(:,p,k) = phi;
    mlens(:,p,k) = faslen;

end
end
end

%%
for i = 1:42
    subplot(N, length(force_conditions), i)
    box off; 
    xlim([0 8])
end

subplot(N,6,31); ylabel('Pennation (deg)')
subplot(N,6,37); ylabel('Length (mm)')
set(gcf,'units','normalized','position',[0 .1 .6 .8])


%% the average figure
if ishandle(10), close(10); end

figure(10)
for j = 1:length(force_conditions)
    
    mt = mts(:,:,j);
    mT = mTs(:,:,j);
    mMG = mMGs(:,:,j);
    mLG = mLGs(:,:,j);
    mSO = mSOs(:,:,j);
    mTA = mTAs(:,:,j);
    
    mtsu = mtsus(:,:,j);
    mphi = mphis(:,:,j) - mphis(1,:,j) +  mean(mphis(1,:,j),2);
    mlen = mlens(:,:,j) - mlens(1,:,j) +  mean(mlens(1,:,j),2);

    subplot(N, length(force_conditions), j)
    plot(mean(mt,2), mean(mT,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(mt,2), mean(mT,2)+std(mT,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(mt,2), mean(mT,2)-std(mT,1,2),'linewidth',1,'color',color(1,:))
    title(titles{j})
     ylim([0 60])
     
    subplot(N, length(force_conditions), j+length(force_conditions)*1)
    plot(mean(mt,2), mean(mMG,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(mt,2), mean(mMG,2)+std(mMG,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(mt,2), mean(mMG,2)-std(mMG,1,2),'linewidth',1,'color',color(1,:))
     ylim([0 60])
     
    subplot(N, length(force_conditions), j+length(force_conditions)*2)
    plot(mean(mt,2), mean(mLG,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(mt,2), mean(mLG,2)+std(mLG,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(mt,2), mean(mLG,2)-std(mLG,1,2),'linewidth',1,'color',color(1,:))
     ylim([0 60])
     
    subplot(N, length(force_conditions), j+length(force_conditions)*3)
    plot(mean(mt,2), mean(mSO,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(mt,2), mean(mSO,2)+std(mSO,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(mt,2), mean(mSO,2)-std(mSO,1,2),'linewidth',1,'color',color(1,:))
     ylim([0 60])
     
    subplot(N, length(force_conditions), j+length(force_conditions)*4)
    plot(mean(mt,2), mean(mTA,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(mt,2), mean(mTA,2)+std(mTA,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(mt,2), mean(mTA,2)-std(mTA,1,2),'linewidth',1,'color',color(1,:))
    ylim([0 60])
 
        subplot(N, length(force_conditions), j+length(force_conditions)*5)
    plot(mean(mtsu,2), mean(mphi,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(mtsu,2), mean(mphi,2)+std(mphi,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(mtsu,2), mean(mphi,2)-std(mphi,1,2),'linewidth',1,'color',color(1,:))
        ylim([10 40])
        
        subplot(N, length(force_conditions), j+length(force_conditions)*6)
    plot(mean(mtsu,2), mean(mlen,2),'linewidth',2,'color',color(1,:)); hold on
    plot(mean(mtsu,2), mean(mlen,2)+std(mlen,1,2),'linewidth',1,'color',color(1,:))
    plot(mean(mtsu,2), mean(mlen,2)-std(mlen,1,2),'linewidth',1,'color',color(1,:))
        ylim([50 100])
    xlabel('Time (s)')
    
end

%%
for i = 1:42
    subplot(N, length(force_conditions), i)
    box off; 
    xlim([0 8])
end

subplot(N,6,1); ylabel('Torque (%MVC)')
subplot(N,6,7); ylabel('MG (%MVC)')
subplot(N,6,13); ylabel('LG (%MVC)')
subplot(N,6,19); ylabel('SOL (%MVC)')
subplot(N,6,25); ylabel('TA (%MVC)')
subplot(N,6,31); ylabel('Pennation (deg)')
subplot(N,6,37); ylabel('Length (mm)')
% set(gcf,'units','normalized','position',[0 .1 .6 .8])
set(gcf,'units','normalized','position',[0 0 1 .99])

