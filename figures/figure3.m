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
    
    
    x1 = linspace(0,100, length(mean(ts,2)))';

    subplot(N, length(force_conditions), j)
    plot(x1, mean(Ts,2),'linewidth',2,'color',color(1,:)); hold on
   
    coord_combine = [[x1 mean(Ts,2)+std(Ts,1,2)]; flipud([x1 mean(Ts,2)-std(Ts,1,2)])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');
    
%     plot(mean(ts,2), mean(Ts,2)+std(Ts,1,2),'linewidth',1,'color',color(1,:))
%     plot(mean(ts,2), mean(Ts,2)-std(Ts,1,2),'linewidth',1,'color',color(1,:))
    ylim([-5 120])

    title(titles{j})

    subplot(N, length(force_conditions), j+length(force_conditions)*1)
    plot(x1, mean(MGs,2),'linewidth',2,'color',color(1,:)); hold on
    
    coord_combine = [[x1 mean(MGs,2)+std(MGs,1,2)]; flipud([x1 mean(MGs,2)-std(MGs,1,2)])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');
    
%     plot(mean(ts,2), mean(MGs,2)+std(MGs,1,2),'linewidth',1,'color',color(1,:))
%     plot(mean(ts,2), mean(MGs,2)-std(MGs,1,2),'linewidth',1,'color',color(1,:))
    ylim([0 60])

    
    subplot(N, length(force_conditions), j+length(force_conditions)*2)
    plot(x1, mean(LGs,2),'linewidth',2,'color',color(1,:)); hold on
    
    coord_combine = [[x1 mean(LGs,2)+std(LGs,1,2)]; flipud([x1 mean(LGs,2)-std(LGs,1,2)])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');
        
%     plot(mean(ts,2), mean(LGs,2)+std(LGs,1,2),'linewidth',1,'color',color(1,:))
%     plot(mean(ts,2), mean(LGs,2)-std(LGs,1,2),'linewidth',1,'color',color(1,:))
    ylim([0 60])

    subplot(N, length(force_conditions), j+length(force_conditions)*3)
    plot(x1, mean(SOs,2),'linewidth',2,'color',color(1,:)); hold on
    
    coord_combine = [[x1 mean(SOs,2)+std(SOs,1,2)]; flipud([x1 mean(SOs,2)-std(SOs,1,2)])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');
    
%     plot(mean(ts,2), mean(SOs,2)+std(SOs,1,2),'linewidth',1,'color',color(1,:))
%     plot(mean(ts,2), mean(SOs,2)-std(SOs,1,2),'linewidth',1,'color',color(1,:))
    ylim([0 60])

    subplot(N, length(force_conditions), j+length(force_conditions)*4)
    plot(x1, mean(TAs,2),'linewidth',2,'color',color(1,:)); hold on
    
    coord_combine = [[x1 mean(TAs,2)+std(TAs,1,2)]; flipud([x1 mean(TAs,2)-std(TAs,1,2)])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');
    
%     plot(mean(ts,2), mean(TAs,2)+std(TAs,1,2),'linewidth',1,'color',color(1,:))
%     plot(mean(ts,2), mean(TAs,2)-std(TAs,1,2),'linewidth',1,'color',color(1,:))
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

mainfolder = 'C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
% mainfolder = 'C:\Users\u0167448\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
subfolders = dir(mainfolder);

foldernames = {'3011', '0812', '1312','1612','1601','1701','1901a','1901b'};
filenames = {'*slow*.mp4','*medium*.mp4','*fast*.mp4','*asym*.mp4','*sine_020*.mp4','*sine_1020*.mp4'}; 

% force_conditions = {'slow','medium','fast','asym','sine_020','sine_1020'};

% fignames = {'sine 0-20', 'sine 0-10', 'passive 5 deg/s','passive 30 deg/s','passive 120 deg/s','ramp: asymmetric','ramp: slow','ramp: medium','ramp: fast'};
participants = foldernames;

dcolor = [color(2,:); .5 .5 .5; color(1,:)];

ix = [length(Qs) 1 5];
% ix = 1;
% ix = 5;
% ix = 6;

v = 120; % number of contraction
msdphis = nan(v, 8, 6, length(Qs));
msdlens = nan(v, 8, 6, length(Qs));

mtsus = nan(101, 8, 6, length(Qs));
mphis = nan(101, 8, 6, length(Qs));
mlens = nan(101, 8, 6, length(Qs));

for p = 1:8
    m = 0;
    foldername = foldernames{p};
    
    disp(foldername)
    
for i = 1:length(Qs)
    if ismember(i,ix)
        m = m+1;
    end
    
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
    
    if k < 5
        f = 0.125;
        T = 8;
        Ts = (0:10)*T;
    else
        f = 1.5;
        T = 1/1.5;
        Ts = (0:v)*T;
    end
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
        tnew = mod(tall(id),1/f);
        
%         figure(100)
%         plot(tnew, pen); hold on

        % sort
        [tnew2, is] = sort(tnew);

        % resample
        ts = linspace(0, T, 101);
        pen_rs(n,:) = interp1(tnew2, pen(is), ts,[],'extrap');
        len_rs(n,:) = interp1(tnew2, len(is), ts,[],'extrap');
% 
%         figure(100)
%         plot(tnew2, pen(is), '-', ts, pen_rs(n,:),'.'); 
%         pause
        


        % moving average
%         phi = movmean(pen(is), 100);
%         faslen = movmean(len(is), 100);

        % moving s.d.
%         sd_phi = movstd(pen(is), 100);
%         sd_faslen = movstd(len(is), 100);
        
%         msd_phi(n,1) = mean(sd_phi);
%         msd_faslen(n,1) = mean(sd_faslen);
%         keyboard
    
    
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
        subplot(N, length(force_conditions), k+length(force_conditions)*5)
    %     plot(ts, Fdat.Region.PEN(is)*180/pi,'.','color',[.5 .5 .5 .5]); hold on
        plot(x2, phi, 'linewidth',2,'color',dcolor(m,:)); hold on

        coord_combine = [[x2;phi+sd_phi] fliplr([x2;phi-sd_phi])];

        h = fill(coord_combine(1,:),coord_combine(2,:),'b');
        set(h,'FaceColor',dcolor(m,:),'FaceAlpha',.3,'LineStyle','none');

        ylim([10 40])

        subplot(N, length(force_conditions), k+length(force_conditions)*6)  
    %     plot(ts, Fdat.Region.FL(is),'.','color',[.5 .5 .5 .5]); hold on
        plot(x2, faslen, '-','linewidth',2,'color',dcolor(m,:)); hold on
        coord_combine = [[x2;faslen+sd_faslen] fliplr([x2;faslen-sd_faslen])];
              h = fill(coord_combine(1,:),coord_combine(2,:),'b');
        set(h,'FaceColor',dcolor(m,:),'FaceAlpha',.3,'LineStyle','none');

        ylim([50 100])
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
end

%%
for i = 1:42
    subplot(N, length(force_conditions), i)
    box off; 
    xlim([0 100])
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
    
    subplot(N, length(force_conditions), j)
    plot(x1, mean(mT,2),'linewidth',2,'color',color(1,:)); hold on
    
    coord_combine = [[x1 mean(mT,2)+std(mT,1,2)]; flipud([x1 mean(mT,2)-std(mT,1,2)])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');
    
%     plot(mean(mt,2), mean(mT,2)+std(mT,1,2),'linewidth',1,'color',color(1,:))
%     plot(mean(mt,2), mean(mT,2)-std(mT,1,2),'linewidth',1,'color',color(1,:))
    title(titles{j})
     ylim([0 60])
     
    subplot(N, length(force_conditions), j+length(force_conditions)*1)
    plot(x1, mean(mMG,2),'linewidth',2,'color',color(1,:)); hold on
    
    coord_combine = [[x1 mean(mMG,2)+std(mMG,1,2)]; flipud([x1 mean(mMG,2)-std(mMG,1,2)])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');
    
%     plot(mean(mt,2), mean(mMG,2)+std(mMG,1,2),'linewidth',1,'color',color(1,:))
%     plot(mean(mt,2), mean(mMG,2)-std(mMG,1,2),'linewidth',1,'color',color(1,:))
     ylim([0 60])
     
    subplot(N, length(force_conditions), j+length(force_conditions)*2)
    plot(x1, mean(mLG,2),'linewidth',2,'color',color(1,:)); hold on
        
    coord_combine = [[x1 mean(mLG,2)+std(mLG,1,2)]; flipud([x1 mean(mLG,2)-std(mLG,1,2)])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');
    
%     plot(mean(mt,2), mean(mLG,2)+std(mLG,1,2),'linewidth',1,'color',color(1,:))
%     plot(mean(mt,2), mean(mLG,2)-std(mLG,1,2),'linewidth',1,'color',color(1,:))
     ylim([0 60])
     
    subplot(N, length(force_conditions), j+length(force_conditions)*3)
    plot(x1, mean(mSO,2),'linewidth',2,'color',color(1,:)); hold on
        
    coord_combine = [[x1 mean(mSO,2)+std(mSO,1,2)]; flipud([x1 mean(mSO,2)-std(mSO,1,2)])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');
    
%     plot(mean(mt,2), mean(mSO,2)+std(mSO,1,2),'linewidth',1,'color',color(1,:))
%     plot(mean(mt,2), mean(mSO,2)-std(mSO,1,2),'linewidth',1,'color',color(1,:))
     ylim([0 60])
     
    subplot(N, length(force_conditions), j+length(force_conditions)*4)
    plot(x1, mean(mTA,2),'linewidth',2,'color',color(1,:)); hold on
    
    coord_combine = [[x1 mean(mTA,2)+std(mTA,1,2)]; flipud([x1 mean(mTA,2)-std(mTA,1,2)])];
    h = fill(coord_combine(:,1),coord_combine(:,2),'b');
    set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');
    
%     plot(mean(mt,2), mean(mTA,2)+std(mTA,1,2),'linewidth',1,'color',color(1,:))
%     plot(mean(mt,2), mean(mTA,2)-std(mTA,1,2),'linewidth',1,'color',color(1,:))
    ylim([0 60])
 

    for m = 1:3
        i = ix(m);
        mtsu = mtsus(:,:,j,i);
        mphi = mphis(:,:,j,i) - mphis(1,:,j,i) +  mean(mphis(1,:,j,i),2);
        mlen = mlens(:,:,j,i) - mlens(1,:,j,i) +  mean(mlens(1,:,j,i),2);
    
        subplot(N, length(force_conditions), j+length(force_conditions)*5)
        plot(x2(:), mean(mphi,2),'linewidth',2,'color',dcolor(m,:)); hold on

        coord_combine = [[x2(:) mean(mphi,2)+std(mphi,1,2)]; flipud([x2(:) mean(mphi,2)-std(mphi,1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolor(m,:),'FaceAlpha',.3,'LineStyle','none');

    %     plot(mean(mtsu,2), mean(mphi,2)+std(mphi,1,2),'linewidth',1,'color',color(1,:))
    %     plot(mean(mtsu,2), mean(mphi,2)-std(mphi,1,2),'linewidth',1,'color',color(1,:))
            ylim([10 40])

            subplot(N, length(force_conditions), j+length(force_conditions)*6)
        plot(x2(:), mean(mlen,2),'linewidth',2,'color',dcolor(m,:)); hold on

        coord_combine = [[x2(:) mean(mlen,2)+std(mlen,1,2)]; flipud([x2(:) mean(mlen,2)-std(mlen,1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolor(m,:),'FaceAlpha',.3,'LineStyle','none');

    %     plot(mean(mtsu,2), mean(mlen,2)+std(mlen,1,2),'linewidth',1,'color',color(1,:))
    %     plot(mean(mtsu,2), mean(mlen,2)-std(mlen,1,2),'linewidth',1,'color',color(1,:))
            ylim([50 100])
        xlabel('Contraction cycle (%)')
    end
    
end

%%
for i = 1:42
    subplot(N, length(force_conditions), i)
    box off; 
    xlim([0 100])
end

subplot(N,6,1); ylabel('Torque (%MVC)')
subplot(N,6,7); ylabel('MG (%MVC)')
subplot(N,6,13); ylabel('LG (%MVC)')
subplot(N,6,19); ylabel('SOL (%MVC)')
subplot(N,6,25); ylabel('TA (%MVC)')
subplot(N,6,31); ylabel('Pennation (deg)')
subplot(N,6,37); ylabel('Length (mm)')
% set(gcf,'units','normalized','position',[0 .1 .6 .8])
set(gcf,'units','normalized','position',[0 0 .5 .99])

%% plot the variability with respect to the cycle-average
if ishandle(100), close(100); end

dcolors = dcolor + .1;

figure(100)
for j = 1:length(force_conditions)
   
    if j < 5
        n = 1:10;
    else
        n = 1:120;
    end

    for m = [2, 1, 3]
        i = ix(m);
        
        subplot(2,6,j);
        coord_combine = [[n(:) mean(msdphis(n,:,j,i),2)+std(msdphis(n,:,j,i),1,2)]; flipud([n(:) mean(msdphis(n,:,j,i),2)-std(msdphis(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
        ylim([0 5])
        box off; hold on
        plot(mean(msdphis(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(msdphis(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
        title(titles{j})
        xlim([min(n) max(n)])
        
        subplot(2,6,j+6);
        coord_combine = [[n(:) mean(msdlens(n,:,j,i),2)+std(msdlens(n,:,j,i),1,2)]; flipud([n(:) mean(msdlens(n,:,j,i),2)-std(msdlens(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
        ylim([0 10])
        box off; hold on
        plot(mean(msdlens(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(msdlens(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
        
        xlabel('Cycle #');
        
        xlim([min(n) max(n)])
        
    end
    
    m = 2;
    i = ix(m);

    subplot(2,6,j);
    plot(mean(msdphis(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
     plot(max(n), mean(msdphis(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on

    subplot(2,6,j+6);
    plot(mean(msdlens(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
     plot(max(n), mean(msdlens(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
end

subplot(261); ylabel('Angle deviation (deg)');
subplot(267); ylabel('Length deviation (mm)');

%%
set(gcf,'units','normalized','position',[.2 .2 .5 .5])

%% final variability wrt Q value
if ishandle(200), close(200); end
figure(200)

X = cell(1,length(Qs));

for i = 1:length(Qs)
    X{i} = num2str(Qs(i));
end

x = reordercats(categorical(X),X);

for j = 1:length(force_conditions)
    if j < 5
        n = 1:10;
    else
        n = 1:120;
    end
       
        subplot(2,6,j);
        bar(x,     squeeze(mean(msdphis(length(n),:,j,:),2)),'EdgeColor',color(1,:),'BarWidth',0.7); hold on
        errorbar(1:length(x), squeeze(mean(msdphis(length(n),:,j,:),2)), squeeze(std(msdphis(length(n),:,j,:),1,2)),'marker','none','linestyle','none','color',color(1,:),'CapSize',2);

        ylim([0 12])
        box off; hold on
        title(titles{j})

%         xlim([0 length(x)])
        
        subplot(2,6,j+6);
        bar(x, squeeze(mean(msdlens(length(n),:,j,:),2)),'EdgeColor',color(1,:),'BarWidth',0.7); hold on
        errorbar(1:length(x),   squeeze(mean(msdlens(length(n),:,j,:),2)), squeeze(std(msdlens(length(n),:,j,:),1,2)),'marker','none','linestyle','none','color',color(1,:),'CapSize',2);
        
        ylim([0 20])
        box off; hold on
       
        xlabel('Q');
        
%         xlim([0 length(x)])
        
end

%%
subplot(261); ylabel('Angle deviation (deg)');
subplot(267); ylabel('Length deviation (mm)');

set(gcf,'units','normalized','position',[.2 .2 .5 .5])