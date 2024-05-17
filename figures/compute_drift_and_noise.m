clear all; close all; clc
% typical example participant, torques, angles, EMGs

Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('RampTarget.mat','tnew','rampTarget')
load('MVC_EMG.mat');

%% time series
force_conditions = {'sine_020','sine_1020'};
image_qualities = {'low', 'high'};
i = 2;

titles = {'0-20 %MVC sine','10-20 %MVC sine'};

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

    T = 1/1.5;
    Ts = (0:120)*T;
    ts = linspace(0, T, 101);
    
    for n = 1:length(Ts)-1
        id = (tcut >= Ts(n)) & (tcut <= Ts(n+1));
        tmod = mod(tcut(id),T);
        
        T_pc(n,:) = interp1(tmod, Tcut(id), ts,[], 'extrap');
        MG_pc(n,:) = interp1(tmod, MGrel_cut(id), ts,[], 'extrap');
        LG_pc(n,:) = interp1(tmod, LGrel_cut(id), ts,[], 'extrap');
        TA_pc(n,:) = interp1(tmod, TArel_cut(id), ts,[], 'extrap');
        SO_pc(n,:) = interp1(tmod, SOrel_cut(id), ts,[], 'extrap');
    end
    
    mu =    [mean(T_pc); mean(MG_pc); mean(LG_pc); mean(SO_pc); mean(TA_pc)];
    sigma = [std(T_pc); std(MG_pc); std(LG_pc); std(SO_pc); std(TA_pc)];
        
%     figure(p)
%     color = get(gca,'colororder');
    
%     for kk = 1:size(mu,1)
    
%         x1 = 0:100;
% 
%         subplot(N, M, j+M*(kk-1))
%         plot(x1, mu(kk,:),'linewidth',2,'color',color(1,:)); hold on
% 
%         coord_combine = [[x1; mu(kk,:)+sigma(kk,:)], fliplr([x1; mu(kk,:)-sigma(kk,:)])];
%         h = fill(coord_combine(1,:),coord_combine(2,:),'b');
%         set(h,'FaceColor',color(1,:),'FaceAlpha',.3,'LineStyle','none');

%         ylim([-2 40])
%     end
    
    % save the cycle average
    mus(:,:,p,j) = mu;
end
end

%% add ultrasound
N = 7;
Qs = [nan, 0, 10.^(-4:0), 1000, inf];
color = get(gca,'colororder');

% mainfolder = 'C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
mainfolder = 'C:\Users\u0167448\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
subfolders = dir(mainfolder);

foldernames = {'3011', '0812', '1312','1612','1601','1701','1901a','1901b'};
filenames = {'*sine_020*.mp4','*sine_1020*.mp4'}; 

% force_conditions = {'slow','medium','fast','asym','sine_020','sine_1020'};

% fignames = {'sine 0-20', 'sine 0-10', 'passive 5 deg/s','passive 30 deg/s','passive 120 deg/s','ramp: asymmetric','ramp: slow','ramp: medium','ramp: fast'};
participants = foldernames;

%% Figure 3: cycle average
dcolor = [.5 .5 .5; color(2,:)+[0 .2 .2]; color(1,:)];

ix = [length(Qs) 1 5];
% ix = 1;
% ix = 5;
% ix = 6;

iis = [1:4, 6:length(Qs), 5];

% iis = length(Qs)

T = 1/1.5;
Ts = (1:120)*T;
tall = 0:.03:Ts(end);

msdphis = nan(length(Ts)-1, 8, 6, length(Qs));
msdlens = nan(length(Ts)-1, 8, 6, length(Qs));

noise_phi = nan(length(Ts)-1, 8, 6, length(Qs));
drift_phi = nan(length(Ts)-1, 8, 6, length(Qs));
noise_len = nan(length(Ts)-1, 8, 6, length(Qs));
drift_len = nan(length(Ts)-1, 8, 6, length(Qs));

mtsus = nan(101, 8, 6, length(Qs));
mphis = nan(101, 8, 6, length(Qs));
mlens = nan(101, 8, 6, length(Qs));

for p = 1:8
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
        end

    
%     SNR_phi(p,k,i) = snr(Fdat.Region.PEN*180/pi);
%     SNR_len(p,k,i) = snr(Fdat.Region.FL);
    
    fs = 100/3;
    
    T = 1/1.5;
    Ts = (1:120)*T;
    tall = 0:.03:Ts(end);

    id = tall > T;
    
    S = Fdat.Region.PEN(id)*180/pi - mean(Fdat.Region.PEN(id)*180/pi);
    
    Y = fft(S);
    L = length(S);
    f1 = fs*(0:(L/2))/L;
    
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

%     plot(f1,P1) 
%     title('Single-Sided Amplitude Spectrum of S(t)')
%     xlabel('f (Hz)')
%     ylabel('|P1(f)|')

    
%     [P1,f1] = pspectrum(Fdat.Region.PEN*180/pi - mean(Fdat.Region.PEN*180/pi),fs);
%     [P2,f2] = pspectrum(Fdat.Region.FL - mean(Fdat.Region.FL),fs);
    
%     figure(p)
%     subplot(1,2,k)
%     plot(f1, P1); hold on
    
    SRN_phi(p,k,i) = sum(P1(f1 > 2));
    DFT_phi(p,k,i) = sum(P1(f1 < 0.5));
%     DFT_phi(p,k,i) = mean(S((end-266):end)) - mean(S(1:267));
%     SRN_len(p,k,i) = sum(P2(f2 > 2));

end
end
end

%%
figure(10)
subplot(121)
semilogx(Qs, squeeze(SRN_phi(:,1,:))); hold on
semilogx(Qs, squeeze(mean(SRN_phi(:,1,:))),'o')

subplot(122)
semilogx(Qs, squeeze(SRN_phi(:,2,:))); hold on
semilogx(Qs, squeeze(mean(SRN_phi(:,2,:))),'o')

%%
figure(11)
subplot(121)
semilogx(Qs, squeeze(DFT_phi(:,1,:))); hold on
semilogx(Qs, squeeze(mean(DFT_phi(:,1,:))),'o')

subplot(122)
semilogx(Qs, squeeze(DFT_phi(:,2,:))); hold on
semilogx(Qs, squeeze(mean(DFT_phi(:,2,:))),'o')