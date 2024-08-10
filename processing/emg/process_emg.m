clear all; close all; clc

conditions = {'slow_low', 'slow_high', 'medium_low', 'medium_high', 'fast_low', 'fast_high', 'asym_low', 'asym_high', ...
    'sine_020', 'sine_1020','pas_005', 'pas_30','pas_120'};

datafolder = % figshare
codefolder = % github

for j = 1:length(conditions)
% 
% EMG = struct('MG.raw', nan(8, 8201),'LG.raw', nan(8, 8201),'TA.raw', nan(8, 8201),'SO.raw', nan(8, 8201));

for i = 1:8

    cd([datafolder, '\p',num2str(i),'\MAT data'])
disp(i)

% first look for version 2, then look for version 1
if exist([conditions{j},'_02.mat'],'file')
    load([conditions{j},'_02.mat'])
elseif exist([conditions{j},'_01.mat'],'file')
    % load
    load([conditions{j},'_01.mat'])
else
    disp('Did not find file')
end

if exist([conditions{j},'_01.mat'],'file') || exist([conditions{j},'_02.mat'],'file')
    % band-pass filter
    fs = 2000; fc = [30 500];
    N = 2; Wn = fc / (.5*fs);
    [b,a] = butter(N, Wn);
    
    MG_filt = filtfilt(b,a,MG.values);
    LG_filt = filtfilt(b,a,LG.values);
    TA_filt = filtfilt(b,a,TA.values);
    SO_filt = filtfilt(b,a,SOL.values);

    % low-pass filter
    fs = 2000; fc = 10;
    N = 2; Wn = fc / (.5*fs);
    [b,a] = butter(N, Wn,'low');
    
    MG_filtfilt = filtfilt(b,a,abs(MG_filt));
    LG_filtfilt = filtfilt(b,a,abs(LG_filt));
    TA_filtfilt = filtfilt(b,a,abs(TA_filt));
    SO_filtfilt = filtfilt(b,a,abs(SO_filt));
    
    % resample
    fsnew = 100; 
    tnew = 0:(1/fsnew):82;
    
    EMG.MG.raw(i,:) = interp1(MG.times, MG_filt, tnew);
    EMG.LG.raw(i,:) = interp1(LG.times, LG_filt, tnew);
    EMG.TA.raw(i,:) = interp1(TA.times, TA_filt, tnew);
    EMG.SO.raw(i,:) = interp1(SOL.times, SO_filt, tnew);
    
    EMG.MG.filt(i,:) = interp1(MG.times, MG_filtfilt, tnew);
    EMG.LG.filt(i,:) = interp1(LG.times, LG_filtfilt, tnew);
    EMG.TA.filt(i,:) = interp1(TA.times, TA_filtfilt, tnew);
    EMG.SO.filt(i,:) = interp1(SOL.times, SO_filtfilt, tnew);
    
end

end

cd([codefolder,'RUB_ultrasound_study\data\emg'])
save([conditions{j},'_EMG.mat'],'tnew','EMG')
disp(['Saved: ', conditions{j},'_EMG.mat']);

end
