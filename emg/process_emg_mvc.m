clear all; close all; clc

dates = {'3011', '0812', '1312', '1612', '1601', '1701', '1901a', '1901b'};
datafolder = 'C:\Users\timvd\OneDrive - University of Calgary\8. Ultrasound comparison - TBD\data';

j = 1;

EMG = struct('MG', nan(8, 8201),'LG', nan(8, 8201),'TA', nan(8, 8201),'SO', nan(8, 8201));

MVCnames = {'MVC', 'MVCdf'};

MVC = struct('MG',nan(8,2),'LG',nan(8,2),'SO',nan(8,2),'TA',nan(8,2));

for k = 1:2
for i = 1:length(dates)

cd([datafolder, '\Test', dates{i}])
disp(i)

% assume 5 possible version
nfiles = 0;

for m = 1:5
    if exist([MVCnames{k},'_0', num2str(m),'.mat'],'file')
        nfiles = nfiles + 1;

        load([MVCnames{k},'_0', num2str(m),'.mat'])
        disp([MVCnames{k},'_0', num2str(m),'.mat'])
    end
end

if nfiles > 0
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
    
    EMG.MG(i,:) = interp1(MG.times, MG_filtfilt, tnew);
    EMG.LG(i,:) = interp1(LG.times, LG_filtfilt, tnew);
    EMG.TA(i,:) = interp1(TA.times, TA_filtfilt, tnew);
    EMG.SO(i,:) = interp1(SOL.times, SO_filtfilt, tnew);

    [MVC.MG(i,k), MGloc] = max(EMG.MG(i,:));
    [MVC.LG(i,k), LGloc] = max(EMG.LG(i,:));
    [MVC.TA(i,k), TAloc] = max(EMG.TA(i,:));
    [MVC.SO(i,k), SOloc] = max(EMG.SO(i,:));
    
    figure(k*10+1)
    subplot(4,2,i)
    plot(MG.times, MG.values); hold on
    plot(MG.times, abs(MG_filt),'--')
    plot(MG.times, MG_filtfilt,'--')
    plot(tnew, EMG.MG(i,:),'.')
    plot(tnew(MGloc), MVC.MG(i,k),'ro')
    ylim([-2 2])
    set(gcf,'name','MG')
    
    figure(k*10+2)
    subplot(4,2,i)
    plot(LG.times, LG.values); hold on
    plot(LG.times, LG_filtfilt,'--')
    plot(tnew, EMG.LG(i,:),'.')
    plot(tnew(LGloc), MVC.LG(i,k),'ro')
    ylim([-2 2])
    set(gcf,'name','LG')
    
    figure(k*10+3)
    subplot(4,2,i)
    plot(SOL.times, SOL.values); hold on
    plot(SOL.times, SO_filtfilt,'--')
    plot(tnew, EMG.SO(i,:),'.')
    plot(tnew(SOloc), MVC.SO(i,k),'ro')
    ylim([-2 2])
    set(gcf,'name','SOL')
    
    figure(k*10+4)
    subplot(4,2,i)
    plot(TA.times, TA.values); hold on
    plot(TA.times, TA_filtfilt,'--')
    plot(tnew, EMG.TA(i,:),'.')
    plot(tnew(TAloc), MVC.TA(i,k),'ro')
    ylim([-2 2])
    set(gcf,'name','TA')
end
end

end

cd('C:\Users\timvd\Documents\RUB_ultrasound_study\emg\summary_data')
save('MVC_EMG.mat','tnew','MVC')

