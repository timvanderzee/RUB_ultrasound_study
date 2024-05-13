clear all; close all; clc

dt = .03;

%% UTT
load('exampleVideoTA_Verheul_2022.mat')

N = length(Fdat.Region.PEN);
t = 0:dt:((N-1)*dt);

for i = 1:N
    gamma(i) = Fdat.geofeatures(i).gamma;
end

Wn = 2 / (.5*30);
[b,a] = butter(2, Wn);

fgamma = filtfilt(b,a, gamma);

figure(1)
subplot(211)
plot(t, Fdat.Region.PEN*180/pi - gamma,'linewidth',2); hold on
ylim([5 35])

subplot(212)
plot(t, Fdat.Region.FL,'linewidth',2); hold on
ylim([25 100])

%% DL_track
Ldata = readmatrix('exampleVideoTA_Verheul_2022_reduced.xlsx', 'Sheet','Fasc_length');
Pdata = readmatrix('exampleVideoTA_Verheul_2022_reduced.xlsx', 'Sheet','Pennation');
Ldata = Ldata(2:end,2:end);
Pdata = Pdata(2:end,2:end);

pixpermm =  13.6;

N = size(Ldata,1);
t = 0:dt:((N-1)*dt);

for i = 1:size(Ldata,1)
    FL(i) = median(Ldata(i,:),'omitnan') / pixpermm;
    PEN(i) = median(Pdata(i,:),'omitnan');
end

subplot(211)
plot(t, PEN)

subplot(212)
plot(t, FL)

%% Hybrid Track
cd('C:\Users\u0167448\Documents\hybrid-muscle-tracking-main\TrackingResults')
load('exampleVideoTA_trackingResults.mat')

N = length(Hybrid);
for i = 1:601
    FL(i) = Hybrid(i).length;
    PEN(i) = Hybrid(i).angle;
end

t = 0:dt:((N-1)*dt);

subplot(211)
plot(t, PEN,'linewidth',2,'color',[.5 .5 .5])

subplot(212)
plot(t, FL,'linewidth',2,'color',[.5 .5 .5])


%% make nice
for i = 1:2
    subplot(2,1,i)
    xlim([0 max(t)]); hold on
    box off
end

xlabel('Time (s)'); ylabel('Length (mm)')

subplot(211)
ylabel('Pennation (deg)')

figure(1)
legend('UltraTimTrack','DL Track','location','best'); legend boxoff