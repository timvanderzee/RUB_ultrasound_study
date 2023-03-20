clear all; close all; clc

cd('C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\timtrack')

% session = 'Test1312';

dates = {'3011', '0812', '1312', '1612', '1601', '1701', '1901a', '1901b'};

conditions = {'slow_low', 'slow_high', 'medium_low', 'medium_high', 'fast_low', 'fast_high', 'asym_low', 'asym_high', ...
    'sine_020', 'sine_1020'};


for i = 1:length(conditions)

    faslen = nan(8, 2667);
    thickness = nan(8, 2667);
    phi = nan(8, 2667);
    t = nan(1, 2667);
    
    for j = 1:length(dates)
    
    if exist([conditions{i},'_geofeatures_Test',dates{j},'.mat'],'file')
        load([conditions{i},'_geofeatures_Test',dates{j},'.mat']);
    
    for k = 1:length(geofeatures)
        faslen_raw(k) = geofeatures(k).faslen;
        thickness_raw(k) = geofeatures(k).thickness;
        phi_raw(k) = geofeatures(k).phi;
    end

    t = linspace(0,80,2667);

    fs = 1/mean(diff(t));
    fc = 10;
    N = 2; Wn = fc ./ (.5*fs);
    [b,a] = butter(N, Wn, 'low');

    % filter
    thickness(j,:) = filtfilt(b,a,thickness_raw);
    phi(j,:) = filtfilt(b,a,phi_raw);
%     faslen = filtfilt(b,a,faslen);

    faslen(j,:) = thickness(j,:)./sind(phi(j,:));
    
%     t = lins

    end
    end

cd(['C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\timtrack\all_subjects'])
save([conditions{i},'_timtrack.mat'],'t','faslen','thickness','phi')
end

