clear all; close all; clc

conditions = {'slow_low', 'slow_high', 'medium_low', 'medium_high', 'fast_low', 'fast_high', 'asym_low', 'asym_high', ...
    'sine_020', 'sine_1020'};
dates = {'3011', '0812', '1312', '1612', '1601', '1701', '1901a', '1901b'};

corfac = ones(length(dates));
corfac(2) = 0.0887 / 0.065; % accident

for i = 1:length(conditions)

    faslen = nan(8, 2667);
    phi = nan(8, 2667);
    t = nan(1, 2667);
    
    for j = 1:length(dates)
    
    
    cd(['C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\ultratrack\individual_subjects\', dates{j}])
    flag = 0;
        if exist([conditions{i},'_01_Test',dates{j},'.mat'],'file')
            load([conditions{i},'_01_Test',dates{j},'.mat'])
            flag = 1;
        end
        
        if exist([conditions{i},'_02_Test',dates{j},'.mat'],'file')
            load([conditions{i},'_02_Test',dates{j},'.mat'])
            flag = 1;
        end
        
        if flag
            faslen(j,:) = Fdat.Region.FL / 10 * corfac(j); % pixels to cm
            phi(j,:) = Fdat.Region.PEN * 180/pi; % pixels to cm
        end
    end
    
cd('C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\ultratrack\all_subjects')
save([conditions{i},'_ultratrack.mat'],'t','faslen','phi')    
    
end






