clear all; close all; clc
Qs = [nan, 0, 10.^(-4:0), 1000, inf];

mainfolder = 'C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
% mainfolder = 'C:\Users\u0167448\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
subfolders = dir(mainfolder);

foldernames = {'3011', '0812', '1312','1612','1601','1701','1901a','1901b'};
filenames = {'*slow*.mp4','*medium*.mp4','*fast*.mp4','*asym*.mp4','*sine_020*.mp4','*sine_1020*.mp4'}; 

ix = [length(Qs) 1 5];

pen0 = nan(length(filenames), 3, 8);
len0 = nan(length(filenames), 3, 8);

for p = 1:8
    m = 0;
    foldername = foldernames{p};
    disp(foldername)
    
for i = ix
    m = m+1;
    
for k = 1:length(filenames)
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

    % get the first length and pennation
    pen0(k,m,p) = mean(Fdat.Region.PEN(1))*180/pi;
    len0(k,m,p) = mean(Fdat.Region.FL(1));
    
%     figure(p)
%     subplot(1,3,m)
%     plot(Fdat.Region.PEN(1:20)*180/pi); hold on
end
end
end

%%
if ishandle(10), close(10); end
figure(10)
for p = 1:8
    subplot(121);
    plot(pen0(:,3,p)); hold on
    
    subplot(122);
    plot(len0(:,3,p)); hold on
end

%%
cd('C:\Users\timvd\Documents\RUB_ultrasound_study\figures')
save('L0.mat','pen0','len0')
