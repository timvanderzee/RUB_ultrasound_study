clear all; close all; clc
datafolder = 'C:\Users\timvd\OneDrive - University of Calgary\8. Ultrasound comparison - TBD\data\';
codefolder = 'C:\Users\timvd\Documents\ultrasound-automated-algorithm\';

addpath(genpath(codefolder))
filename = 'sine_020';

% TimTrack parameters
load('parms.mat')
% parms.cut_image = 0;

parms.ROI = [239   936; 50   608]; % [0812]
parms.extrapolation = 0;
parms.fas.thetares = .5;

parms.apo.deep.cut(1) = .35;
parms.apo.super.cut(1) = .03;


%% get folder names
folder = dir(datafolder);

for i = 1:(length(folder)-2)
    foldernames{i} = folder(i+2).name;
end


%% first frame: test settings
parms.show = 1;
parms.fas.show = 1;
parms.redo_ROI = 0;

for i = 1:length(foldernames)
    cd([datafolder,foldernames{i},'\videos'])
    
    if exist([filename,'_02.mp4'])
        fullfilename = [filename,'_02.mp4'];
    else
        fullfilename = [filename,'_01.mp4'];
    end

 
    v = VideoReader(fullfilename);
    f = readFrame(v);

    I = rgb2gray(f);
    
    if i == 5 || i == 6 || i == 7
       ROI(:,:,i) = [1 size(I,2); 1 size(I,1)];
    else
       ROI(:,:,i) = [239   936; 44   608];
    end
    
    if i == 7 
        ROI(2,end,i) = 608;
    else
        ROI(2,end,i) = 500;
    end
    
    parms.ROI = ROI(:,:,i);

    %profile on
    Icut = I(parms.ROI(2,1):parms.ROI(2,2), parms.ROI(1,1):parms.ROI(1,2));

    subplot(4,2,i)
    [geofeatures, ~, parms] = auto_ultrasound(Icut,parms);
end

%% all frames
parms.show = 0;
parms.fas.show = 0;
parms.redo_ROI = 0;

for i = 8 %1:(length(foldernames)-1)
    cd([datafolder,foldernames{i},'\videos'])
    
    if exist([filename,'_02.mp4'])
        fullfilename = [filename,'_02.mp4'];
    else
        fullfilename = [filename,'_01.mp4'];
    end
    
    parms.ROI = ROI(:,:,i);

    % read video
    v = VideoReader(fullfilename);
    k = 0;

    clear geofeatures
    while hasFrame(v)

        k = k+1;
        disp(k)

        f = readFrame(v);
        I = rgb2gray(f);

        Icut = I(parms.ROI(2,1):parms.ROI(2,2), parms.ROI(1,1):parms.ROI(1,2));

        [geofeatures(k), ~, parms] = auto_ultrasound(Icut,parms);

    end

    cd(['C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\timtrack\individual_subjects\', foldernames{i}(5:end)])
    save([filename, '_timtrack_',foldernames{i},'.mat'], 'geofeatures','parms')
end

return
%% extract fascicle length and pennation
for i = 1:length(geofeatures)
    faslen(i) = geofeatures(i).faslen;
    phi(i) = geofeatures(i).phi;
    thickness(i) = geofeatures(i).thickness;
    
end

fs = 3.333;
fc = .5;
N = 2; Wn = fc ./ (.5*fs);
[b,a] = butter(N, Wn, 'low');
% filter
thickness_filt = filtfilt(b,a,thickness);
phi_filt = filtfilt(b,a,phi);
faslen_filt = filtfilt(b,a,faslen);

faslen_filt2 = thickness_filt./sind(phi_filt);

if ishandle(2), close(2); end; figure(2)
subplot(221);
plot(phi); hold on
plot(phi_filt)

subplot(222);
plot(thickness); hold on
plot(thickness_filt)

subplot(223);
plot(faslen); hold on
plot(faslen_filt); 
plot(faslen_filt2,'--');