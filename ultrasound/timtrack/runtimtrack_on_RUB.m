clear all; close all; clc
addpath(genpath('C:\Users\timvd\Documents\ultrasound-automated-algorithm'))

dates = {'3011', '1612','1601'};
for i = 2:3
date = dates{i};

cd(['C:\Users\timvd\OneDrive - University of Calgary\8. Ultrasound comparison - TBD\data\Test',date,'\videos']);
files = dir('*.mp4');

%%
filename = files(1).name;
v = VideoReader(filename);
f = readFrame(v);
I = rgb2gray(f);
% 
load('parms.mat')
% parms.cut_image = 0;

parms.ROI = [239   936; 50   500]; % [0812]

parms.apo.deep.cut(1) = .35;
parms.apo.super.cut(1) = .03;

% figure(1)
[geofeatures, parms] = do_TimTrack(I,parms);


%%
k = 0;
parms.show = 0;
parms.fas.show = 0;

for j = 1 %2:length(files)
    clear geofeatures
    
    cd(['C:\Users\timvd\OneDrive - University of Calgary\8. Ultrasound comparison - TBD\data\Test',date,'\videos']);
    filename = files(j).name;
    disp(filename)
    
    % read video
    v = VideoReader(filename);
    k = 0;
    
    while hasFrame(v)
%         for i = 1:9
%             f = readFrame(v);
%         end

        k = k+1;
        disp(k)

        f = readFrame(v);
        I = rgb2gray(f);

        %profile on
        Icut = I(parms.ROI(2,1):parms.ROI(2,2), parms.ROI(1,1):parms.ROI(1,2));

        [geofeatures(k), ~, parms] = auto_ultrasound(Icut,parms);
        
        cd(['C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\timtrack\individual_subjects\', date])
        save([filename(1:end-7), '_geofeatures_Test',date,'.mat'], 'geofeatures','parms')
    %     profile viewer
    end
end
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