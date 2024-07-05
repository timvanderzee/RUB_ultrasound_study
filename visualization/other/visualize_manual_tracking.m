clear all; close all; clc

vidObj = VideoReader('exampleVideoTA_Verheul_2022.mp4');

frames = read(vidObj);


load('Manual_Tracking_BR.mat')

id = FasData.digitizedFrames;
sid = sort(id);


ZeroPadL = 200 * ones(vidObj.Height, vidObj.Width, 'uint8');
ZeroPadR = 200 * ones(vidObj.Height, vidObj.Width, 'uint8');

%%
figure(1)

set(gcf,'units','normalized','position',[.1 0 .5 1])
% fas = line('Color', 'red');
TZ = load('Manual_Tracking_TZ.mat');
BR = load('Manual_Tracking_BR.mat');
PT = load('Manual_Tracking_PT.mat');

for i = 1:length(sid)
    
    frameImg = rgb2gray(frames(:,:,:,sid(i)));
    cImage = [ZeroPadL, frameImg, ZeroPadR];
    
    subplot(311)
    imshow(cImage);
    line('Xdata', round(BR.FasData.pts{i}(:,1)),'Ydata', round(BR.FasData.pts{i}(:,2)), 'Color', 'red','linewidth',2);
    title('Brent')
    
    subplot(312)
    imshow(cImage);    
    line('Xdata', round(PT.FasData.pts{i}(:,1)),'Ydata', round(PT.FasData.pts{i}(:,2)), 'Color', 'blue','linewidth',2);
    title('Paolo')
    
    subplot(313)
    imshow(cImage);
    line('Xdata', round(TZ.FasData.pts{i}(:,1)),'Ydata', round(TZ.FasData.pts{i}(:,2)), 'Color', 'green','linewidth',2);
    title('Tim')
    
    drawnow
    
  F(i) = getframe(gcf) ;
end

%%
% create the video writer with 1 fps
writerObj = VideoWriter('Manual_tracking.mp4','MPEG-4');
writerObj.FrameRate = vidObj.FrameRate/6;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
% convert the image to a frame
frame = F(i) ;    
writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

