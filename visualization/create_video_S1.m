clear all; close all; clc

vidObj = VideoReader('exampleVideoTA_Verheul_2022.mp4');

frames = read(vidObj);

load('Manual_Tracking_BR.mat')

id = FasData.digitizedFrames;
sid = sort(id);


ZeroPadL = 200 * ones(vidObj.Height, vidObj.Width, 'uint8');
ZeroPadR = 200 * ones(vidObj.Height, vidObj.Width, 'uint8');


%% UTT
load('exampleVideoTA_Verheul_2022_tracked_Q=0001.mat')

for i = 1:length(Fdat.Region.Fascicle.fas_x)
    fas_x(i,:) = -(Fdat.Region.Fascicle.fas_x_end{i} - vidObj.Width)  + vidObj.Width;
    fas_y(i,:) = Fdat.Region.Fascicle.fas_y_end{i};
    
    FL_UTT(i) = sqrt(diff(fas_x(i,:))^2 + diff(fas_y(i,:)).^2);
   
    
    sup_x(i,:) = -(Fdat.Region.sup_x{i}-vidObj.Width)  + vidObj.Width;
    deep_x(i,:) = -(Fdat.Region.deep_x{i}-vidObj.Width)  + vidObj.Width;
    sup_y(i,:) = Fdat.Region.sup_y{i};
    deep_y(i,:) = Fdat.Region.deep_y{i};
    
    super_apo   = [sup_x(i,:); sup_y(i,:)];
    deep_apo    = [deep_x(i,:); deep_y(i,:)];
    
    super_coef(i,:)  = polyfit(super_apo(1,:), super_apo(2,:), 1);
    deep_coef(i,:)   = polyfit(deep_apo(1,:), deep_apo(2,:), 1);
    
    FA_UTT(i) =  atan2d(abs(diff(fas_y(i,:))),abs(diff(fas_x(i,:))));
    DA_UTT(i) = atan2d(-deep_coef(i,1),1);
%     DA_UTT2(i) = atan2d(-diff(deep_apo(2,:)), diff(deep_apo(1,:)));
    
end

%% HybridTrack
load('exampleVideoTA_trackingResults.mat');
for i = 1:length(Fdat.Region.Fascicle.fas_x)
    fas_x_H(i,:) = [FPT(i).apo1int(1) Hybrid(i).apo2int(1)] + vidObj.Width;
    fas_y_H(i,:) = [FPT(i).apo1int(2) Hybrid(i).apo2int(2)];
end

%%
pixpermm =  13.6;

figure(2)
plot(FL_UT); hold on
plot(Fdat.Region.FL * 13.6,'--')


%%
figure(1)

set(gcf,'units','normalized','position',[.1 0 .5 1])

TZ = load('Manual_Tracking_TZ.mat');
BR = load('Manual_Tracking_BR.mat');
PT = load('Manual_Tracking_PT.mat');

for i = 1:length(sid)-1
    
    frameImg = rgb2gray(frames(:,:,:,sid(i)));
    cImage = [ZeroPadL, frameImg, ZeroPadR];
    
    manual_x = mean(sort([BR.FasData.pts{i}(:,1) PT.FasData.pts{i}(:,1) TZ.FasData.pts{i}(:,1)],1),2);
    manual_y = mean(sort([BR.FasData.pts{i}(:,2) PT.FasData.pts{i}(:,2) TZ.FasData.pts{i}(:,2)],1),2);

    imshow(cImage);
    line('Xdata', BR.FasData.pts{i}(:,1),'Ydata', BR.FasData.pts{i}(:,2), 'Color', [.3 .3 .3],'linewidth',1, 'marker','o','linestyle','--')
    line('Xdata', PT.FasData.pts{i}(:,1),'Ydata', PT.FasData.pts{i}(:,2), 'Color', [.5 .5 .5],'linewidth',1, 'marker','o','linestyle','--');
    line('Xdata', TZ.FasData.pts{i}(:,1),'Ydata', TZ.FasData.pts{i}(:,2), 'Color', [.7 .7 .7],'linewidth',1, 'marker','o','linestyle','--');
    
    line('Xdata', BR.ApoData.pts{i}(:,1),'Ydata', BR.ApoData.pts{i}(:,2), 'Color', [.3 .3 .3],'linewidth',1, 'marker','o','linestyle','--')
    line('Xdata', PT.ApoData.pts{i}(:,1),'Ydata', PT.ApoData.pts{i}(:,2), 'Color', [.5 .5 .5],'linewidth',1, 'marker','o','linestyle','--');
    line('Xdata', TZ.ApoData.pts{i}(:,1),'Ydata', TZ.ApoData.pts{i}(:,2), 'Color', [.7 .7 .7],'linewidth',1, 'marker','o','linestyle','--');
    
    line('Xdata',manual_x,'Ydata', manual_y, 'Color', [.5 .5 .5],'linewidth',2, 'marker','o');
    
    
    line('Xdata', fas_x(sid(i),:),'Ydata', fas_y(sid(i),:), 'Color', 'red','linewidth',2, 'marker','o');
    line('Xdata', [1 vidObj.Width*3],'Ydata', polyval(super_coef(sid(i),:), [1 vidObj.Width*3]), 'Color', 'red','linewidth',1,'linestyle','--')
    line('Xdata', [1 vidObj.Width*3],'Ydata', polyval(deep_coef(sid(i),:), [1 vidObj.Width*3]), 'Color', 'red','linewidth',1,'linestyle','--');
    
    
    line('Xdata', fas_x_H(sid(i),:),'Ydata', fas_y_H(sid(i),:), 'Color', 'blue','linewidth',2, 'marker','o');
    
    drawnow
%         pause 
    F(i) = getframe(gcf) ;
    

end

%%
FL_man = mean([BR.FasData.FLength(1:100); PT.FasData.FLength(1:100); TZ.FasData.FLength(1:100)]);
As = [BR.FasData.FAngle(1:100); PT.FasData.FAngle(1:100); TZ.FasData.FAngle(1:100)];
As(As<0) = As(As<0) + 180;
DA_man = -mean([BR.ApoData.Angle(1:100); PT.ApoData.Angle(1:100); TZ.ApoData.Angle(1:100)]);
FA_man = mean(As);
PA_man = FA_man + DA_man;

PA_UTT = FA_UTT + DA_UTT;



if ishandle(3), close(3); end
figure(3)
subplot(221)
plot(sid(1:end-1), FL_man)
hold on
plot(FL_UTT)

subplot(222)
plot(sid(1:end-1), FA_man)
hold on
plot(FA_UTT)
plot(Fdat.Region.ANG,'--')

subplot(223)
plot(sid(1:end-1), DA_man);
hold on
plot(DA_UTT)

subplot(224)
plot(sid(1:end-1), PA_man);
hold on
plot(PA_UTT)
plot(Fdat.Region.PEN,'--')

RMSE1 = sqrt(mean((FL_UTT(sid(1:end-1))' - FL_man(:)).^2,'omitnan')) / pixpermm
RMSE2 = sqrt(mean((Fdat.Region.PEN(sid(1:end-1))' - PA_man(:)).^2,'omitnan'))

RMSE2 = sqrt(mean((PA_UTT(sid(1:end-1))' - PA_man(:)).^2,'omitnan'))

%%
cd([mainfolder, 'RUB_ultrasound_study\data\ultrasound\TA_video'])

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

