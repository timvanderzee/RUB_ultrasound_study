clear all; clc

addpath(genpath('C:\Users\timvd\Documents\UltraTimTrack'))
cd('C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\3011')

example_video = 'slow_high_01.mp4';
v = VideoReader(example_video);

%% read
frames = read(v, [1 500]);
cutframes = frames(44:606, 233:939,:,:);

figure(1)
imshow(cutframes(:,:,:,1))

%%
id = 1:30:size(cutframes,4);
frames_ds = cutframes(:,:,:,id);

cd('C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing');

vidObj = VideoWriter('the_5_frames_video.mp4',   'MPEG-4'  );
open(vidObj);

for i = 1:10
    figure(10)
    subplot(2,5,i)
    imshow(frames_ds(:,:,:,i));
    
    writeVideo(vidObj,frames_ds(:,:,:,i));
    
    
end

   close(vidObj);

%% save as image
cd('C:\Users\timvd\Documents\RUB_ultrasound_study\10_frames')
iptsetpref('ImshowBorder','tight');
if ishandle(10), close(10); end
for i = 1:10
    figure(10)
    % subplot(2,5,i)
    imshow(frames_ds(:,:,:,i));

    saveas(gcf,['frame_',num2str(i),'.jpg'])

end

%% do TimTrack
load('TimTrack_parms.mat')
parms.fas.show = 1;
parms.show = 1;
if ishandle(11), close(11); end
figure(11)
for i = 1:10
    handles.geofeatures(i) = auto_ultrasound(rgb2gray(frames_ds(:,:,:,i)),parms);
end

%% get ROI
handles.NumFrames = 10;
handles.vidWidth = size(frames_ds,2);
[handles] = lowpass_ROI(handles);

%% do Ultratrack
frames = cutframes;
n = id(5);

for i = n
    im1 = rgb2gray(frames(:,:,:,i));
    im2 = rgb2gray(frames(:,:,:,i+1));
    show_UltraTrack(im1,im2, handles)
    drawnow
end
iptsetpref('ImshowBorder','tight');

%%


function[handles] = lowpass_ROI(handles)
    i = 1;
    n = handles.vidWidth;

    frames = 1:handles.NumFrames;
    APO_TT = nan(5,length(frames));
    
    for f = frames
        handles.Region(i).ROIy{f} = round([polyval(handles.geofeatures(f).super_coef, 1) polyval(handles.geofeatures(f).deep_coef, [1 n]) polyval(handles.geofeatures(f).super_coef, [n 1])])';
        handles.Region(i).ROIx{f} = [1 1 n n 1]';
            
        APO_TT(:,f) = handles.Region(i).ROIy{f};
    end
        
    y = APO_TT;

    for f = frames
        handles.Region.ROIy{f} = round(y(:,f));
    end
end

function[handles] = show_UltraTrack(im1,im2, handles)

    %% Optical flow and state estimation
    % setup current and new image
    frame_no = 1;
    i = 1;
            
    points = detectMinEigenFeatures(im1,'FilterSize',11, 'MinQuality', 0.005);
    points = double(points.Location);
    [inPoints] = inpolygon(points(:,1),points(:,2), handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no});
    points = points(inPoints,:);

    %% Initialize point tracker and run opticflow
    %calculate block size according to ROI 
    width       = floor(max(abs(diff(handles.Region.ROIx{frame_no}))) * 0.0725); %0.20
    height      = floor(max(abs(diff(handles.Region.ROIy{frame_no}))) * 0.08); %0.40; thickness changes?

    handles.BlockSize = [width height]; %save as width and height for later comparison

    % Ensure vidWidth and vidHeight are both odd numbers
    if mod(width, 2) == 0
        width = width + 1; % Increment by 1 to make it odd
    end

    if mod(height, 2) == 0
        height = height + 1; % Increment by 1 to make it odd
    end

    pointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',[height width]);
    initialize(pointTracker,points,im1);

    % Compute the flow and new roi
    [pointsNew, isFound] = step(pointTracker, im2);

    figure(20)
    imshow(im1); hold on
    quiver(points(isFound,1),points(isFound,2), (pointsNew(isFound,1)-points(isFound,1)),(pointsNew(isFound,2)-points(isFound,2)), 2,'linewidth',1);
    hold off
    
        
end
