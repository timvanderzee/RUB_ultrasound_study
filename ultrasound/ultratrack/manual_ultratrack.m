clear all; close all; clc
addpath(genpath('C:\Users\timvd\Documents\ultrasound-automated-algorithm'))

dates = {'3011', '1612','1601'};
i = 1;
date = dates{i};

cd(['C:\Users\timvd\OneDrive - University of Calgary\8. Ultrasound comparison - TBD\data\Test',date,'\videos']);
files = dir('*.mp4');

filename = files(1).name;
v = VideoReader(filename);
% frames = read(v);
%%
% load('parms.mat')
% % parms.cut_image = 0;
% 
% parms.ROI = [239   936; 50   500]; % [0812]
% 
% parms.apo.deep.cut(1) = .35;
% parms.apo.super.cut(1) = .03;
% 
% % figure(1)
% [geofeatures, parms] = do_TimTrack(I,parms);

load(['asym_high_geofeatures_Test',date,'.mat'])


%% find points within ROI

if ishandle(2), close(2); end; figure(2);

xpos = nan(2,500,2);
ypos = nan(2,500,2);

xpos(:,1,1) = geofeatures(1).apo_intersect(:,1);
ypos(:,1,1) = geofeatures(1).apo_intersect(:,2);

xpos(:,1,2) = geofeatures(500).apo_intersect(:,1);
ypos(:,1,2) = geofeatures(500).apo_intersect(:,2);

% next images
j = 0;

id = nan(2,500);
id(1,:) = 1:500;
id(2,:) = 500:-1:1;

k = 1;
% for k = 1:2
for j = 1:500
     disp(j)
     
%     f = frames(:,:,:,id(k,j));
    f = readFrame(v);
    I = rgb2gray(f);
    Ic = I(parms.ROI(2,1):parms.ROI(2,2), parms.ROI(1,1):parms.ROI(1,2));
    x = 1:size(Ic,2);
    
    % optional: fix to edges
    APOROI(:,1) = [1; 1; size(Ic,2); size(Ic,2); 1];

    % make sure ROIs are on aponeuroses
    APOROI([2,3],2) = polyval(geofeatures(j).deep_coef, APOROI([2,3],1));
    APOROI([1,4,5],2) = polyval(geofeatures(j).super_coef, APOROI([1,4,5],1));

    [M, APOROIx,APOROIy] = roipoly(Ic, APOROI(:,1),APOROI(:,2));

    if j == 1
        % first image
        points = detectMinEigenFeatures(Ic,'FilterSize',11, 'MinQuality', 0.005);
        points = double(points.Location);

        [inPoints] = inpolygon(points(:,1),points(:,2), APOROIx, APOROIy);
        points = points(inPoints,:);

        pointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',[21 71]);
        initialize(pointTracker,points,Ic);   

    else

        [pointsNew, isFound] = step(pointTracker, Ic);
        [w,~] = estimateGeometricTransform2D(points(isFound,:), pointsNew(isFound,:), 'affine', 'MaxDistance',50);
        [xpos(:,j,k),ypos(:,j,k)] = transformPointsForward(w, xpos(:,j-1,k), ypos(:,j-1,k));

        points = detectMinEigenFeatures(Ic,'FilterSize',11, 'MinQuality', 0.005);
        points = double(points.Location);
        [inPoints] = inpolygon(points(:,1),points(:,2), APOROIx, APOROIy);
        points = points(inPoints,:);
        setPoints(pointTracker, points);
    end
% end
end

%%

if ishandle(10), close(10); end
x = xpos(1,:,1);

for i = 1:500
    xcor(i) = x(i) - 0.5 * (x(i) - x(1));
end

figure(10);
subplot(121); 
plot(x); hold on
plot(xcor)

xpos_cor = xpos;
xpos_cor(1,:,1) = xcor;

faslen = sqrt(diff(xpos(:,:,k)).^2 + diff(ypos(:,:,k)).^2) / 564 * 5;
faslen_cor = sqrt(diff(xpos_cor(:,:,k)).^2 + diff(ypos(:,:,k)).^2) / 564 * 5;

subplot(122);
plot(faslen); hold on
plot(faslen_cor)

%%

for k = 1:2
faslen(:,k) = sqrt(diff(xpos(:,:,k)).^2 + diff(ypos(:,:,k)).^2) / 564 * 5;
 
end

L = [faslen(:,1) flip(faslen(:,2))];

if ishandle(10), close(10); end
figure(10)
plot(L); hold on

w = [(500:-1:1)' (1:500)'];

Lw = sum(w .* L,2) / 501;
plot(Lw)


