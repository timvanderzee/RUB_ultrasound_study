clear all; close all; clc
addpath(genpath('C:\Users\timvd\Documents\ultrasound-automated-algorithm'))

dates = {'3011', '1612','1601'};
i = 1;
date = dates{i};

cd(['C:\Users\timvd\OneDrive - University of Calgary\8. Ultrasound comparison - TBD\data\Test',date,'\videos']);
files = dir('*.mp4');

filename = files(1).name;
v = VideoReader(filename);

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

xpos = geofeatures(1).apo_intersect(:,1);
ypos = geofeatures(1).apo_intersect(:,2);

% next images
j = 0;
while hasFrame(v)
     j = j+1;
     disp(j)
     
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
        [xpos(:,j),ypos(:,j)] = transformPointsForward(w, xpos(:,j-1), ypos(:,j-1));

        points = detectMinEigenFeatures(Ic,'FilterSize',11, 'MinQuality', 0.005);
        points = double(points.Location);
        [inPoints] = inpolygon(points(:,1),points(:,2), APOROIx, APOROIy);
        points = points(inPoints,:);
        setPoints(pointTracker, points);
    end
    
    % timtracks deep intersection estimate, assuming superficial from flow
    fas_coef(1) = -tand(geofeatures(j).alpha);
    Mx = xpos(2,j);
    My = ypos(2,j);
    fas_coef(2) =  My - Mx * fas_coef(1);
    x_measure(1) = round(fzero(@(x) polyval(geofeatures(j).deep_coef(:) - fas_coef(:),x),0));
    y_measure(1) = polyval(geofeatures(j).deep_coef(:), x_measure(1));
    
    % low-gain drift correction for superficial
    x_measure(2) = xpos(2,1);
    y_measure(2) = polyval(geofeatures(j).super_coef(:), x_measure(2));
    
    alpha = [.5 .05];
    
    for i = 1:2
        x_apriori = xpos(i,j);
        y_apriori = ypos(i,j);

        
        x_aposteriori = x_apriori + alpha(i) * (x_measure(i) - x_apriori);
        y_aposteriori = y_apriori + alpha(i) * (y_measure(i) - y_apriori);

        xpos(i,j) = x_aposteriori;
        ypos(i,j) = y_aposteriori;
    end
%     figure(2)
%     imshow(Ic); hold on
%     plot(points(:,1),points(:,2),'r.')
% 
% %     plot(x, polyval(geofeatures.deep_coef,x),'r-')
% %     plot(x, polyval(geofeatures.super_coef,x),'r-')
%     
%     plot(xpos(:,j),ypos(:,j),'go-')
%     
%     drawnow
%     
%     hold off
    

end

%% plot time-series
for i = 1:j
    TT_faslen(i) = geofeatures(i).faslen / 564 * 5;
    TT_alpha(i) = geofeatures(i).alpha;
end

HT_faslen = sqrt(diff(xpos).^2 + diff(ypos).^2) / 564 * 5;
HT_alpha = -atan2d(diff(ypos), diff(xpos));

if ishandle(3), close(3); end; figure(3)


load('asym_high_01_Test3011.mat')


subplot(121);
plot(TT_faslen); hold on
plot(HT_faslen);

plot(Fdat.Region.FL/10)

subplot(122);
plot(TT_alpha); hold on
plot(HT_alpha);

plot(Fdat.Region.PEN/pi*180)

