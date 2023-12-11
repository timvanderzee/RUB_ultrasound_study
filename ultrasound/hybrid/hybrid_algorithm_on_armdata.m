clear all; close all; clc
addpath(genpath('C:\Users\timvd\Documents\ultrasound data\arm'))
v = VideoReader('Arm_Image06.wmv');

i = 0;
while hasFrame(v)
    i = i+1;
    disp(i)
    Is(:,:,i) = rgb2gray(readFrame(v));
end

%% do TimTrack
load('parms.mat')
parms.cut_image = 0;
% 
parms.apo.deep.cut(1) = .35;
% parms.apo.super.cut(1) = .03;

parms.apo.deep.maxangle = 20;
parms.apo.super.maxangle = 20;
parms.extrapolation = 1;

figure(1)
parms.show = 1;

[geofeatures, parms] = do_TimTrack(Is,parms);

%% save
cd('C:\Users\timvd\Documents\ultrasound data\arm')
save('Arm_Image06_geofeatures.mat');

%% cut
Is = Is(parms.ROI(2,1):parms.ROI(2,2), parms.ROI(1,1):parms.ROI(1,2),:);
apox = [1 size(Is,2)];

%% convert TimTrack to x,y points
for j = 1:length(geofeatures)
    xpos_tt(:,j) = geofeatures(j).apo_intersect(:,1);
    ypos_tt(:,j) = geofeatures(j).apo_intersect(:,2);

    sapo_tt(:,j) = polyval(geofeatures(j).super_coef,apox,1)';
    dapo_tt(:,j) = polyval(geofeatures(j).deep_coef,apox,1)';
end

%% TimTrack video
show_video(Is, xpos_tt, ypos_tt, dapo_tt, sapo_tt)
% show_video(Is, xpos_tt, ypos_tt, dapo_tt, sapo_tt,'timtrack_armvideo.gif')

%% optic flow
if ishandle(2), close(2); end; figure(2);

% fascicle intersection horizontal and vertical positions
xpos = geofeatures(1).apo_intersect(:,1);
ypos = geofeatures(1).apo_intersect(:,2);

% aponeurosis vertical positions
sapo = polyval(geofeatures(1).super_coef, apox)';
dapo = polyval(geofeatures(1).deep_coef, apox)';

% next images
% j = 0;
for j = 1:size(Is,3)
%      j = j+1;
     disp(j)
     
%     f = readFrame(v);
%     I = rgb2gray(f);
%     Ic = I(parms.ROI(2,1):parms.ROI(2,2), parms.ROI(1,1):parms.ROI(1,2));

%     Ic = rgb2gray(Is(:,:,:,j));
    Ic = Is(:,:,j);
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
        [~,sapo(:,j)] = transformPointsForward(w, apox', sapo(:,j-1));
        [~,dapo(:,j)] = transformPointsForward(w, apox', dapo(:,j-1));

        points = detectMinEigenFeatures(Ic,'FilterSize',11, 'MinQuality', 0.005);
        points = double(points.Location);
        [inPoints] = inpolygon(points(:,1),points(:,2), APOROIx, APOROIy);
        points = points(inPoints,:);
        setPoints(pointTracker, points);
       
    end
end

%% optic flow video
show_video(Is, xpos, ypos, dapo, sapo,'opticflow_armvideo.gif')

%% optic flow vs. timtrack
close all; 

figure(1)
subplot(221);
plot(xpos(1,:)); hold on
title('Horizontal deep intersect')
plot(xpos_tt(1,:)); hold on

subplot(222);
plot(xpos(2,:)); hold on
title('Horizontal superficial intersect')
plot(xpos_tt(2,:)); hold on

subplot(223);
plot(ypos(1,:)); hold on
title('Vertical deep intersect')
plot(ypos_tt(1,:)); hold on

subplot(224);
plot(ypos(2,:)); hold on
title('Vertical superficial intersect')
plot(ypos_tt(2,:)); hold on

figure(2);
subplot(221); 
plot(dapo(1,:)); hold on
plot(dapo_tt(1,:));
title('Vertical deep left intersect')

subplot(222); 
plot(sapo(1,:)); hold on
plot(sapo_tt(1,:));
title('Vertical superficial left')

subplot(223); 
plot(dapo(2,:)); hold on
plot(dapo_tt(2,:));
title('Vertical deep right')

subplot(224); 
plot(sapo(2,:)); hold on
plot(sapo_tt(2,:));
title('Vertical superficial right')

%% difference vectors
clear dxpos dypos dsapo ddapo
for i = 1:2
    dxpos(i,:) = diff(xpos(i,:));
    dypos(i,:) = diff(ypos(i,:));

    dsapo(i,:) =  diff(sapo(i,:));
    ddapo(i,:) =  diff(dapo(i,:));
end

%% state estimation vertical positions aponeuroses
alpha = .1;
sapo_est = sapo_tt;
dapo_est = dapo_tt;

for j = 1:length(dxpos)
    for i = 1:2
        s_apriori = sapo_est(i,j) + dsapo(i,j);
        d_apriori = dapo_est(i,j) + ddapo(i,j);
        
        s_aposteriori = s_apriori + alpha * (sapo_tt(i,j) - s_apriori);
        d_aposteriori = d_apriori + alpha * (dapo_tt(i,j) - d_apriori);
        
        sapo_est(i,j+1) = s_aposteriori;
        dapo_est(i,j+1) = d_aposteriori;
        
    end
end

%% show results
subplot(221); plot(dapo_est(1,:));
subplot(222); plot(sapo_est(1,:));
subplot(223); plot(dapo_est(2,:));
subplot(224); plot(sapo_est(2,:));

%% re-fit the first-orders
for j = 1:length(xpos)
    deep_coef(j,:) = polyfit(apox(:), dapo_est(:,j),1);
    super_coef(j,:) = polyfit(apox(:), sapo_est(:,j),1); 
    
    deep_coef_tt(j,:) = geofeatures(j).deep_coef;
    super_coef_tt(j,:) = geofeatures(j).super_coef;
end

figure(3)
subplot(221); plot(deep_coef_tt(:,1)); hold on
plot(deep_coef(:,1)); title('Deep coef slope')

subplot(222); plot(super_coef_tt(:,1)); hold on
plot(super_coef(:,1)); title('Super coef slope')

subplot(223); plot(deep_coef_tt(:,2)); hold on
plot(deep_coef(:,2)); title('Deep coef offset')

subplot(224); plot(super_coef_tt(:,2)); hold on
plot(super_coef(:,2)); title('Super coef offset')

%% state estimation superficial intersection point
alpha = .05;
ypos_est = ypos;
xpos_est = xpos;

for j = 1:length(dxpos)        
        % 'measurement' is value in first frame
        x_measure = xpos(2,1);

        % from optic flow
        x_apriori = xpos_est(2,j) + dxpos(2,j);
        
        % state estimation
        x_aposteriori = x_apriori + alpha * (x_measure - x_apriori);

        % store
        xpos_est(2,j+1) = x_aposteriori;
        
        % assure point is on superficial aponeurosis
        ypos_est(2,j+1) = polyval(super_coef(j+1,:), xpos_est(2,j+1));
end

figure(1)
subplot(222); plot(xpos_est(2,:)); hold on
subplot(224); plot(ypos_est(2,:)); hold on

%% state estimation deep
alpha = .1;

for j = 1:length(dxpos)
    
    % timtracks deep intersection estimate, assuming superficial point
    fas_coef(1) = -tand(geofeatures(j).alpha);
    Mx = xpos_est(2,j);
    My = ypos_est(2,j);
    fas_coef(2) =  My - Mx * fas_coef(1);
    x_measure = round(fzero(@(x) polyval(deep_coef(j,:)' - fas_coef(:),x),0));

    % from optic flow
    x_apriori = xpos_est(1,j) + dxpos(1,j);

    % state estimation
    x_aposteriori = x_apriori + alpha(1) * (x_measure - x_apriori);

    % store
    xpos_est(1,j+1) = x_aposteriori;
    
    % assure point is on deep aponeurosis
    ypos_est(1,j+1) =  polyval(deep_coef(j+1,:), xpos_est(1,j+1));

end

figure(1)
subplot(221); plot(xpos_est(1,:)); hold on
subplot(223); plot(ypos_est(1,:)); hold on

%% evaluate results
show_video(Is, xpos_est, ypos_est, dapo_est, sapo_est)

%% show video function 
function[] = show_video(Is, xpos_est, ypos_est, dapo_est, sapo_est, GIF_filename)
    [n,m,~,~] = size(Is);
    apox = [1 m];
    p = 1;
    
    close all; figure(1)
    for i = 1:length(xpos_est)-1
        faslen(i) = sqrt(diff(xpos_est(:,i)).^2 + diff(ypos_est(:,i)).^2);
    end

    subplot(5,1,5);
    plot(1:length(xpos_est)-1, faslen); 
    h = line('xdata',[],'ydata',[],'marker','.','markersize',10);
    xlabel('Sample number'); ylabel('Fascicle length (pixels)'); box off
    
    for i = 1:length(xpos_est)-1
       

        data_padded = [ones(n,round(m/2),p) Is(:,:,i) ones(n,round(m/2),p)];     % zero padding

        subplot(5,1,1:4)
        imshow(data_padded,'xdata',[-round(m/2) round(1.5*m)], 'ydata',[1 n]); hold on

        plot(xpos_est(:,i),ypos_est(:,i),'ro-','linewidth',2)
        plot(apox, dapo_est(:,i), 'g', 'linewidth',2);
        plot(apox, sapo_est(:,i), 'b', 'linewidth',2);

        hold off; drawnow
        
        subplot(5,1,5); set(h, 'xdata', i, 'ydata', faslen(i)); hold off
        
        
        if nargin > 5
            if i == 1
                gif([GIF_filename,'.gif'])
            end
            gif
        end
    end
end