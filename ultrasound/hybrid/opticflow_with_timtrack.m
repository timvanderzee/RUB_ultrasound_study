clear all; close all; clc
datafolder = 'C:\Users\timvd\OneDrive - University of Calgary\8. Ultrasound comparison - TBD\data\';
TimTrackfolder = 'C:\Users\timvd\Documents\ultrasound-automated-algorithm\';
codefolder = 'C:\Users\timvd\Documents\RUB_ultrasound_study\';

addpath(genpath(TimTrackfolder))
addpath(genpath(codefolder))
filename = 'sine_020';

%% get folder names
folder = dir(datafolder);

for i = 1:(length(folder)-2)
    foldernames{i} = folder(i+2).name;
end

%%
for i = 1 %:length(foldernames)
    cd(['C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\timtrack\individual_subjects\', foldernames{i}(5:end)])
    load([filename, '_timtrack_',foldernames{i},'.mat'], 'geofeatures','parms')
    
    %% load and video
    cd([datafolder,foldernames{i},'\videos'])
    
    if exist([filename,'_02.mp4'])
        fullfilename = [filename,'_02.mp4'];
    else
        fullfilename = [filename,'_01.mp4'];
    end
    
    v = VideoReader(fullfilename);
    
    % next images
    j = 0;
    while hasFrame(v)
        j = j+1; 
        disp(j)
        
        f = readFrame(v);
        I = rgb2gray(f);
        Ic = I(parms.ROI(2,1):parms.ROI(2,2), parms.ROI(1,1):parms.ROI(1,2));
        
        x = 1:size(Ic,2);
        apox = [1 size(Ic,2)];
        
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

            % fascicle intersection horizontal and vertical positions
            xpos = geofeatures(1).apo_intersect(:,1);
            ypos = geofeatures(1).apo_intersect(:,2);

            % aponeurosis vertical positions
            sapo = polyval(geofeatures(1).super_coef, apox)';
            dapo = polyval(geofeatures(1).deep_coef, apox)';
            
        else

            [pointsNew, isFound] = step(pointTracker, Ic);
            [w,~] = estimateGeometricTransform2D(points(isFound,:), pointsNew(isFound,:), 'affine', 'MaxDistance',50);
            
%             rot_est(j) = asind((w.T(1,2) - w.T(2,1)) / 2); % deg

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
    
    cd(['C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\hybrid\individual_subjects\', foldernames{i}(5:end)])
    save([filename, '_opticflow_with_timtrack_',foldernames{i},'.mat'], 'geofeatures','parms','xpos', 'ypos', 'dapo', 'sapo')
    
end

