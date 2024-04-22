clear all; clc

Qs = [nan, 0, 10.^(-4:0), 1000, inf];
foldernames = {'3011', '0812', '1312','1612','1601','1701','1901a','1901b'};
mainfolder = 'C:\Users\timvd\OneDrive - KU Leuven\8. Ultrasound comparison - TBD\UltraTimTrack_testing\';
filenames = {'*pas_005*.mp4','*pas_30*.mp4','*pas_120*.mp4'};
% filenames = {'*slow*.mp4','*medium*.mp4','*fast*.mp4','*asym*.mp4','*sine_020*.mp4','*sine_1020*.mp4'}; 

for p = 1
    foldername = foldernames{p};

    for j = 5 %1:length(Qs)
            Qval = Qs(j);
            
        for k = 3 %1:length(filenames)

            % get names
            cd([mainfolder foldername]);
            files = dir(filenames{k});
            vidname = files.name(1:end-4);

            filename = [vidname,'_analyzed_Q=',strrep(num2str(Qval),'.','')];

            cd([mainfolder foldername,'\analyzed\mat']);
            disp([mainfolder foldername,'\analyzed\mat'])
            disp(filename);
            
            if exist([filename,'.mat'],'file')
                load([filename,'.mat']);
            else
                disp(['Does not exist: ', filename])
            end


            %% recreate a priori estimates
            x_apriori(1) = Fdat.Region.Fascicle.alpha{1};
            x_apost(1)=  Fdat.Region.Fascicle.alpha{1};
            alpha_prev = Fdat.Region.Fascicle.alpha{1};
            xflow(1) = Fdat.Region.Fascicle.alpha{1};
            y(1) = Fdat.Region.Fascicle.alpha{1};

            xkal(1) = Fdat.Region.Fascicle.alpha{1};

            P_apost = 0;
            P_apriori = P_apost;
            
            N = length(Fdat.Region.FL);
            
            if isfinite(Qval)
            
            for i = 1:N
                X(i) = Fdat.Region.Fascicle.alpha{i};
            end

            for i = 2:N

                fas_prev = [Fdat.Region.Fascicle.fas_x{i-1}' Fdat.Region.Fascicle.fas_y{i-1}'];
                alpha_prev = Fdat.Region.Fascicle.alpha{i-1};

                % transform
                w = Fdat.Region.warp(:,:,i-1);
                fas_new = transformPointsForward(w, fas_prev);

                % Estimate the change in fascicle angle from the change in points
                dalpha = abs(atan2d(diff(fas_new(:,2)), diff(fas_new(:,1)))) - abs(atan2d(diff(fas_prev(:,2)), diff(fas_prev(:,1))));

                % a priori state estimate
                x_minus = alpha_prev + dalpha;

                % estimate covariances
                dx = abs(dalpha);
                dx(dx<0.005) = 0;
                f.Q = getQ(Qval, dx);
                f.R = Fdat.R(1);

                % apriori estimate from optical flow
                f.x_minus = x_minus;

                % measurement from Hough transform
                f.y = Fdat.geofeatures(i).alpha;

                % previous state covariance
                f.P_prev = P_apost(i-1);

                % run kalman filter
                F = run_kalman_filter(f);

                P_apriori(i) = F.P_minus;

                % ascribe
                alpha = F.x_plus;
                Kgain(i) = F.K;
                P_apost(i) = F.P_plus;

                xflow(i) = xflow(i-1) + dalpha;
                x_apriori(i) = x_minus;

                x_apost(i) = alpha;
                y(i) = Fdat.geofeatures(i).alpha;
            end


            %% smoothing                      
            for i = (N-1):-1:1
                xprev = x_apriori(i+1);
                Pprev = P_apriori(i+1);
                Pnew = P_apost(i);

                C(i) = Pnew / Pprev;
                C(isnan(C)) = 1;
                
                Fdat.Region.Fascicle.alpha{i} = x_apost(i) + C(i) * (Fdat.Region.Fascicle.alpha{i+1}  - xprev);
            
                fasx2 = Fdat.Region.Fascicle.fas_x{i}(2);
                fasy2 = Fdat.Region.Fascicle.fas_y{i}(2);
                
                % fit the current aponeurosis
                ROI         = [Fdat.Region.ROIx{i} Fdat.Region.ROIy{i}];
                super_apo   = ROI([1,4],:);
                deep_apo    = ROI([2,3],:);
                super_coef  = polyfit(super_apo(:,1), super_apo(:,2), 1);
                deep_coef   = polyfit(deep_apo(:,1), deep_apo(:,2), 1);
                
                % get the deep attachment point from the superficial point and the angle
                fas_coef(1) = -tand(Fdat.Region.Fascicle.alpha{i});
                fas_coef(2) =  fasy2 - fas_coef(1) * fasx2;
                fasx1 = (fas_coef(2) - deep_coef(2)) / (deep_coef(1) - fas_coef(1));
                fasy1 = deep_coef(2) + fasx1*deep_coef(1);

                %% update
                % state and dependent variables
                Fdat.Region.Fascicle.fas_x{i}   = [fasx1 fasx2];
                Fdat.Region.Fascicle.fas_y{i}   = [fasy1 fasy2];
                
            end
            end
            
            for i = 1:N
                Fdat.Region.Fascicle.fas_p_minus{i} = P_apriori(i);
                Fdat.Region.Fascicle.X_minus{i} = x_apriori(i);
            end
            
            faspen = nan(1,N);
            faslen = nan(1,N);
            
            
            ImDepth = 50.7; % mm

            if p < 6
                vidHeight = 564;
            else
                vidHeight = 780;
            end

            scalar = ImDepth/vidHeight;

            for i = 1:N
                % calculate the length and pennation for the current frame
                faspen(i) = atan2(abs(diff(Fdat.Region.Fascicle.fas_y{i})),...
                    abs(diff(Fdat.Region.Fascicle.fas_x{i})));
                
                faslen(i) = scalar*sqrt(diff(Fdat.Region.Fascicle.fas_y{i}).^2 +...
                    diff(Fdat.Region.Fascicle.fas_x{i}).^2);
            end
            
            Fdat.Region.FL = faslen;
            Fdat.Region.PEN = faspen;
            
            cd([mainfolder foldername,'\analyzed\mat']);
            save([filename,'_v2.mat'],'Fdat','TrackingData');
            
            
        end
    end
end


%%
function [K] = run_kalman_filter(k)
% this assumes we already have the aposteriori state estimate (k.x_minus),
% the measurement (k.y) and the process- and measurement noise covariances (k.R and k.Qvalue)

% a posteriori variance estimate
K.P_minus = k.P_prev + k.Q;

% kalman xshiftcor
K.K = K.P_minus / (K.P_minus + k.R);

if (K.K < 0) || (K.K > 1)
    disp('warning: kalman gain outside 0-1 interval');

    K.K(K.K<0) = 0;
    K.K(K.K>1) = 1;
end

% estimated state
K.x_plus = k.x_minus + K.K * (k.y - k.x_minus);

% estimated variance
K.P_plus = (1-K.K) * K.P_minus;
end

function[Q] = getQ(Qval, dx)

% Optical flow error is proportional to flow
Q = Qval  * dx^2;
end
