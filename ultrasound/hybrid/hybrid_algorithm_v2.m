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
for s = 1
    
    %% load timtrack
    cd(['C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\timtrack\individual_subjects\', foldernames{s}(5:end)])
    load([filename, '_timtrack_',foldernames{s},'.mat'], 'geofeatures','parms')

    apox = [1 parms.ROI(1,2)-parms.ROI(1,1)];
    
    % convert TimTrack to x,y points
    for j = 1:length(geofeatures)
        xpos_tt(:,j) = geofeatures(j).apo_intersect(:,1);
        ypos_tt(:,j) = geofeatures(j).apo_intersect(:,2);

        sapo_tt(:,j) = polyval(geofeatures(j).super_coef,apox,1)';
        dapo_tt(:,j) = polyval(geofeatures(j).deep_coef,apox,1)';
    end

    %% load optic flow
   cd(['C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\hybrid\individual_subjects\', foldernames{s}(5:end)])
   load([filename, '_opticflow_with_timtrack_',foldernames{s},'.mat'], 'geofeatures','parms','xpos', 'ypos', 'dapo', 'sapo')

    %% optic flow vs. timtrack
    figure(1)
    set(gcf,'name','Fascicle-aponeurosis intersection locations')
    subplot(221);
    plot(xpos(1,:)); hold on
    title('Horizontal deep intersection')
    plot(xpos_tt(1,:)); hold on

    subplot(222);
    plot(xpos(2,:)); hold on
    title('Horizontal superficial intersection')
    plot(xpos_tt(2,:)); hold on

    subplot(223);
    plot(ypos(1,:)); hold on
    title('Vertical deep intersection')
    plot(ypos_tt(1,:)); hold on

    subplot(224);
    plot(ypos(2,:)); hold on
    title('Vertical superficial intersection')
    plot(ypos_tt(2,:)); hold on

    figure(2);
    set(gcf,'name','Aponeurosis locations (at image borders)')
    subplot(221); 
    plot(dapo(1,:)); hold on
    plot(dapo_tt(1,:));
    title('Vertical deep left intersection')

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

    for j = 1:2
        figure(j)
    for i = 1:4
        subplot(2,2,i)
        box off
        xlabel('Sample #')
        ylabel('Position (pixels)')
    end

    legend('Optic flow','TimTrack','location','best')
    legend box off
    end

    %% calculate difference vectors
    clear dxpos dypos dsapo ddapo
    for i = 1:2
        dxpos(i,:) = diff(xpos(i,:));
        dypos(i,:) = diff(ypos(i,:));

        dsapo(i,:) =  diff(sapo(i,:));
        ddapo(i,:) =  diff(dapo(i,:));
    end

    %% state estimation vertical positions aponeuroses
    % do state estimation on vertical position

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
    % use estimates to recompute aponeurosis fits
    for j = 1:length(xpos)
        deep_coef(j,:) = polyfit(apox(:), dapo_est(:,j),1);
        super_coef(j,:) = polyfit(apox(:), sapo_est(:,j),1); 

        deep_coef_tt(j,:) = geofeatures(j).deep_coef;
        super_coef_tt(j,:) = geofeatures(j).super_coef;
    end

    if ishandle(3), close(3); end; figure(3)
    subplot(221); plot(deep_coef_tt(:,1)); hold on
    plot(deep_coef(:,1)); title('Deep aponeurosis slope')

    subplot(222); plot(super_coef_tt(:,1)); hold on
    plot(super_coef(:,1)); title('Superficial aponeurosis slope')

    subplot(223); plot(deep_coef_tt(:,2)); hold on
    plot(deep_coef(:,2)); title('Deep aponeurosis offset')

    subplot(224); plot(super_coef_tt(:,2)); hold on
    plot(super_coef(:,2)); title('Superficial aponeurosis offset')


    for i = 1:4
        subplot(2,2,i)
        box off
        xlabel('Sample #')
        ylabel('Position (pixels)')
    end

    legend('TimTrack','State estimation','location','best')
    legend box off


    %% state estimation superficial aponeurosis intersection point
    % state estimation for horizontal point
    % get vertical point from aponeurosis

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

    %% state estimation deep aponeurosis intersection point
    % state estimation for horizontal point
    % get vertical point from aponeurosis

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

    %% create video
    cd([datafolder,foldernames{s},'\videos'])
    
    if exist([filename,'_02.mp4'])
        fullfilename = [filename,'_02.mp4'];
    else
        fullfilename = [filename,'_01.mp4'];
    end
    
    v = VideoReader(fullfilename);
    
    % next images
    i = 0;
    
    
    figure(100)
    set(gcf,'position', [0 0 800 850])
    
    while hasFrame(v)
        i = i+1; 
        disp(i)
        
        f = readFrame(v);
        I = rgb2gray(f);
        Ic = I(parms.ROI(2,1):parms.ROI(2,2), parms.ROI(1,1):parms.ROI(1,2));
    
        [n,m] = size(Ic);
            
        p = 1;
            
        data_padded = [ones(n,round(m/2),p) Ic ones(n,round(m/2),p)];     % zero padding
        
        figure(100)
        subplot(311); imshow(data_padded,'xdata',[-round(m/2) round(1.5*m)], 'ydata',[1 n]); hold on
        plot(xpos_tt(:,i),ypos_tt(:,i),'ro-','linewidth',2)
        plot(apox, dapo_tt(:,i), 'g', 'linewidth',2);
        plot(apox, sapo_tt(:,i), 'b', 'linewidth',2);
        hold off; drawnow; title('TimTrack')
         
        subplot(312); imshow(data_padded,'xdata',[-round(m/2) round(1.5*m)], 'ydata',[1 n]); hold on
        plot(xpos(:,i),ypos(:,i),'ro-','linewidth',2)
        plot(apox, dapo(:,i), 'g', 'linewidth',2);
        plot(apox, sapo(:,i), 'b', 'linewidth',2);
        hold off; drawnow;  title('Optic flow')
         
        subplot(313); imshow(data_padded,'xdata',[-round(m/2) round(1.5*m)], 'ydata',[1 n]); hold on
        plot(xpos_est(:,i),ypos_est(:,i),'ro-','linewidth',2)
        plot(apox, dapo_est(:,i), 'g', 'linewidth',2);
        plot(apox, sapo_est(:,i), 'b', 'linewidth',2);
        hold off; drawnow;  title('State estimator')
        
        if i == 1
            gif([filename,'.gif'])
        end
        gif

    
    end
end

%% evaluate results
% show_video(Is, xpos_est, ypos_est, dapo_est, sapo_est,[vidname 'estimator_video.gif'])


%% functions
function[] = show_video(Is, xpos_est, ypos_est, dapo_est, sapo_est, GIF_filename)
    [n,m,~,~] = size(Is);
    apox = [1 m];
    p = 1;
    
    close all; figure(1)


    for i = 1:length(xpos_est)-1

        data_padded = [ones(n,round(m/2),p) rgb2gray(Is(:,:,:,i)) ones(n,round(m/2),p)];     % zero padding
        imshow(data_padded,'xdata',[-round(m/2) round(1.5*m)], 'ydata',[1 n]); hold on

        plot(xpos_est(:,i),ypos_est(:,i),'ro-','linewidth',2)
        plot(apox, dapo_est(:,i), 'g', 'linewidth',2);
        plot(apox, sapo_est(:,i), 'b', 'linewidth',2);

        hold off; drawnow
        
        if nargin > 5
            if i == 1
                gif([GIF_filename,'.gif'])
            end
            gif
        end
    end
end