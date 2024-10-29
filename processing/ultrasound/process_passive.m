close all; 
force_conditions = {'pas_005', 'pas_30','pas_120'};

fs = 100;
n = 2667;
tus = 0:.03:((n-1)*.03);

angle_rs = nan(n, 8, 3);

for j = 1:length(force_conditions)
    
    % load
    load([force_conditions{j},'_summary.mat'],'angle', 'tnew')

    % down-sample
    for p = 1:8
        angle_rs(: ,p,j) = interp1(tnew, angle(p,:)', tus,[],'extrap');
    end
   
    % get max angle
    max_angle(:,j) = max(angle,[],2);
    min_angle(:,j) = min(angle,[],2);
end

% mean(max_angle)

%% ultrasound
filenames = {'*pas_005*','*pas_30*','*pas_120*'};

nn = 2667;
pen = nan(nn, 8, 3, 3);
len = nan(nn, 8, 3, 3);

for p = 1:8
    
    for k = 1:length(filenames)
        cd([datafolder, '\fascicle_tracking_estimates\p', num2str(p)])
        files = dir(filenames{k});

        for i = 1:length(files)
            filename = files(i).name;

            if exist(filename,'file')
                load(filename);

                n = min(2667, length(Fdat.Region.PEN));

                % save
                pen(1:n,p,k,i) = Fdat.Region.PEN(1:n);
                len(1:n,p,k,i) = Fdat.Region.FL(1:n);
            else
                disp('Does not exist')
            end
            
%             figure(100+p);
%             subplot(1,3,k);
%             plot(angle_rs(:,p,k), pen(:,p,k,i),'.'); hold on
%             


        end
    end
end


%% bin
mpen = nan(50,8,3,3);
spen = nan(50,8,3,3);
mlen = nan(50,8,3,3);
slen = nan(50,8,3,3);

for p = 1:8
    for a = 1:3
        for k = 1:3

                as = floor(min(angle_rs(:,p,k))):1:(ceil(max(angle_rs(:,p,k)))-1);
            for i = 1:length(as)-1
                id = (angle_rs(:,p,k) > as(i)) & (angle_rs(:,p,k) < as(i+1));

                mpen(i,p,k,a) = mean(pen(id,p,k,a));
                spen(i,p,k,a) = std(pen(id,p,k,a));

                mlen(i,p,k,a) = mean(len(id,p,k,a));
                slen(i,p,k,a) = std(len(id,p,k,a));
            end

        end
    end
end


%% average over angles
for p = 1:8
    for a = 1:3
        for k = 1:3
            mspen(a,p,k) = mean(spen(:,p,k,a),'omitnan');
            mslen(a,p,k) = mean(slen(:,p,k,a),'omitnan');
        end
    end
end

%% plot
figure(100)
for a = 1:3
    errorbar(1:3, squeeze(mean(mspen(a,:,:),2)), squeeze(std(mspen(a,:,:),1,2))); hold on
end

%% save
if ~isfolder([codefolder, '\data\ultrasound\tracking'])
    mkdir([codefolder, '\data\ultrasound\tracking'])
end

cd([codefolder, '\data\ultrasound\tracking'])
save('passive_summary.mat')
