close all; 

v = 10; % number of contraction
msdphis = nan(v, 8, 6, 3,2);
msdlens = nan(v, 8, 6, 3,2);

mtsus = nan(101, 8, 6, 3,2);
mphis = nan(101, 8, 6, 3,2);
mlens = nan(101, 8, 6, 3,2);

for p = 1:8
%     disp(['p',num2str(p)])
for j = 1:2
    if j == 1
        filenames = {'*slow_low*','*medium_low*','*fast_low*','*asym_low*'}; 
        
    else
        filenames = {'*slow_high*','*medium_high*','*fast_high*','*asym_high*'}; 
%         fs = 100/3;
    end
    
for k = 1:length(filenames)
    cd([codefolder, '\data\ultrasound\GM\p', num2str(p)])
    files = dir(filenames{k});

    for i = 1:length(files)
        filename = files(i).name;

        if exist(filename,'file')
            load(filename,'Fdat');
        else
            disp(['Not found: ', filename,'.mat'])
        end
        
        if length(Fdat.Region.FL) > 2000
            fs = 100/3;
        else
            fs = 25;
        end

        T = 8;
        Ts = (0:10)*T;
        dt = 1/fs;

        if j == 1
            tall = 0:dt:(Ts(end)-dt);

        else
            tall = 0:dt:Ts(end);
        end

        msd_phi = nan(length(Ts)-1,1);
        msd_faslen = nan(length(Ts)-1,1);
        pen_rs = nan(length(Ts)-1, 101);
        len_rs = nan(length(Ts)-1, 101);

        for n = 1:length(Ts)-1
            id = (tall >= Ts(n)) & (tall <= Ts(n+1));

            % cut
            pen = Fdat.Region.PEN(id);
            len = Fdat.Region.FL(id);
            tnew = mod(tall(id),T);

            % unique
            [tnew,IA,IC] = unique(tnew);
            upen = pen(IA);
            ulen = len(IA);

            % sort
            [tnew2, is] = sort(tnew);

            % resample
            ts = linspace(0, T, 101);
            pen_rs(n,:) = interp1(tnew2, upen(is), ts,[],'extrap');
            len_rs(n,:) = interp1(tnew2, ulen(is), ts,[],'extrap');    

            % mean and standard deviation
            phi = mean(pen_rs,'omitnan');
            sd_phi = std(pen_rs,1,'omitnan');
            faslen = mean(len_rs,'omitnan');
            sd_faslen = std(len_rs,1,'omitnan');

            msd_phi(n,1) = mean(sd_phi,'omitnan');
            msd_faslen(n,1) = mean(sd_faslen,'omitnan');
        end

        x2 = 0:100;

        % save
        mtsus(:,p,k,i,j) = ts;
        mphis(:,p,k,i,j) = phi;
        mlens(:,p,k,i,j) = faslen;

        msdphis(1:length(Ts)-1,p,k,i,j) = msd_phi;
        msdlens(1:length(Ts)-1,p,k,i,j) = msd_faslen;

    end
end
end
end

%% save
cd([codefolder, '\data\ultrasound\tracking'])
save('ramp_summary.mat','msdphis','mtsus','mphis','mlens','msdlens')
