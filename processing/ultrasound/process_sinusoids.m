close all; 
filenames = {'*sine_1020*','*sine_020*'}; 

T = 1/1.5;
Ts = (1:120)*T;
tall = 0:.03:Ts(end);

msdphis = nan(length(Ts)-1, 8, 6, 3);
msdlens = nan(length(Ts)-1, 8, 6, 3);

noise_phi = nan(length(Ts)-2, 8, 6, 3);
drift_phi = nan(length(Ts)-2, 8, 6, 3);
noise_len = nan(length(Ts)-2, 8, 6, 3);
drift_len = nan(length(Ts)-2, 8, 6, 3);

mtsus = nan(101, 8, 6, 3);
mphis = nan(101, 8, 6, 3);
mlens = nan(101, 8, 6, 3);

for p = 1:8
%     disp(['p',num2str(p)])
    for k = 1:length(filenames)
        cd([codefolder, '\data\ultrasound\GM\p', num2str(p)])
        files = dir(filenames{k});

        for i = 1:length(files)
            filename = files(i).name;

            if exist(filename,'file')
                load(filename);

                msd_phi = nan(length(Ts)-1,1);
                msd_faslen = nan(length(Ts)-1,1);
                pen_rs = nan(length(Ts)-1, 101);
                len_rs = nan(length(Ts)-1, 101);

                % low-pass
                fs = 1/.03;
                fc = 2.5;
                N = 2;
                Wn = fc / (.5*fs);
                [b,a] = butter(N, Wn, 'high');
                hpen = abs(filtfilt(b,a,Fdat.Region.PEN(:)*180/pi));

                for n = 1:length(Ts)-1
                    id = (tall >= Ts(n)) & (tall <= Ts(n+1));

                    % cut
                    pen = Fdat.Region.PEN(id);
                    len = Fdat.Region.FL(id);
                    tnew = mod(tall(id),T); 

                    % sort
                    [tnew2, is] = sort(tnew);

                    % resample
                    ts = linspace(0, T, 101);
                    pen_rs(n,:) = interp1(tnew2, pen(is), ts,[],'extrap');
                    len_rs(n,:) = interp1(tnew2, len(is), ts,[],'extrap');    

                    % mean and standard deviation
                    phi = mean(pen_rs,'omitnan');
                    sd_phi = std(pen_rs,1,'omitnan');
                    faslen = mean(len_rs,'omitnan');
                    sd_faslen = std(len_rs,1,'omitnan');

                    msd_phi(n,1) = mean(sd_phi,'omitnan');
                    msd_faslen(n,1) = mean(sd_faslen,'omitnan');

                end

                % save
                mtsus(:,p,k,i) = ts;
                mphis(:,p,k,i) = phi;
                mlens(:,p,k,i) = faslen;

                % drift estimate
                drift_phi(:,p,k,i) = abs(cumsum(mean(diff(pen_rs),2)));
                drift_len(:,p,k,i) = abs(cumsum(mean(diff(len_rs),2)));

                % noise estimte
                noise_phi(:,p,k,i) = std(diff(pen_rs),1,2);
                noise_len(:,p,k,i) = std(diff(len_rs),1,2);

                % deviation estimates
                msdphis(:,p,k,i) = msd_phi;
                msdlens(:,p,k,i) = msd_faslen;

            else
                disp(['Not found: ', filename,'.mat'])
            end

        end
    end
end


%% save
cd([codefolder, '\data\ultrasound\tracking'])
save('sines_summary.mat')
