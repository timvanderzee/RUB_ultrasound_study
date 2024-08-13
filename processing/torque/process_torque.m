close all;
conditions = {'slow_low', 'slow_high', 'medium_low', 'medium_high', 'fast_low', 'fast_high', 'asym_low', 'asym_high', ...
    'sine_020', 'sine_1020','pas_005', 'pas_30','pas_120'};

torque = nan(8, 8201);
angle = nan(8,8201);

for j = 1:length(conditions)
for i = 1:8

    cd([datafolder, '\dynamometry\p',num2str(i),'\MAT data'])

%     disp(i)

    % first look for version 2, then look for version 1
    if exist([conditions{j},'_03.mat'],'file')
        load([conditions{j},'_03.mat'])
    elseif exist([conditions{j},'_02.mat'],'file')
        load([conditions{j},'_02.mat'])
    elseif exist([conditions{j},'_01.mat'],'file')
        load([conditions{j},'_01.mat'])
    else
        disp('Did not find file')
    end

    if exist([conditions{j},'_01.mat'],'file') || exist([conditions{j},'_02.mat'],'file')
        % filter
        fs = 2000; fc = 10;
        N = 2; Wn = fc / (.5*fs);
        [b,a] = butter(N, Wn,'low');
        Tfilt = filtfilt(b,a, Torque.values);
        Afilt = filtfilt(b,a, Angle.values);

        % resample
        fsnew = 100; 
        tnew = 0:(1/fsnew):82;
        torque(i,:) = interp1(Torque.times, Tfilt, tnew);
        angle(i,:) = interp1(Angle.times, Afilt, tnew);

%         figure(1)
%         subplot(4,2,i)
%         plot(Torque.times, Torque.values); hold on
%         plot(Torque.times, Tfilt,'--')
%         plot(tnew, torque(i,:),'.')
% 
%         figure(2)
%         subplot(4,2,i)
%         plot(Angle.times, Angle.values); hold on
%         plot(Angle.times, Afilt,'--')
%         plot(tnew, angle(i,:),'.')

    end

end


cd([codefolder,'\data\torque'])

save([conditions{j},'_summary.mat'],'tnew','torque','angle')
% disp(['Saved: ', conditions{j},'_summary.mat']);
end
