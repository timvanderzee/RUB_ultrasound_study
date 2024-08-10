clear all; close all; clc
load('MVC_EMG.mat');


%% ramp conditions
force_conditions = {'slow','medium','fast','asym'};
image_qualities = {'low', 'high'};

for i = 1:2
for j = 1:4
    
% load
load([force_conditions{j},'_',image_qualities{i},'_EMG.mat'])

% subtract rest torque, divide by MVC torque, multiply by 100%
TArel = EMG.TA ./ max(MVC.TA,[],2) * 100;

for k = 1:8
    figure(k)
    subplot(2,2,j);

%     plot(tnew, EMG.TA(k,:)); hold on
    plot(tnew, TArel(k,:)); hold on
    title(force_conditions{j})
    
end

end
end


