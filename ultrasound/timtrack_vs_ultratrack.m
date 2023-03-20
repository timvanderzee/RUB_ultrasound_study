clear all; close all; clc


% ramp conditions
conditions = {'slow_high', 'medium_high', 'fast_high', 'asym_high', ...
    'sine_020', 'sine_1020'};
% image_qualities = {'low', 'high'};
q = 2;

for j = 1:length(conditions)
    
    
% load
TTdata = load([conditions{j},'_timtrack.mat']);

TTdata.faslen = TTdata.faslen / 564 * 5; % pixels to cm
TTdata.thickness = TTdata.thickness / 564 * 5; % pixels to cm

UTdata = load([conditions{j},'_ultratrack.mat']);


for i = 1:5
    
    % simple state-estimation
    dUTphi = diff(UTdata.phi(i,:));
    alpha = .1;
    SEphi(1) = TTdata.phi(i,1);

    for ii = 1:length(TTdata.phi(i,:))-1

        phi_apriori = SEphi(ii) + dUTphi(ii);

        phi_aposteriori = phi_apriori + alpha * (TTdata.phi(i,ii) - phi_apriori);

        SEphi(1,ii+1) = phi_aposteriori;
    end
    
    SEdata.phi(i,:) = SEphi;
    SEdata.faslen(i,:) = TTdata.thickness(i,:) ./ sind(SEdata.phi(i,:));
    
    
    figure(i)

    subplot(2,6,j); 
    plot(TTdata.t, TTdata.faslen(i,:)); hold on
    plot(TTdata.t, UTdata.faslen(i,:));  
    plot(TTdata.t, SEdata.faslen(i,:));  
    title(conditions{j}(1:4))
    xlim([0 80]); box off; ylabel('Length (cm)')
    
    subplot(2,6,j+6); 
    plot(TTdata.t, TTdata.phi(i,:)); hold on
    plot(TTdata.t, UTdata.phi(i,:)); 
    plot(TTdata.t, SEdata.phi(i,:)); 
    ylabel('Pennation (deg)')
    xlim([0 80]); box off
    
    figure(100+i)
    subplot(121);
    plot(TTdata.faslen(i,:), UTdata.faslen(i,:),'.'); hold on
    plot([0 10], [0 10],'k-')
    
    subplot(122);
    plot(TTdata.phi(i,:), UTdata.phi(i,:),'.'); hold on
    plot([0 40], [0 40],'k-')
end


end