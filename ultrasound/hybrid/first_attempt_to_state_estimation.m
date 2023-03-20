clear all; close all; clc
% cd('C:\Users\timvd\Documents\RUB_ultrasound_study')
% Tmax    = readmatrix('max_torques.txt');
% Trest   = readmatrix('rest_torques.txt'); 

% ramp conditions
force_conditions = {'slow','medium','fast','asym'};
image_qualities = 'high';

cd('C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\ultratrack\individual_subjects\3011')

for j = 1:4
    

load([force_conditions{j},'_high_timtrack.mat'])

TTfaslen = faslen / 564 * 5; % pixels to cm
TTthickness = thickness / 564 * 5; % pixels to cm
TTphi = phi;

% load
load([force_conditions{j},'_',image_qualities,'_01.mat'],'Fdat')

UTfaslen = Fdat.Region.FL / 10; % pixels to cm
UTphi = Fdat.Region.PEN * 180/pi; % pixels to cm

figure(1)
color = get(gca,'colororder');

t = linspace(0,82,2667);

subplot(2,2,j);
plot(t, UTphi); hold on
plot(t, TTphi(1,:))
xlabel('Time (s)'); ylabel('Angle (deg)'); title('Pennation angle'); box off; xlim([0 82])

figure(2)
subplot(2,2,j);
plot(t, UTfaslen); hold on
plot(t, TTfaslen(1,:));
xlabel('Time (s)'); ylabel('Length (cm)'); title('Fascicle length'); box off; xlim([0 82])


%% state-estimation
dUTphi = diff(UTphi);
alpha = .1;
SEphi(1) = TTphi(1,1);

for i = 1:length(TTphi)-1
    
    phi_apriori = SEphi(i) + dUTphi(i);
    
    phi_aposteriori = phi_apriori + alpha * (TTphi(1,i) - phi_apriori);
    
    SEphi(i+1) = phi_aposteriori;
end

figure(1)
subplot(2,2,j)
plot(t, SEphi)

SElen = TTthickness(1,:)./ sind(SEphi);

figure(2)
subplot(2,2,j)
plot(t, SElen)

coef_phi       = polyfit(t,  SEphi,1);
drift_phi(j) = coef_phi(1);

coef_len       = polyfit(t,  SElen,1);
drift_len(j) = coef_len(1);

% noise
noise_phi(j) = mean(abs(diff(SEphi)));
noise_len(j) = mean(abs(diff(SElen)));

end


%% sine conditions
sine_conditions = {'sine_020','sine_1020'};
cd('C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\ultratrack\individual_subjects\3011')

for j = 1:2
    
        load([sine_conditions{j},'_timtrack.mat'])
    
        TTphi = phi;
    TTfaslen = faslen / 564 * 5; % pixels to cm
thickness = thickness / 564 * 5; % pixels to cm
    
    
% load
load([sine_conditions{j},'_01.mat'],'Fdat')

UTfaslen = Fdat.Region.FL / 10; % pixels to cm
UTphi = Fdat.Region.PEN * 180/pi; % pixels to cm

figure(3)
color = get(gca,'colororder');

t = linspace(0,82,2667);

subplot(1,2,j);
plot(t, UTphi); hold on
plot(t, TTphi(1,:))
xlabel('Time (s)'); ylabel('Angle (deg)'); title('Pennation angle'); box off; xlim([0 82])


figure(4)
subplot(1,2,j);
plot(t, UTfaslen); hold on
plot(t, TTfaslen(1,:));
xlabel('Time (s)'); ylabel('Length (cm)'); title('Fascicle length'); box off; xlim([0 82])

%% state-estimation
dUTphi = diff(UTphi);
alpha = .1;
SEphi(1) = TTphi(1,1);

for i = 1:length(TTphi)-1
    
    phi_apriori = SEphi(i) + dUTphi(i);
    
    phi_aposteriori = phi_apriori + alpha * (TTphi(1,i) - phi_apriori);
    
    SEphi(i+1) = phi_aposteriori;
end

figure(3)
subplot(1,2,j)
plot(t, SEphi)

SElen = TTthickness(1,:)./ sind(SEphi);

figure(4)
subplot(1,2,j)
plot(t, SElen)


coef_phi       = polyfit(t,  SEphi,1);
drift_phi(j+4) = coef_phi(1);

coef_len       = polyfit(t,  SElen,1);
drift_len(j+4) = coef_len(1);

% noise
noise_phi(j+4) = mean(abs(diff(SEphi)));
noise_len(j+4) = mean(abs(diff(SElen)));


end

%%

