clear all; close all; clc
% average over contractions


% ramp conditions
conditions = {'slow_high', 'medium_high', 'fast_high', 'asym_high','sine_020', 'sine_1020'};

tpoints = 0:8:80;
figure(1); 
color = get(gca,'colororder');
close
    
for p = 1:5
    
for i = 2
for j = 1:6
    
    
TTdata = load([conditions{j},'_timtrack.mat']);

TTdata.faslen = TTdata.faslen / 564 * 5; % pixels to cm
TTdata.thickness = TTdata.thickness / 564 * 5; % pixels to cm

UTdata = load([conditions{j},'_ultratrack.mat']);


clear xi

for k = 1:length(tpoints)-1
    x1 = TTdata.faslen(p,TTdata.t >= tpoints(k) & TTdata.t <= tpoints(k+1));
    x2 = TTdata.phi(p,TTdata.t >= tpoints(k) & TTdata.t <= tpoints(k+1));
    Li(k,:) = interp1(linspace(0,1, length(x1)), x1, linspace(0,1,267));
    Pi(k,:) = interp1(linspace(0,1, length(x2)), x2, linspace(0,1,267));
end


uLi = UTdata.faslen(p,:) - UTdata.faslen(p,2) + mean(Li(:,1));
uPi = UTdata.phi(p,:) - UTdata.phi(p,2) + mean(Pi(:,1));

ti = linspace(0,8,267);
figure(p*10)
set(gcf,'name',['Subject: ', num2str(p), ' Fascicle length'])
subplot(3,2,j)
% plot(xi','color',color(j,:)); hold on
plot(ti,mean(Li),'-','linewidth',2,'color',color(1,:)); hold on
plot(ti,mean(Li)+std(Li),'-','linewidth',1,'color',color(1,:))
plot(ti,mean(Li)-std(Li),'-','linewidth',1,'color',color(1,:))

plot(TTdata.t, uLi,'linewidth',2,'color',color(2,:)); 
xlim([0 8]); ylim([2 8]); title(conditions{j}(1:4));
xlabel('Time (s)'); ylabel('Angle (deg)'); box off

figure(p*10+1)
set(gcf,'name',['Subject: ', num2str(p), ' Pennation'])
subplot(3,2,j)
% plot(xi','color',color(j,:)); hold on
plot(ti,mean(Pi),'-','linewidth',2,'color',color(1,:)); hold on
plot(ti,mean(Pi)+std(Pi),'-','linewidth',1,'color',color(1,:))
plot(ti,mean(Pi)-std(Pi),'-','linewidth',1,'color',color(1,:))

plot(TTdata.t, uPi,'linewidth',2,'color',color(2,:))
xlim([0 8]); ylim([10 50]); title(conditions{j}(1:4));
xlabel('Time (s)'); ylabel('Angle (deg)'); box off

% save contraction-average
mLi(j,:,p) = mean(Li);
mPi(j,:,p) = mean(Pi);

muLi(j,:,p) = interp1(TTdata.t, uLi, ti);
muPi(j,:,p) = interp1(TTdata.t, uPi, ti);

end
end
end

%%
if ishandle(100), close(100); end 
if ishandle(101), close(101); end 

for j = 1:6

figure(100); subplot(3,2,j); 
plot(ti, mean(mLi(j,:,:),3)); hold on;
plot(ti, mean(muLi(j,:,:),3))
ylim([3.5 7])
xlabel('Time (s)'); ylabel('Length (cm)'); box off

figure(101); subplot(3,2,j); 
plot(ti, mean(mPi(j,:,:),3)); hold on
plot(ti, mean(muPi(j,:,:),3))
ylim([15 30]); 
xlabel('Time (s)'); ylabel('Angle (deg)'); box off
end