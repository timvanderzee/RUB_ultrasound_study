figure(5)
color = get(gca,'colororder');
dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
dcolors = dcolor;

n = 1:118;
N = 1:119;

ix = [2 1 3];
id = ix;

algos = {'UT','TT','UTT'};

%% sines: effect of contraction number
load('sines_summary.mat','msdphis','msdlens')
dcolor = [color(6,:); color(6,:); color(2,:)+[0 .2 .2]; color(2,:)+[0 .2 .2]; color(4,:); color(4,:)];
        
phi_10 = squeeze(reshape(msdphis(10,:,1:2,ix), [1 16 3]));
phi_last = squeeze(reshape(msdphis(end,:,1:2,ix), [1 16 3]));
phi = reshape([phi_10; phi_last], 16, 6);

len_10 = squeeze(reshape(msdlens(10,:,1:2,ix), [1 16 3]));
len_last = squeeze(reshape(msdlens(end,:,1:2,ix), [1 16 3]));
len = reshape([len_10; len_last], 16, 6);

xlabels = {'short','long','short','long','short','long'};

figure(5)

for i = 1:3
    j = (i-1)*2+1;
    subplot(321)
    plot([j j+1], len(:,j:(j+1)),'color',[.8 .8 .8],'linewidth',.5); hold on

    subplot(232)
    plot([j j+1], phi(:,j:(j+1)),'color',[.8 .8 .8],'linewidth',.5); hold on
end


subplot(121);
violinplot(len, xlabels,'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

subplot(122);
violinplot(phi, xlabels,'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

set(gcf,'units','normalized','position',[.2 .7 .4 .2])

% x = squeeze(reshape(msdphis(end,:,1:2,ix), [1 16 3]));
% 
% subplot(222);
% violinplot(x, xlabels, 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

%% ramps: effect of image quality
figure(6)
load('ramp_summary.mat','msdphis','msdlens')

xlabels = {'high','low','high','low','high','low'};

len_low = squeeze(reshape(msdlens(end,:,1:4,ix,1), [1 32 3]));
len_high = squeeze(reshape(msdlens(end,:,1:4,ix,2), [1 32 3]));
len = reshape([len_high; len_low], 32, 6);

phi_low = squeeze(reshape(msdphis(end,:,1:4,ix,1), [1 32 3]));
phi_high = squeeze(reshape(msdphis(end,:,1:4,ix,2), [1 32 3]));
phi = reshape([phi_high; phi_low], 32, 6);

for i = 1:3
    j = (i-1)*2+1;
    subplot(121)
    plot([j j+1], len(:,j:(j+1)),'color',[.8 .8 .8],'linewidth',.5); hold on

    subplot(122)
    plot([j j+1], phi(:,j:(j+1)),'color',[.8 .8 .8],'linewidth',.5); hold on
end

subplot(121);
violinplot(len, xlabels,'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

subplot(122);
violinplot(phi, xlabels,'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

set(gcf,'units','normalized','position',[.2 .4 .4 .2])
%% all trials: effect of image-to-image dissimiliarity
figure(7)
% close all

load('ramp_summary.mat','msdphis','msdlens')
len_small_ramp = squeeze(msdlens(end,:,1,ix,2));
len_large_ramp = squeeze(msdlens(end,:,3,ix,2));
phi_small_ramp = squeeze(msdphis(end,:,1,ix,2));
phi_large_ramp = squeeze(msdphis(end,:,3,ix,2));

load('sines_summary.mat','msdphis','msdlens')
len_small_sine = squeeze(msdlens(10,:,1,ix));
len_large_sine = squeeze(msdlens(10,:,2,ix));
phi_small_sine = squeeze(msdphis(10,:,1,ix));
phi_large_sine = squeeze(msdphis(10,:,2,ix));

load('passive_summary.mat','mspen','mslen')
len_small_pas = squeeze(mslen(ix,:,1))';
len_large_pas = squeeze(mslen(ix,:,3))';
phi_small_pas = squeeze(mspen(ix,:,1))';
phi_large_pas = squeeze(mspen(ix,:,3))';

% len_small_sine = [];
% len_large_sine = [];
% phi_small_sine = [];
% phi_large_sine = [];
% 
% len_small_ramp = [];
% len_large_ramp = [];
% phi_small_ramp = [];
% phi_large_ramp = [];
% 
% len_small_pas = [];
% len_large_pas = [];
% phi_small_pas = [];
% phi_large_pas = [];

len = reshape([len_small_sine;len_small_ramp; len_small_pas;len_large_sine;len_large_ramp;len_large_pas], 24, 6);
phi = reshape([phi_small_sine;phi_small_ramp; phi_small_pas;phi_large_sine;phi_large_ramp;phi_large_pas], 24, 6);

xlabels = {'small','large','small','large','small','large'};

for i = 1:3
    j = (i-1)*2+1;
    subplot(121)
    plot([j j+1], len(:,j:(j+1)),'color',[.8 .8 .8],'linewidth',.5); hold on

    subplot(122)
    plot([j j+1], phi(:,j:(j+1)),'color',[.8 .8 .8],'linewidth',.5); hold on
end

subplot(121);
violinplot(len, xlabels,'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

subplot(122);
violinplot(phi, xlabels,'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

set(gcf,'units','normalized','position',[.2 .1 .4 .2])

%% make nice
for j = 5:7
    figure(j)
    for i = 1:2
    subplot(1,2,i)
    box off

    if mod(i,2)
    ylabel('Overall var. (mm)')
    else

    ylabel('Overall var. (deg)')
    end

    end


    if j == 5
        subplot(121); title({'Fascicle length', 'Sequence duration'})
        subplot(122); title({'Fascicle angle', 'Sequence duration'})

    elseif j == 6
    subplot(121); title('Image quality')
    subplot(122); title('Image quality')

    elseif j == 7
    subplot(121); title('Image-to-image dissimiarity')
    subplot(122); title('Image-to-image dissimiarity')
    end
end

%% save
cd('C:\Users\u0167448\OneDrive\8. Ultrasound comparison - TBD\figures\svg\output from MATLAB')
abc = 'abc';

for j = 5:7
    figure(j)
    plot2svg(['Fig6', abc(j-4),'.svg'],gcf)
end


