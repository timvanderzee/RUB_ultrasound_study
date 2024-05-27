clear all; close all; clc
% typical example participant, torques, angles, EMGs


figure(1)
color = get(gca,'colororder');
dcolor = [color(2,:)+[0 .1 .1];.6 .6 .6; color(1,:)];
dcolors = dcolor + .1;
n = 1:118;
N = 1:119;

ix = [9 1 5];
id = ix;

%% passive
load('variability_passive.mat','pen_sd')


if ishandle(1), close(1); end; figure(1)
titles = {'5 deg/s','30 deg/s','120 deg/s'};

for j = 1:3

    subplot(1,3,j);
    categories = repmat({'TT','UT', 'UTT'},8,1);
    violinplot(pen_sd(:,:,j)', categories, 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false);          

                
end

%%
figure(1)
for j = 1:3
    subplot(1,3,j)
        ylim([0 5])
    box off; hold on
    title(titles{j})
end

subplot(1,3,1)
ylabel('Angle deviation (deg)')
set(gcf,'units','normalized','position',[.2 .5 .4 .2])

%% ramp
load('cycle_averages_ramps.mat','msdphis')


figure(2)
titles = {'Slow','Medium','Fast','Asymmetric'};

for j = 1:4

    subplot(1,4,j);
    violinplot(squeeze(msdphis(end,:,j,ix)), categories, 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false);          

                
end

subplot(141)
ylabel('Angle deviation (deg)')
%%
figure(2)
for j = 1:4
    subplot(1,4,j)
        ylim([0 10])
    box off; hold on
    title(titles{j})
end


set(gcf,'units','normalized','position',[.2 .2 .4 .2])



