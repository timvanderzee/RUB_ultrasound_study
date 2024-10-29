figure(11)
color = get(gca,'colororder');
dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
dcolors = dcolor;

n = 1:118;
N = 1:119;

ix = [2 1 3];
id = ix;

algos = {'UT','TT','UTT'};

%% all trials: effect of image-to-image dissimiliarity
% close all
figure(11)
dcolor = [color(6,:); color(6,:); color(2,:)+[0 .2 .2]; color(2,:)+[0 .2 .2]; color(4,:); color(4,:)];

for k = 1:3
if k == 1
load('ramp_summary.mat','msdphis','msdlens')
len_small = squeeze(msdlens(end,:,1,ix,2));
len_large = squeeze(msdlens(end,:,3,ix,2));
phi_small = squeeze(msdphis(end,:,1,ix,2));
phi_large = squeeze(msdphis(end,:,3,ix,2));

elseif k == 2
load('sines_summary.mat','msdphis','msdlens')
len_small = squeeze(msdlens(10,:,1,ix));
len_large = squeeze(msdlens(10,:,2,ix));
phi_small = squeeze(msdphis(10,:,1,ix));
phi_large = squeeze(msdphis(10,:,2,ix));

elseif k == 3
load('passive_summary.mat','mspen','mslen')
len_small = squeeze(mslen(ix,:,1))';
len_large = squeeze(mslen(ix,:,3))';
phi_small = squeeze(mspen(ix,:,1))';
phi_large = squeeze(mspen(ix,:,3))';
end

len = reshape([len_small;len_large], 8, 6);
phi = reshape([phi_small;phi_large], 8, 6);

xlabels = {'small','large','small','large','small','large'};

d = (k-1)*2+1;
for i = 1:3
    j = (i-1)*2+1;
    
    subplot(3,2,d)
    plot([j j+1], len(:,j:(j+1)),'color',[.8 .8 .8],'linewidth',.5); hold on

    subplot(3,2,d+1)
    plot([j j+1], phi(:,j:(j+1)),'color',[.8 .8 .8],'linewidth',.5); hold on
end

    subplot(3,2,d)
violinplot(len, xlabels,'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

    subplot(3,2,d+1)
violinplot(phi, xlabels,'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on
end

%% make nice
for i = 1:6
subplot(3,2,i)
box off

if mod(i,2)
ylabel('Overall var. (mm)')
else

ylabel('Overall var. (deg)')
end

end
set(gcf,'units','normalized','position',[.2 .2 .4 .5])

subplot(321); title({'Fascicle length', 'Sinusoidal trials'})
subplot(322); title({'Fascicle angle', 'Sinusoidal trials'})

subplot(323); title('Ramp trials')
subplot(324); title('Ramp trials')

subplot(325); title('Passive trials')
subplot(326); title('Passive trials')

%%
% cd('C:\Users\u0167448\OneDrive\8. Ultrasound comparison - TBD\figures\svg\output from MATLAB')
% abc = 'abc';
% 
% plot2svg('FigS2.svg',gcf)
