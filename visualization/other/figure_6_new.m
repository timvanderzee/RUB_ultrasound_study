clear all; close all; clc
% typical example participant, torques, angles, EMGs

figure(1)
color = get(gca,'colororder');
% dcolor = [color(2,:)+[0 .1 .1];.6 .6 .6; color(1,:)];
% dcolors = dcolor + .1;
dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
dcolors = dcolor;

n = 1:118;
N = 1:119;

ix = [1 9 5];
id = ix;

algos = {'UT','TT','UTT'};


%% ramp
load('cycle_averages_ramps.mat','msdphis')
rmsdphis = msdphis;

load('cycle_averages_sines_v2.mat','msdphis')

titles = {'Small','Large','Slow','Medium','Fast','Asymmetric'};
close all

quals = {'low','high'};

if ishandle(1), close(1); end; figure(1)

js = [1,2,4,5,7,8];

for j = 1:6
     if j < 3
         X = msdphis;
         dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
         J = j;
     else
         X = rmsdphis;
         J = j-2;
        dcolor = [color(6,:); color(6,:); color(2,:)+[0 .2 .2]; color(2,:)+[0 .2 .2]; color(4,:); color(4,:)];
     end
     
     x = [];
     disp(titles{j})
    
for m = 1:3
     i = ix(m);
     
     for k = 1:size(X,5)   

        x = [x X(end,:,J,i,k)'];
            
        disp([algos{m},' - ', quals{k},': ', num2str(round(mean(X(end,:,J,i,k),2),2)), ' - ', num2str(round(std(X(end,:,J,i,k),1,2),2))])
        disp(' ')
    end    
end

    subplot(3,3, js(j));
    for p = 1:8
        plot(1:size(x,2), x(p,:),'-','color',[.8 .8 .8],'linewidth',.5); hold on
    end
    
    violinplot(x, [], 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

    % connect the dots of one participant

    
    disp(' ')

%     disp(' ')
            

end

%% statistics
j = 1;
y = squeeze(rmsdphis(end,:,j,ix,2));

[~,p1] = ttest(y(:,1), y(:,3)) % UT vs UTT
[~,p2] = ttest(y(:,2), y(:,3)) % TT vs UTT


%% passive
load('variability_passive.mat','mspen')
         dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
figure(1)
titles = {'5 deg/s','30 deg/s','120 deg/s'};
% figure(2)

for j = 1:3
    disp(titles{j})
    disp(' ')
    
    subplot(3,3,j*3);
    categories = repmat({'TT','UT', 'UTT'},8,1);
    
   for p = 1:8
        plot(1:3, mspen(:,p,j),'-','color',[.8 .8 .8],'linewidth',.5); hold on
    end
    
    
    violinplot(mspen(:,:,j)', [], 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false);          

   for m = [2,1,3]
        i = ix(m);
            
        disp([algos{m},': ', num2str(round(mean(mspen(m,:,j),2),2)), ' - ', num2str(round(std(mspen(m,:,j),1,2),2))])
        disp(' ')
    end           
end

%%
%% statistics
j = 1;
y = squeeze(mspen(:,:,j))';

[~,p1] = ttest(y(:,1), y(:,3)) % UT vs UTT
[~,p2] = ttest(y(:,2), y(:,3)) % TT vs UTT

% subplot(241)
% ylabel('Angle deviation (deg)')
%%

titles = {'Wide: 0-20 %MVC','Narrow: 10-20 %MVC','Slow: 5 deg/s','Slow: 20 %MVC/s','Medium: 33 %MVC/s','Medium: 30 deg/s','Fast: 100 %MVC/s','Slow up - fast down','Fast: 120 deg/s'};
figure(1)
for j = 1:9
    subplot(3,3,j)
    ylim([0 22])
    box off
    title(titles{j})
end

subplot(331); ylabel('Overall variability (deg)')
subplot(334); ylabel('Overall variability (deg)')
subplot(337); ylabel('Overall variability (deg)')

set(gcf,'units','normalized','position',[.2 .2 .4 .7])





