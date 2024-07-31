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
load('cycle_averages_ramps.mat','msdphis','msdlens')
rmsdphis = msdphis;
rmsdlens = msdlens;

load('cycle_averages_sines_v2.mat','msdphis','msdlens')

% titles = {'Small','Large','Slow','Medium','Fast','Asymmetric'};
close all

quals = {'low','high'};

if ishandle(1), close(1); end; figure(1)
h=tiledlayout(3,24,'TileSpacing','tight','Padding','tight')

js = [1,2,4,5,7,8];

for j = 1:12
        dcolor = [color(6,:); color(6,:); color(2,:)+[0 .2 .2]; color(2,:)+[0 .2 .2]; color(4,:); color(4,:)];
        
     if j < 3
%          dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
         X = msdlens;
         J = j;
         ymax = 60;
         
     elseif j > 2 && j < 5
         X = msdphis;
         J = j-2;
         ymax = 8;
         
     elseif j > 4 && j < 9 
         X = rmsdlens;
         J = j-4;

        ymax = 10;
        
     elseif j > 8
        X = rmsdphis;
        J = j-8;
        ymax = 22;
     end
     
     x = [];
%      disp(titles{j})
    
for m = 1:3
     i = ix(m);
     
     for k = 1:size(X,5)   

         if j < 5
             x = [x X(10,:,J,i,k)' X(end,:,J,i,k)'];

         else
              x = [x X(end,:,J,i,k)'];
        end
            
        disp([algos{m},' - ', quals{k},': ', num2str(round(mean(X(end,:,J,i,k),2),2)), ' - ', num2str(round(std(X(end,:,J,i,k),1,2),2))])
        disp(' ')
    end    
end

%     subplot(3,3, js(j));
    if j < 5
        nexttile([1 6]);
    else
        nexttile([1 3])
    end

    for p = 1:8
        plot(1:size(x,2), x(p,:),'-','color',[.8 .8 .8],'linewidth',.5); hold on
    end
    
    if j > 4
        xlabels = {'low','high','low','high','low','high'};
    else
        xlabels = {'10','final','10','final','10','final'};
    end
    
    violinplot(x, xlabels, 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

    ylim([0 ymax]); 
    box off
    
    if sum(j == [1 5])
        ylabel('Overall var. (mm)');
    elseif sum(j == [3 9])
        ylabel('Overall var. (deg)');
    else
        set(gca,'YTicklabels',[],'YColor',[1 1 1])
    end
    
%     if j < 9
%         set(gca,'XTicklabels',[])
%     end
        
    % connect the dots of one participant

    
    disp(' ')

%     disp(' ')
            

end

%% statistics
% j = 1;
% y = squeeze(rmsdphis(end,:,j,ix,2));
% 
% [~,p1] = ttest(y(:,1), y(:,3)) % UT vs UTT
% [~,p2] = ttest(y(:,2), y(:,3)) % TT vs UTT


%% passive
load('variability_passive.mat','mspen','mslen')
         dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
figure(1)
titles = {'5 deg/s','30 deg/s','120 deg/s'};
% figure(2)

for j = 1:6
    if j < 4
        X = mslen;
        J = j;
        ymax = 10;
    else
        X = mspen;
        J = j-3;
        ymax = 5;
    end
    
    disp(titles{J})
    disp(' ')
    
%     subplot(3,3,j*3);
    nexttile([1 4]);
    categories = repmat({'TT','UT', 'UTT'},8,1);
    
   for p = 1:8
        plot(1:3, X(:,p,J),'-','color',[.8 .8 .8],'linewidth',.5); hold on
    end
    
    
    violinplot(X(:,:,J)', {'UT','TT','UTT'}, 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false);          
    box off
    
    if j == 1
        ylabel('Overall var. (mm)');
    elseif j == 4   
        ylabel('Overall var. (deg)');
    else
        set(gca,'YTicklabels',[],'YColor',[1 1 1])
    end
    
    ylim([0 ymax])
    
   for m = [2,1,3]
        i = ix(m);
            
        disp([algos{m},': ', num2str(round(mean(X(m,:,J),2),2)), ' - ', num2str(round(std(X(m,:,J),1,2),2))])
        disp(' ')
    end           
end

%%
figure(1)
 set(gcf,'units','normalized','position',[.2 .2 .4 .5])
% % 
 cd('C:\Users\u0167448\OneDrive\8. Ultrasound comparison - TBD\figures\svg\output from MATLAB')
% % export_fig test.svg
exportgraphics(gcf,'fig6_out9.png','Resolution',1000)
%% statistics
% j = 1;
% y = squeeze(mspen(:,:,j))';
% 
% [~,p1] = ttest(y(:,1), y(:,3)) % UT vs UTT
% [~,p2] = ttest(y(:,2), y(:,3)) % TT vs UTT

% subplot(241)
% ylabel('Angle deviation (deg)')
%%

% titles = {'Wide: 0-20 %MVC','Narrow: 10-20 %MVC','Slow: 5 deg/s','Slow: 20 %MVC/s','Medium: 33 %MVC/s','Medium: 30 deg/s','Fast: 100 %MVC/s','Slow up - fast down','Fast: 120 deg/s'};
% figure(1)
% for j = 1:9
%     subplot(3,3,j)
%     ylim([0 22])
%     box off
%     title(titles{j})
% end
% 
% subplot(331); ylabel('Overall variability (deg)')
% subplot(334); ylabel('Overall variability (deg)')
% subplot(337); ylabel('Overall variability (deg)')
% 
% set(gcf,'units','normalized','position',[.2 .2 .4 .7])





