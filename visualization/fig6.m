clear all; close all; clc

figure(1)
color = get(gca,'colororder');
dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
dcolors = dcolor;

n = 1:118;
N = 1:119;

ix = [1 9 5];
ix = [2 1 3];
id = ix;

algos = {'UT','TT','UTT'};

%% ramp
load('ramp_summary.mat','msdphis','msdlens')
rmsdphis = msdphis;
rmsdlens = msdlens;

load('sines_summary.mat','msdphis','msdlens')

% titles = {'Small','Large','Slow','Medium','Fast','Asymmetric'};
close all

quals = {'low','high'};

if ishandle(1), close(1); end; figure(1)
h=tiledlayout(3,24,'TileSpacing','tight','Padding','tight')

js = [1,2,4,5,7,8];
titles = {'Sine: small range', 'Sine: large range', 'Sine: small range', 'Sine: large range', 'Ramp: slow', 'Ramp: moderate', 'Ramp: fast', 'Ramp: asymmetric', 'Ramp: slow', 'Ramp: moderate', 'Ramp: fast', 'Ramp: asymmetric'};
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
        title(titles{j})
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

    disp(' ')   

end




%% passive
load('passive_summary.mat','mspen','mslen')
         dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
figure(1)
titles = {'Passive: slow', 'Passive: moderate', 'Passive: fast','Passive: slow', 'Passive: moderate', 'Passive: fast'};
ip = [2, 1, 3];
for j = 1:6
    if j < 4
        X = mslen(ip,:,:);
        J = j;
        ymax = 10;
    else
        X = mspen(ip,:,:);
        J = j-3;
        ymax = 5;
    end
    
    disp(titles{J})
    disp(' ')
    
    nexttile([1 4]);
    categories = repmat({'TT','UT', 'UTT'},8,1);
    
   for p = 1:8
        plot(1:3, X(:,p,J),'-','color',[.8 .8 .8],'linewidth',.5); hold on
    end
    
    title(titles{j})
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

