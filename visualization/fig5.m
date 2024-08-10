clear all; close all; clc
load('sines_summary.mat', 'drift_len','drift_phi','noise_len','noise_phi')

n = 1:118;
ix = [2 1 3];

drift = drift_len(:,:,2,:);
drift(:,:,2,:) = drift_phi(:,:,2,:);

noise = noise_len(:,:,2,:);
noise(:,:,2,:) = noise_phi(:,:,2,:);

figure(1)
color = get(gca,'colororder');
dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
M = 2;

for j = 1:M

    for m = [2,1,3]
        i = ix(m);
        
        subplot(2,4, 5 + (j-1)*2);
        coord_combine = [[n(:) mean(drift(n,:,j,i),2)+std(drift(n,:,j,i),1,2)]; flipud([n(:) mean(drift(n,:,j,i),2)-std(drift(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolor(m,:).^.1,'LineStyle','none');
        box off; hold on
        plot(mean(drift(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(drift(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
        xlim([min(n) max(n)])
        
        set(gca,'XTicklabels',[])
        
        if j == 1
            ylim([0 180])
        else
            ylim([0 30])
        end
        
        subplot(2,4, 1 + (j-1)*2);
        coord_combine = [[n(:) mean(noise(n,:,j,i),2)+std(noise(n,:,j,i),1,2)]; flipud([n(:) mean(noise(n,:,j,i),2)-std(noise(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolor(m,:).^.1,'LineStyle','none');
        box off; hold on
        plot(mean(noise(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(noise(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on        
        xlim([min(n) max(n)])
        
        if j == 1
            ylim([0 10])
        else
            ylim([0 4])
        end
        
    end
    
    m = 2;
    i = ix(m);

    subplot(2,4,5 + (j-1)*2);
    plot(mean(drift(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
    plot(max(n), mean(drift(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on

    subplot(2,4, 1+ (j-1)*2);
    plot(mean(noise(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
    plot(max(n), mean(noise(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
end


%% violinplots
id = ix;

for j = 1:M
       
    subplot(2,4,6 +(j-1)*2);
    violinplot(squeeze(drift(length(n),:,j,ix)), [], 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on

    if j == 1
        ylim([0 180])
    else
        ylim([0 30])
    end

    box off; hold on
    set(gca,'YTicklabels',[],'XTicklabels',[])

    subplot(2,4,2 + (j-1)*2);
    violinplot(squeeze(mean(noise(:,:,j,ix))), {'UT','TT','UTT'}, 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false);
    set(gca,'YTicklabels',[])

    if j == 1
        ylim([0 10])
    else
        ylim([0 4])
    end

    box off; hold on

        
end

%% make nice
figure(1)
set(gcf,'units','normalized','position',[.2 .2 .4 .4])

subplot(241); 
ylabel('Cycle-cycle var. (mm)');
title('Fascicle length')

subplot(243)
ylabel('Cycle-cycle var. (deg)');
title('Fascicle angle')

subplot(245); 
ylabel('Cum. deviation (mm)')
xlabel('Contraction cycle #')

subplot(246)
xlabel('Tracking algorithm')

subplot(247) 
xlabel('Contraction cycle #')
ylabel('Cum. deviation (deg)')

subplot(248)
xlabel('Tracking algorithm')
