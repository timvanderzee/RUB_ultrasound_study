clear all; close all; clc
% typical example participant, torques, angles, EMGs
load('cycle_averages_sines_v2.mat')
load('L0.mat')

%% stats
% drift_phi


%%

dcolor = [color(2,:)+[0 .1 .1];.6 .6 .6; color(1,:)];
dcolors = dcolor + .1;
n = 1:118;

ix = [9 2 5];

figure(1)

for j = 1:M

    for m = [2,1,3]
        i = ix(m);
        
        subplot(2,4, 1+ (j-1)*2);
        coord_combine = [[n(:) mean(drift_phi(n,:,j,i),2)+std(drift_phi(n,:,j,i),1,2)]; flipud([n(:) mean(drift_phi(n,:,j,i),2)-std(drift_phi(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
        ylim([0 65])
        box off; hold on
        plot(mean(drift_phi(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(drift_phi(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
        xlim([min(n) max(n)])
        
        subplot(2,4, 5 + (j-1)*2);
        coord_combine = [[n(:) mean(noise_phi(n,:,j,i),2)+std(noise_phi(n,:,j,i),1,2)]; flipud([n(:) mean(noise_phi(n,:,j,i),2)-std(noise_phi(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
        ylim([0 7])
        box off; hold on
        plot(mean(noise_phi(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(noise_phi(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
       
        
        xlim([min(n) max(n)])
        
    end
    
    m = 2;
    i = ix(m);

    subplot(2,4,1 + (j-1)*2);
    plot(mean(drift_phi(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
     plot(max(n), mean(drift_phi(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on

    subplot(2,4, 5+ (j-1)*2);
    plot(mean(noise_phi(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
     plot(max(n), mean(noise_phi(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
end


% X = cell(1,length(Qs));
% 
% for i = 1:length(Qs)
%     X{i} = num2str(Qs(i));
% end
% 
% x = reordercats(categorical(X),X);

%%
% if ishandle(10), close(10); end
% figure(10)
% violinplot(squeeze(noise_phi(length(n),:,j,ix)), [], 'ViolinColor', dcolor); hold on

%%
id = ix;

for j = 1:M
    
%     for i = 1:length(id)
       
        subplot(2,4,2 +(j-1)*2);
        violinplot(squeeze(drift_phi(length(n),:,j,ix)), [], 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on
%         errorbar(i, squeeze(mean(drift_phi(length(n),:,j,id(i)),2)), squeeze(std(drift_phi(length(n),:,j,id(i)),1,2)),'marker','none','linestyle','none','color',dcolor(i,:),'CapSize',2);
%         plot(i, drift_phi(length(n),:,j,id(i)),'ko','markerfacecolor',[1 1 1])
        
        ylim([0 65])
        box off; hold on
%         title(titles{j})

%         xlim([0 length(x)])
        
        subplot(2,4,6 + (j-1)*2);
        violinplot(squeeze(noise_phi(length(n),:,j,ix)), [], 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false);
%         errorbar(i, squeeze(mean(noise_phi(length(n),:,j,id(i)),2)), squeeze(std(noise_phi(length(n),:,j,id(i)),1,2)),'marker','none','linestyle','none','color',dcolor(i,:),'CapSize',2);
%         plot(i, noise_phi(length(n),:,j,id(i)),'ko','markerfacecolor',[1 1 1])
                
        ylim([0 7])
        box off; hold on
%     end
        
%         xlim([0 length(x)])
        
end

%% make nice

figure(1)
set(gcf,'units','normalized','position',[.2 .2 .4 .4])

subplot(241); 
ylabel('Cum. deviation (deg)');
title('All cycles')

subplot(242); 
title('Last cycle')

subplot(243)
title('All cycles')

subplot(244)
title('Last cycle')

subplot(245); 
ylabel('Irregularity (deg)')
xlabel('Cycle #')
xlim([0 20])

subplot(246)
xlabel('Algorithm')

subplot(247); 
xlabel('Cycle #')
xlim([0 20])

subplot(248)
xlabel('Algorithm')

%% statistics
drift = squeeze(drift_phi(length(n),:,j,ix));

[h,p] = ttest(drift(:,1), drift(:,3))


