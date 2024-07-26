clear all; close all; clc
% typical example participant, torques, angles, EMGs
load('cycle_averages_sines_v2.mat')
load('L0.mat')

dcolor = [color(2,:)+[0 .1 .1];.6 .6 .6; color(1,:)];
dcolor = [color(6,:); color(2,:)+[0 .2 .2]; color(4,:)];
dcolors = dcolor;

n = 1:118;
N = 2:119;

ix = [1 9 5];

drift = drift_len;
drift(:,:,2,:) = drift_phi(:,:,1,:);

noise = noise_len;
noise(:,:,2,:) = noise_phi(:,:,1,:);

figure(1)

for j = 1:M

    for m = [2,1,3]
        i = ix(m);
        
        subplot(2,4, 1+ (j-1)*2);
        coord_combine = [[n(:) mean(drift(n,:,j,i),2)+std(drift(n,:,j,i),1,2)]; flipud([n(:) mean(drift(n,:,j,i),2)-std(drift(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
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
        
        subplot(2,4, 5 + (j-1)*2);
        coord_combine = [[n(:) mean(noise(n,:,j,i),2)+std(noise(n,:,j,i),1,2)]; flipud([n(:) mean(noise(n,:,j,i),2)-std(noise(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
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

    subplot(2,4,1 + (j-1)*2);
    plot(mean(drift(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
     plot(max(n), mean(drift(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on

    subplot(2,4, 5+ (j-1)*2);
    plot(mean(noise(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
     plot(max(n), mean(noise(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
end


%%
id = ix;

for j = 1:M
    
%     for i = 1:length(id)
       
        subplot(2,4,2 +(j-1)*2);
        violinplot(squeeze(drift(length(n),:,j,ix)), [], 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false); hold on
%         errorbar(i, squeeze(mean(drift_phi(length(n),:,j,id(i)),2)), squeeze(std(drift_phi(length(n),:,j,id(i)),1,2)),'marker','none','linestyle','none','color',dcolor(i,:),'CapSize',2);
%         plot(i, drift_phi(length(n),:,j,id(i)),'ko','markerfacecolor',[1 1 1])
        
        if j == 1
            ylim([0 180])
        else
            ylim([0 30])
        end
        
        box off; hold on
                set(gca,'YTicklabels',[],'XTicklabels',[])
%         title(titles{j})

%         xlim([0 length(x)])
        
        subplot(2,4,6 + (j-1)*2);
        violinplot(squeeze(mean(noise(:,:,j,ix))), {'UT','TT','UTT'}, 'ViolinColor', dcolor,    'ShowMean', true, 'ShowMedian', false);
        set(gca,'YTicklabels',[])
%         errorbar(i, squeeze(mean(noise_phi(length(n),:,j,id(i)),2)), squeeze(std(noise_phi(length(n),:,j,id(i)),1,2)),'marker','none','linestyle','none','color',dcolor(i,:),'CapSize',2);
%         plot(i, noise_phi(length(n),:,j,id(i)),'ko','markerfacecolor',[1 1 1])
                
        
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
ylabel('Cum. deviation (mm)');

subplot(243)
ylabel('Cum. deviation (deg)');

subplot(245); 
ylabel('Cycle-cycle var. (mm)')
xlabel('Contraction cycle #')

subplot(246)
xlabel('Tracking algorithm')

subplot(247) 
xlabel('Contraction cycle #')
ylabel('Cycle-cycle var. (deg)')

subplot(248)
xlabel('Tracking algorithm')


