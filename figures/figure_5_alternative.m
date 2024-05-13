clear all; close all; clc
% typical example participant, torques, angles, EMGs
load('cycle_averages_sines.mat')
load('L0.mat')

%% stats
% drift_phi


%%

dcolor = [color(2,:)+[0 .1 .1];.6 .6 .6; color(1,:)];
dcolors = dcolor + .1;
n = 1:118;

ix = [9 1 5];

figure(1)

for j = 1:M

    for m = [2,1,3]
        i = ix(m);
        
        subplot(2,4, 1+ (j-1)*2);
        coord_combine = [[n(:) mean(drift_phi(n,:,j,i),2)+std(drift_phi(n,:,j,i),1,2)]; flipud([n(:) mean(drift_phi(n,:,j,i),2)-std(drift_phi(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
        ylim([-2 60])
        box off; hold on
        plot(mean(drift_phi(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(drift_phi(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
        xlim([min(n) max(n)])
        
        subplot(2,4, 5 + (j-1)*2);
        coord_combine = [[n(:) mean(noise_phi(n,:,j,i),2)+std(noise_phi(n,:,j,i),1,2)]; flipud([n(:) mean(noise_phi(n,:,j,i),2)-std(noise_phi(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
        ylim([0 5])
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


X = cell(1,length(Qs));

for i = 1:length(Qs)
    X{i} = num2str(Qs(i));
end

x = reordercats(categorical(X),X);

for j = 1:M
       
        subplot(2,4,2 +(j-1)*2);
        bar(x,     squeeze(mean(drift_phi(length(n),:,j,:),2)),'EdgeColor',color(1,:),'BarWidth',0.7); hold on
        errorbar(1:length(x), squeeze(mean(drift_phi(length(n),:,j,:),2)), squeeze(std(drift_phi(length(n),:,j,:),1,2)),'marker','none','linestyle','none','color',color(1,:),'CapSize',2);

        ylim([-2 60])
        box off; hold on
%         title(titles{j})

%         xlim([0 length(x)])
        
        subplot(2,4,6 + (j-1)*2);
        bar(x, squeeze(mean(noise_phi(length(n),:,j,:),2)),'EdgeColor',color(1,:),'BarWidth',0.7); hold on
        errorbar(1:length(x),   squeeze(mean(noise_phi(length(n),:,j,:),2)), squeeze(std(noise_phi(length(n),:,j,:),1,2)),'marker','none','linestyle','none','color',color(1,:),'CapSize',2);
        
        ylim([0 5])
        box off; hold on
       
        
%         xlim([0 length(x)])
        
end

%% make nice

figure(1)
set(gcf,'units','normalized','position',[.2 .2 .4 .4])

subplot(241); 
ylabel('Drift (deg)');

subplot(242); 

subplot(245); 
ylabel('Noise (deg)')
xlabel('Cycle #')

subplot(247); 
xlabel('Cycle #')

subplot(248)
% xlabel('Q value')