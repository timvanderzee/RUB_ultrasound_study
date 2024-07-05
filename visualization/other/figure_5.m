clear all; close all; clc
% typical example participant, torques, angles, EMGs
load('cycle_averages_sines.mat')
load('L0.mat')

%% Figure 5: variability with respect to the cycle-average
if ishandle(2), close(2); end

dcolor = [color(2,:)+[0 .1 .1];.6 .6 .6; color(1,:)];
dcolors = dcolor + .1;
n = 1:119;

ix = [9 1 5];

figure(2)

for j = 1:M

    for m = [2,1,3]
        i = ix(m);
        
        subplot(4,2, 1+ (j-1)*4);
        coord_combine = [[n(:) mean(msdphis(n,:,j,i),2)+std(msdphis(n,:,j,i),1,2)]; flipud([n(:) mean(msdphis(n,:,j,i),2)-std(msdphis(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
        ylim([0 5])
        box off; hold on
        plot(mean(msdphis(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(msdphis(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
        xlim([min(n) max(n)])
        
        subplot(4,2, 3 + (j-1)*4);
        coord_combine = [[n(:) mean(msdlens(n,:,j,i),2)+std(msdlens(n,:,j,i),1,2)]; flipud([n(:) mean(msdlens(n,:,j,i),2)-std(msdlens(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
        ylim([0 10])
        box off; hold on
        plot(mean(msdlens(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(msdlens(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
       
        
        xlim([min(n) max(n)])
        
    end
    
    m = 2;
    i = ix(m);

    subplot(4,2,1 + (j-1)*4);
    plot(mean(msdphis(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
     plot(max(n), mean(msdphis(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on

    subplot(4,2, 3+ (j-1)*4);
    plot(mean(msdlens(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
     plot(max(n), mean(msdlens(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
end


% set(gcf,'units','normalized','position',[.2 .2 .5 .5])

% final variability wrt Q value
% if ishandle(200), close(200); end
% % figure(200)

X = cell(1,length(Qs));

for i = 1:length(Qs)
    X{i} = num2str(Qs(i));
end

x = reordercats(categorical(X),X);

for j = 1:M
       
        subplot(4,2,2 +(j-1)*4);
        bar(x,     squeeze(mean(msdphis(length(n),:,j,:),2)),'EdgeColor',color(1,:),'BarWidth',0.7); hold on
        errorbar(1:length(x), squeeze(mean(msdphis(length(n),:,j,:),2)), squeeze(std(msdphis(length(n),:,j,:),1,2)),'marker','none','linestyle','none','color',color(1,:),'CapSize',2);

        ylim([0 5])
        box off; hold on
%         title(titles{j})

%         xlim([0 length(x)])
        
        subplot(4,2,4 + (j-1)*4);
        bar(x, squeeze(mean(msdlens(length(n),:,j,:),2)),'EdgeColor',color(1,:),'BarWidth',0.7); hold on
        errorbar(1:length(x),   squeeze(mean(msdlens(length(n),:,j,:),2)), squeeze(std(msdlens(length(n),:,j,:),1,2)),'marker','none','linestyle','none','color',color(1,:),'CapSize',2);
        
        ylim([0 10])
        box off; hold on
       
        
%         xlim([0 length(x)])
        
end

subplot(4,2,1); ylabel('\DeltaAngle (deg)');
title('Temporal deviation')

subplot(423); ylabel('\DeltaAngle (deg)');
subplot(422); title('Steady deviation');

subplot(4,2,5); ylabel('\DeltaLength (mm)');
subplot(427); ylabel('\DeltaLength (mm)');
set(gcf,'units','normalized','position',[.2 .1 .4 .8])
xlabel('Cycle #')
subplot(428)
xlabel('Q value')