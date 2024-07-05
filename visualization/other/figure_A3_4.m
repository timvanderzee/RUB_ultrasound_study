clear all; close all; clc

load('cycle_averages_ramps.mat')
load('L0.mat')

%% Figure 4: participant average
if ishandle(10), close(10); end
dcolor = [.5 .5 .5; [229 197 184]/255; color(1,:)];
dcolors = [[.9 .9 .9]*255; 241 224 217; 203 235 255]/255;

ix = [1 9 5];

figure(10)
for j = 1:M
    
    % mean and sd across participants
    MU = mean(mus(:,:,:,j),3);
    SIGMA = std(mus(:,:,:,j),1,3);
    
    for kk = 1:size(MU,1)
    
        x1 = 0:100;

        subplot(N, M, j+M*(kk-1))
        coord_combine = [[x1; MU(kk,:)+SIGMA(kk,:)], fliplr([x1; MU(kk,:)-SIGMA(kk,:)])];
        h = fill(coord_combine(1,:),coord_combine(2,:),'b');
        set(h,'FaceColor',dcolors(end,:),'FaceAlpha',.3,'LineStyle','none');
        plot(x1, MU(kk,:),'linewidth',2,'color',dcolor(end,:)); hold on
        
        ylim([-2 60])
    end

    m = 0;
    for i = iis
        if ismember(i,ix)
            m = m+1;
        
       
            mtsu = mtsus(:,:,j,i);
            mphi = mphis(:,:,j,i) - repmat(reshape(pen0(j,m,:), [1 8]),101,1); % 
            mlen = mlens(:,:,j,i) - repmat(reshape(len0(j,m,:), [1 8]),101,1);  %+  mean(mlens(1,:,j,i),2);

            subplot(N, M, j+M*5)
            coord_combine = [[x2(:) mean(mphi,2)+std(mphi,1,2)]; flipud([x2(:) mean(mphi,2)-std(mphi,1,2)])];
            h = fill(coord_combine(:,1),coord_combine(:,2),'b');
            set(h,'FaceColor',dcolors(m,:),'LineStyle','none');
                ylim([-2 15]); hold on
            plot(x2(:), mean(mphi,2),'linewidth',2,'color',dcolor(m,:)); hold on
            
                subplot(N, M, j+M*6)
            coord_combine = [[x2(:) mean(mlen,2)+std(mlen,1,2)]; flipud([x2(:) mean(mlen,2)-std(mlen,1,2)])];
            h = fill(coord_combine(:,1),coord_combine(:,2),'b');
            set(h,'FaceColor',dcolors(m,:),'LineStyle','none');
                ylim([-30 2])
            xlabel('Contraction cycle (%)')
            hold on; 
            plot(x2(:), mean(mlen,2),'linewidth',2,'color',dcolor(m,:)); hold on

        end
    end
    
end

%% make nice

for i = 1:(N*M)
    subplot(N, M, i);
%         axis tight
    xlim([0 100])
    box off
end


% set(gcf,'units','normalized','position',[.2 0 .2 .9])

set(gcf,'units','normalized','position',[0 0 .4 .99])

subplot(N,M,1); ylabel('Torque (%MVC)')
subplot(N,M,M+1); ylabel('MG (%MVC)')
subplot(N,M,M*2+1); ylabel('LG (%MVC)')
subplot(N,M,M*3+1); ylabel('SOL (%MVC)')
subplot(N,M,M*4+1); ylabel('TA (%MVC)')
subplot(N,M,M*5+1); ylabel('\DeltaPennation (deg)')
subplot(N,M,M*6+1); ylabel('\DeltaLength (mm)')

return
%% Figure 5: variability with respect to the cycle-average
if ishandle(100), close(100); end
dcolor = [color(2,:)+[0 .1 .1];.6 .6 .6; color(1,:)];

dcolors = dcolor + .1;
n = 1:10;

ix = [9 1 5];

figure(100)

for j = 1:M

    for m = [2,1,3]
        i = ix(m);
        
        subplot(M*2, 2, 1 + (j-1)*4);
        coord_combine = [[n(:) mean(msdphis(n,:,j,i),2)+std(msdphis(n,:,j,i),1,2)]; flipud([n(:) mean(msdphis(n,:,j,i),2)-std(msdphis(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
        ylim([0 5])
        box off; hold on
        plot(mean(msdphis(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(msdphis(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
        xlim([min(n) max(n)])
        
        subplot(M*2, 2, 3 + (j-1)*4);
        coord_combine = [[n(:) mean(msdlens(n,:,j,i),2)+std(msdlens(n,:,j,i),1,2)]; flipud([n(:) mean(msdlens(n,:,j,i),2)-std(msdlens(n,:,j,i),1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(m,:).^.1,'LineStyle','none');
        ylim([0 10])
        box off; hold on
        plot(mean(msdlens(:,:,j,i),2),'color',dcolor(m,:),'linewidth',2); hold on
        plot(max(n), mean(msdlens(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on
        
%         xlabel('Cycle #');
        
        xlim([min(n) max(n)])
        
    end
    
    m = 2;
    i = ix(m);

    subplot(M*2, 2, 1 + (j-1)*4);
    plot(mean(msdphis(:,:,j,i),2),'--','color',dcolor(m,:),'linewidth',1); hold on
     plot(max(n), mean(msdphis(length(n),:,j,i),2),'.','color',dcolor(m,:),'markersize',10); hold on

    subplot(M*2, 2, 3 + (j-1)*4);
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
       
        subplot(M*2, 2, 2 + (j-1)*4);
        bar(x,     squeeze(mean(msdphis(length(n),:,j,:),2)),'EdgeColor',color(1,:),'BarWidth',0.7); hold on
        errorbar(1:length(x), squeeze(mean(msdphis(length(n),:,j,:),2)), squeeze(std(msdphis(length(n),:,j,:),1,2)),'marker','none','linestyle','none','color',color(1,:),'CapSize',2);

        ylim([0 5])
        box off; hold on
%         title(titles{j})

%         xlim([0 length(x)])
        
        subplot(M*2, 2, 4 + (j-1)*4);
        bar(x, squeeze(mean(msdlens(length(n),:,j,:),2)),'EdgeColor',color(1,:),'BarWidth',0.7); hold on
        errorbar(1:length(x),   squeeze(mean(msdlens(length(n),:,j,:),2)), squeeze(std(msdlens(length(n),:,j,:),1,2)),'marker','none','linestyle','none','color',color(1,:),'CapSize',2);
        
        ylim([0 10])
        box off; hold on
       
%         xlabel('Q');
        
%         xlim([0 length(x)])
        
end

% subplot(4,2,1); ylabel('\DeltaAngle (deg)');
% title('Temporal deviation')
% 
% subplot(423); ylabel('\DeltaAngle (deg)');
% 
% 
% subplot(422); title('Steady deviation');
% 
% subplot(4,2,5); ylabel('\DeltaLength (mm)');
% subplot(427); ylabel('\DeltaLength (mm)');
set(gcf,'units','normalized','position',[.2 0 .4 1])