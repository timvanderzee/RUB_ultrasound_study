clear all; close all; clc
% typical example participant, torques, angles, EMGs
load('cycle_averages_sines.mat')
load('L0.mat')

%% Figure 4: participant average
if ishandle(10), close(10); end
ix = [1 9 5];

ms = [2 1 3];
M = 2;

dcolor = [[229 197 184]/255; .5 .5 .5; color(1,:)];
dcolors = [241 224 217; [.9 .9 .9]*255; 203 235 255]/255;
% dcolors(dcolors>1) = 1;
x2 = 0:100;

figure(1)


for j = 1:M
    
    % mean and sd across participants
    MU = mean(mus(:,:,:,j),3);
    SIGMA = std(mus(:,:,:,j),1,3);
    
    for kk = 1:size(MU,1)
    
        x1 = 0:100;

        subplot(N, M, j+M*(kk-1))
        coord_combine = [[x1; MU(kk,:)+SIGMA(kk,:)], fliplr([x1; MU(kk,:)-SIGMA(kk,:)])];
        h = fill(coord_combine(1,:),coord_combine(2,:),'b');
        set(h,'FaceColor',dcolors(3,:),'LineStyle','none'); hold on
        plot(x1, MU(kk,:),'linewidth',2,'color',color(1,:)); hold on
        
        ylim([0 40])
    end

    m = 0;
    for i = ix
        m = m+1;

        % N x participant x trial x Q
        mtsu = mtsus(:,:,j,i);
        mphi = mphis(:,:,j,i) - repmat(reshape(pen0(j+4,m,:), [1 8]),101,1); % subtract rest
        mlen = mlens(:,:,j,i) - repmat(reshape(len0(j+4,m,:), [1 8]),101,1); % subtract rest

        subplot(N, M, j+M*5)
        coord_combine = [[x2(:) mean(mphi,2)+std(mphi,1,2)]; flipud([x2(:) mean(mphi,2)-std(mphi,1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(ms(m),:),'LineStyle','none'); hold on
        ylim([-2 15])
        plot(x2(:), mean(mphi,2),'linewidth',2,'color',dcolor(ms(m),:)); hold on
        
        subplot(N, M, j+M*6)
        coord_combine = [[x2(:) mean(mlen,2)+std(mlen,1,2)]; flipud([x2(:) mean(mlen,2)-std(mlen,1,2)])];
        h = fill(coord_combine(:,1),coord_combine(:,2),'b');
        set(h,'FaceColor',dcolors(ms(m),:),'LineStyle','none');
        ylim([-15 15])
        xlabel('Contraction cycle (%)')
        hold on
        plot(x2(:), mean(mlen,2),'linewidth',2,'color',dcolor(ms(m),:)); hold on
    end
    
    i = 9;
    m = 2;

    subplot(N, M, j+M*5)
    mphi = mphis(:,:,j,i) - repmat(reshape(pen0(j+4,m,:), [1 8]),101,1); % subtract rest
    plot(x2(:), mean(mphi,2),'--','linewidth',1,'color',dcolor(ms(m),:)); hold on
    
     subplot(N, M, j+M*6)
        mlen = mlens(:,:,j,i) - repmat(reshape(len0(j+4,m,:), [1 8]),101,1); % subtract rest
    plot(x2(:), mean(mlen,2),'--','linewidth',1,'color',dcolor(ms(m),:)); hold on
    
end

%% make nice

for i = 1:(N*M)
    subplot(N, M, i);
%         axis tight
    xlim([0 100])
    box off
end


% set(gcf,'units','normalized','position',[.2 0 .2 .9])

set(gcf,'units','normalized','position',[0 0 .4 .9])

subplot(N,M,1); ylabel('Torque (N-m)')
subplot(N,M,M+1); ylabel('MG (%MVC)')
subplot(N,M,M*2+1); ylabel('LG (%MVC)')
subplot(N,M,M*3+1); ylabel('SOL (%MVC)')
subplot(N,M,M*4+1); ylabel('TA (%MVC)')
subplot(N,M,M*5+1); ylabel('Pennation (deg)')
subplot(N,M,M*6+1); ylabel('Length (mm)')