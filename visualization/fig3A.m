% typical example participant, torques, angles, EMGs
% for all trials except passive

Tmax    = readmatrix('max_torques.txt');
Trest   = readmatrix('rest_torques.txt'); 

load('RampTarget.mat','rampTarget')
load('MVC_EMG.mat', 'MVC', 'tnew');

%% time series
force_conditions = {'sine_020','sine_1020'};
p = 1;
titles = {'Large range (0-20% MVC)','Small range (10-20% MVC)'};

figure(1)
color = get(gca,'colororder');
N = 7;
M = length(force_conditions);

c = [.1 .1; .15 .05];

for j = 1:M
    
% load
load([force_conditions{j},'_summary.mat'], 'torque')

% load
load([force_conditions{j},'_EMG.mat'], 'EMG')

% subtract rest torque, divide by MVC torque, multiply by 100%
MGrel = EMG.MG.raw ./ max(MVC.MG,[],2) * 100;
LGrel = EMG.LG.raw ./ max(MVC.LG,[],2) * 100;
TArel = EMG.TA.raw ./ max(MVC.TA,[],2) * 100;
SOrel = EMG.SO.raw ./ max(MVC.SO,[],2) * 100;

MGrel_filt = EMG.MG.filt ./ max(MVC.MG,[],2) * 100;
LGrel_filt = EMG.LG.filt ./ max(MVC.LG,[],2) * 100;
TArel_filt = EMG.TA.filt ./ max(MVC.TA,[],2) * 100;
SOrel_filt = EMG.SO.filt ./ max(MVC.SO,[],2) * 100;

% subtract rest torque, divide by MVC torque, multiply by 100%
Trel = (torque(p,:) - Trest(p)) ./ Tmax(p) * 100;
t = tnew - 1;

Target = c(j,1)*Tmax(p) - c(j,2)*Tmax(p) *cos(2*pi*1.5*(t-1-0.5/1.5));
Target(t < (1/1.5)) = 0;

% subplot(M*2,1,j*2-1);
subplot(N, M, j)
plot(t, Target,  '--','linewidth',2,'color',[.5 .5 .5]); hold on
plot(t, torque(p,:)-Trest(p),'linewidth',2,'color',color(5,:)); hold on
box off; 
ylim([-5 120])

title(titles{j})

subplot(N, M, j+M)
plot(t, MGrel(p,:),'linewidth',1,'color',[.5 .5 .5]); hold on
plot(t, MGrel_filt(p,:),'linewidth',2,'color',color(5,:)); hold on
box off; 
ylim([-100 100])

subplot(N, M, j+M*2)
plot(t, LGrel(p,:),'linewidth',1,'color',[.5 .5 .5]); hold on
plot(t, LGrel_filt(p,:),'linewidth',2,'color',color(5,:)); hold on
box off; 
ylim([-100 100])

subplot(N, M, j+M*3)
plot(t, SOrel(p,:),'linewidth',1,'color',[.5 .5 .5]); hold on
plot(t, SOrel_filt(p,:),'linewidth',2,'color', color(5,:)); hold on
box off; 
ylim([-100 100])

subplot(N, M, j+M*4)
plot(t, TArel(p,:),'linewidth',1,'color',[.5 .5 .5]); hold on
plot(t, TArel_filt(p,:),'linewidth',2,'color',color(5,:)); hold on
box off; 
ylim([-100 100])

end

sgtitle('Sinusoidal trials')

%% ultrasound
filenames = {'*sine_020*','*sine_1020*'}; 
dcolor = [color(2,:)+[0 .2 .2]; color(6,:); color(4,:)];

for k = 1:length(filenames)
    cd([codefolder, '\data\ultrasound\GM\p1'])
    files = dir(filenames{k});

    for i = 1:length(files)
        filename = files(i).name;

        if exist(filename,'file')
            load(filename);

            % recreate t
            t = 0:.03:((2667-1)*.03);

            subplot(N, M, k+M*5)
            plot(t,Fdat.Region.PEN,'color',dcolor(i,:),'linewidth',2); hold on
            ylim([15 30])
            box off

            subplot(N, M, k+M*6)  
            plot(t,Fdat.Region.FL,'color',dcolor(i,:),'linewidth',2); hold on
            ylim([50 75])
            box off
            xlabel('Time (s)')
        end
    end
end

%% make nice
for i = 1:N*M
    subplot(N, M, i);
    xlim([0 8])
end

figure(1)
set(gcf,'units','normalized','position',[0 0 .4 .99])

subplot(N,M,1); ylabel('Torque (N-m)')
subplot(N,M,M+1); ylabel('MG (%MVC)')
subplot(N,M,M*2+1); ylabel('LG (%MVC)')
subplot(N,M,M*3+1); ylabel('SOL (%MVC)')
subplot(N,M,M*4+1); ylabel('TA (%MVC)')
subplot(N,M,M*5+1); ylabel('Pennation (deg)')
subplot(N,M,M*6+1); ylabel('Length (mm)')

