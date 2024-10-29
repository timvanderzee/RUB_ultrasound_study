dt = 1/30;

Qs = [(.0001:.0001:.001) .002 .003 .004 .005 .01 .1 1 10 100 1000];

N = 601;
t = 0:dt:((N-1)*dt);



%% manual tracking
obs = {'TZ', 'PT', 'BR'};
ls = {'-','--',':'};

Ns = [101; 84];
vids = {'TA_video','calf_raise'};
vidnames = {'exampleVideoTA_Verheul_2022','calf_raise'};

pixpermms =  [13.6 508 / 65];

for j = 1:2
pixpermm = pixpermms(j);
    Fangle = nan(Ns(j),3);
    Flength = nan(Ns(j),3);
    Aangle = nan(Ns(j),3);

    for i = 1:length(obs)
        load([codefolder, '\data\ultrasound\',vids{j},'\manual\Manual_Tracking_',obs{i},'.mat'])

        [ts, id] = sort(t(FasData.digitizedFrames));

        FasData.FAngle(FasData.FAngle<0) = FasData.FAngle(FasData.FAngle<0) + 180;

        Fangle(1:length(id),i) = FasData.FAngle(id)';
        Aangle(1:length(id),i) = ApoData.Angle(id)';
        Flength(1:length(id),i) = FasData.FLength(id)'/pixpermm;
    end

    Pangle = Fangle - Aangle;
    clear RMSE

    for i = 1:length(Qs)
        if Qs(i) < 1
            f = ['%.',num2str(-floor(log10(Qs(i)))),'f'];
        else
            f = '%.0f';
        end
        filename = [vidnames{j},'_tracked_Q=',strrep(num2str(Qs(i), f),'.',''), '.mat'];

        load(filename)
        % load('exampleVideoTA_Verheul_2022.mat')


        ANG = Fdat.Region.ANG;
        PEN = Fdat.Region.PEN;
        FL = Fdat.Region.FL;


        RMSE(i,:) = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))]
    end

if j == 1
    %% Hybrid Track
    load('exampleVideoTA_trackingResults.mat')

    N = 601;
    clear FL PEN
    for i = 1:601
        FL(i) = Hybrid(i).length;
        PEN(i) = Hybrid(i).angle;
    end

    t = 0:dt:((N-1)*dt);

    HRMSE = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))];
else
    Ldata = readmatrix('calf_raise.xlsx', 'Sheet','Fasc_length');
    Pdata = readmatrix('calf_raise.xlsx', 'Sheet','Pennation');
    Ldata = Ldata(2:end,2:end);
    Pdata = Pdata(2:end,2:end);

    N = size(Ldata,1);
    PEN_DLT = nan(1,size(Ldata,1));
    FL_DLT = nan(1,size(Ldata,1));

    for i = 1:size(Ldata,1)
        FL(i) = median(Ldata(i,:),'omitnan') / pixpermm;
        PEN(i) = median(Pdata(i,:),'omitnan');
    end

    t = 0:dt:((N-1)*dt);

    HRMSE = [sqrt(mean((FL(FasData.digitizedFrames)' - mean(Flength,2)).^2,'omitnan')) sqrt(mean((PEN(FasData.digitizedFrames)' - mean(Pangle,2)).^2,'omitnan'))];
end
    %%
    figure(10)
    subplot(2,2,(j-1)*2+1);
    semilogx(Qs, RMSE(:,1),'linewidth',2); hold on
    plot([min(Qs) max(Qs)], [HRMSE(1) HRMSE(1)],'k-')   
    ylabel('RMSD (mm)')

    if j == 1
        ylim([0 6])
    else
        ylim([0 5])
    end

    title('Fascicle length')

    subplot(2,2,(j-1)*2+2);
    semilogx(Qs, RMSE(:,2),'linewidth',2); hold on
    plot([min(Qs) max(Qs)], [HRMSE(2) HRMSE(2)],'k-')
    box off

    if j == 1
        ylim([0 1])
    else
        ylim([0 5])
    end

    ylabel('RMSD (deg)')
    title('Fascicle angle')

end

%%
xt = 10.^(-4:3);
for i = 1:4
    subplot(2,2,i)
    xlabel('Process noise covariance parameter c')
    box off
    xlim([min(Qs) max(Qs)])

    set(gca,'xtick',xt)
end

set(gcf,'units','normalized','position',[.1 .1 .35 .6])

