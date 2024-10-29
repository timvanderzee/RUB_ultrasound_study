clear all; close all; clc

load('passive_summary.mat','mspen','mslen')
ip = [2, 1, 3];

load('ramp_summary.mat','msdphis','msdlens')
ix = [2, 3, 1];

rmsdphis = msdphis;
rmsdlens = msdlens;

load('sines_summary.mat','msdphis','msdlens','noise_phi','drift_phi','noise_len','drift_len')

%% noise and drift
% cycles x participants x range x algos
algos = {'UT','UTT','TT'};

for i = 1:3
    disp([algos{i}, ' drift: ', num2str(round([mean(drift_len(end,:,1,ix(i))) std(drift_len(end,:,1,ix(i)))],1)), ' mm'])
    disp([algos{i}, ' noise: ', num2str(round([mean(mean(noise_len(:,:,1,ix(i)))) std(mean(noise_len(:,:,1,ix(i))))],1)), ' mm'])    

    disp([algos{i}, ' drift: ', num2str(round([mean(drift_phi(end,:,1,ix(i))) std(drift_phi(end,:,1,ix(i)))],1)), ' deg'])
    disp([algos{i}, ' noise: ', num2str(round([mean(mean(noise_phi(:,:,1,ix(i)))) std(mean(noise_phi(:,:,1,ix(i))))],1)), ' deg'])
    
end

id = 2;

% UTT vs UT
[~,p1] = ttest(drift_len(end,:,id,ix(1)), drift_len(end,:,id,ix(2)));
[~,p2] = ttest(mean(noise_len(:,:,id,ix(1))), mean(noise_len(:,:,id,ix(2))));
[~,p3] = ttest(drift_phi(end,:,id,ix(1)), drift_phi(end,:,id,ix(2)));
[~,p4] = ttest(mean(noise_phi(:,:,id,ix(1))), mean(noise_phi(:,:,id,ix(2))));

% UTT vs TT
[~,p5] = ttest(drift_len(end,:,id,ix(3)), drift_len(end,:,id,ix(2)));
[~,p6] = ttest(mean(noise_len(:,:,id,ix(3))), mean(noise_len(:,:,id,ix(2))));
[~,p7] = ttest(drift_phi(end,:,id,ix(3)), drift_phi(end,:,id,ix(2)));
[~,p8] = ttest(mean(noise_phi(:,:,id,ix(3))), mean(noise_phi(:,:,id,ix(2))));

disp('UTT vs UT')
disp('Drift')
disp(['Length: p = ', num2str(p1), ' Angle: p = ', num2str(p3)])
disp('Noise')
disp(['Length: p = ', num2str(p2), ' Angle: p = ', num2str(p4)])

disp(' ')
disp('UTT vs TT')
disp('Drift')
disp(['Length: p = ', num2str(p5), ' Angle: p = ', num2str(p7)])
disp('Noise')
disp(['Length: p = ', num2str(p6), ' Angle: p = ', num2str(p8)])


%% variability: pennation

PEN.sine.low.UT = squeeze(msdphis([10, 119], :, 1, ix(1)))';
PEN.sine.low.TT = squeeze(msdphis([10, 119], :, 1, ix(3)))';
PEN.sine.low.UTT = squeeze(msdphis([10, 119], :, 1, ix(2)))';

PEN.sine.high.UT = squeeze(msdphis([10, 119], :, 2, ix(1)))';
PEN.sine.high.TT = squeeze(msdphis([10, 119], :, 2, ix(3)))';
PEN.sine.high.UTT = squeeze(msdphis([10, 119], :, 2, ix(2)))';

PEN.ramp.slow.UT = [squeeze(rmsdphis(10, :, 1, ix(1),1))' squeeze(rmsdphis(10, :, 1, ix(1),2))'];
PEN.ramp.slow.TT = [squeeze(rmsdphis(10, :, 1, ix(3),1))' squeeze(rmsdphis(10, :, 1, ix(3),2))'];
PEN.ramp.slow.UTT = [squeeze(rmsdphis(10, :, 1, ix(2),1))' squeeze(rmsdphis(10, :, 1, ix(2),2))'];

PEN.ramp.med.UT = [squeeze(rmsdphis(10, :, 2, ix(1),1))' squeeze(rmsdphis(10, :, 2, ix(1),2))'];
PEN.ramp.med.TT = [squeeze(rmsdphis(10, :, 2, ix(3),1))' squeeze(rmsdphis(10, :, 2, ix(3),2))'];
PEN.ramp.med.UTT = [squeeze(rmsdphis(10, :, 2, ix(2),1))' squeeze(rmsdphis(10, :, 2, ix(2),2))'];

PEN.ramp.fast.UT = [squeeze(rmsdphis(10, :, 3, ix(1),1))' squeeze(rmsdphis(10, :, 3, ix(1),2))'];
PEN.ramp.fast.TT = [squeeze(rmsdphis(10, :, 3, ix(3),1))' squeeze(rmsdphis(10, :, 3, ix(3),2))'];
PEN.ramp.fast.UTT = [squeeze(rmsdphis(10, :, 3, ix(2),1))' squeeze(rmsdphis(10, :, 3, ix(2),2))'];

PEN.ramp.asym.UT = [squeeze(rmsdphis(10, :, 4, ix(1),1))' squeeze(rmsdphis(10, :, 4, ix(1),2))'];
PEN.ramp.asym.TT = [squeeze(rmsdphis(10, :, 4, ix(3),1))' squeeze(rmsdphis(10, :, 4, ix(3),2))'];
PEN.ramp.asym.UTT = [squeeze(rmsdphis(10, :, 4, ix(2),1))' squeeze(rmsdphis(10, :, 4, ix(2),2))'];

PEN.pas.slow.UT = squeeze(mspen(ip(1),:,1))';
PEN.pas.slow.TT = squeeze(mspen(ip(2),:,1))';
PEN.pas.slow.UTT = squeeze(mspen(ip(3),:,1))';

PEN.pas.med.UT = squeeze(mspen(ip(1),:,2))';
PEN.pas.med.TT = squeeze(mspen(ip(2),:,2))';
PEN.pas.med.UTT = squeeze(mspen(ip(3),:,2))';

PEN.pas.fast.UT = squeeze(mspen(ip(1),:,3))';
PEN.pas.fast.TT = squeeze(mspen(ip(2),:,3))';
PEN.pas.fast.UTT = squeeze(mspen(ip(3),:,3))';

%% length
LEN.sine.low.UT = squeeze(msdlens([10, 119], :, 1, ix(1)))';
LEN.sine.low.TT = squeeze(msdlens([10, 119], :, 1, ix(3)))';
LEN.sine.low.UTT = squeeze(msdlens([10, 119], :, 1, ix(2)))';

LEN.sine.high.UT = squeeze(msdlens([10, 119], :, 2, ix(1)))';
LEN.sine.high.TT = squeeze(msdlens([10, 119], :, 2, ix(3)))';
LEN.sine.high.UTT = squeeze(msdlens([10, 119], :, 2, ix(2)))';


LEN.ramp.slow.UT = [squeeze(rmsdlens(10, :, 1, ix(1),1))' squeeze(rmsdlens(10, :, 1, ix(1),2))'];
LEN.ramp.slow.TT = [squeeze(rmsdlens(10, :, 1, ix(3),1))' squeeze(rmsdlens(10, :, 1, ix(3),2))'];
LEN.ramp.slow.UTT = [squeeze(rmsdlens(10, :, 1, ix(2),1))' squeeze(rmsdlens(10, :, 1, ix(2),2))'];

LEN.ramp.med.UT = [squeeze(rmsdlens(10, :, 2, ix(1),1))' squeeze(rmsdlens(10, :, 2, ix(1),2))'];
LEN.ramp.med.TT = [squeeze(rmsdlens(10, :, 2, ix(3),1))' squeeze(rmsdlens(10, :, 2, ix(3),2))'];
LEN.ramp.med.UTT = [squeeze(rmsdlens(10, :, 2, ix(2),1))' squeeze(rmsdlens(10, :, 2, ix(2),2))'];

LEN.ramp.fast.UT = [squeeze(rmsdlens(10, :, 3, ix(1),1))' squeeze(rmsdlens(10, :, 3, ix(1),2))'];
LEN.ramp.fast.TT = [squeeze(rmsdlens(10, :, 3, ix(3),1))' squeeze(rmsdlens(10, :, 3, ix(3),2))'];
LEN.ramp.fast.UTT = [squeeze(rmsdlens(10, :, 3, ix(2),1))' squeeze(rmsdlens(10, :, 3, ix(2),2))'];

LEN.ramp.asym.UT = [squeeze(rmsdlens(10, :, 4, ix(1),1))' squeeze(rmsdlens(10, :, 4, ix(1),2))'];
LEN.ramp.asym.TT = [squeeze(rmsdlens(10, :, 4, ix(3),1))' squeeze(rmsdlens(10, :, 4, ix(3),2))'];
LEN.ramp.asym.UTT = [squeeze(rmsdlens(10, :, 4, ix(2),1))' squeeze(rmsdlens(10, :, 4, ix(2),2))'];


LEN.pas.slow.UT = squeeze(mslen(ip(1),:,1))';
LEN.pas.slow.TT = squeeze(mslen(ip(2),:,1))';
LEN.pas.slow.UTT = squeeze(mslen(ip(3),:,1))';

LEN.pas.med.UT = squeeze(mslen(ip(1),:,2))';
LEN.pas.med.TT = squeeze(mslen(ip(2),:,2))';
LEN.pas.med.UTT = squeeze(mslen(ip(3),:,2))';

LEN.pas.fast.UT = squeeze(mslen(ip(1),:,3))';
LEN.pas.fast.TT = squeeze(mslen(ip(2),:,3))';
LEN.pas.fast.UTT = squeeze(mslen(ip(3),:,3))';


%% sinusoidal trials
algo_1 = 'UTT';
algo_2 = 'TT';

y = LEN; % PEN or LEN

close all
modelString = 'Acc ~ (Cycle*Algo) + (Range*Algo) + (1|Subject)';

Acc = [y.sine.low.(algo_1) y.sine.low.(algo_2) y.sine.high.(algo_1) y.sine.high.(algo_2)];
Algo = repmat({algo_1, algo_1,algo_2,algo_2},8,2);
Subject = repmat((1:8)',1,8);
Cycle = repmat([5 120], 8,4);

Range = [repmat({'low'},8,4), repmat({'high'},8,4)];

N = numel(Cycle);

T = table(reshape(Acc,N,1), reshape(Algo,N,1), reshape(Subject,N,1), reshape(Cycle,N,1), reshape(Range, N, 1));
T.Properties.VariableNames = {'Acc' 'Algo' 'Subject','Cycle','Range'};


% do linear regression!
LM_sine = fitlme(T,modelString)

%% symmetrical ramp trials
% algo_1 = 'UTT';
% algo_2 = 'UT';

% y = LEN; % PEN or LEN
close all
modelString = 'Acc ~ (Quality*Algo) + (Speed*Algo) + (1|Subject)';

Acc = [y.ramp.slow.(algo_1) y.ramp.slow.(algo_2) y.ramp.med.(algo_1) y.ramp.med.(algo_2) y.ramp.fast.(algo_1) y.ramp.fast.(algo_2)];
Algo = repmat({algo_1, algo_1,algo_2,algo_2},8,3);
Subject = repmat((1:8)',1,size(Acc,2));
Quality = repmat({'low','high'}, 8,size(Acc,2)/2);
Speed = [repmat(20,8,4) repmat(33,8,4) repmat(100,8,4)];

N = numel(Speed);

T = table(reshape(Acc,N,1), reshape(Algo,N,1), reshape(Subject,N,1), reshape(Quality,N,1), reshape(Speed, N, 1));
T.Properties.VariableNames = {'Acc' 'Algo' 'Subject','Quality','Speed'};


% do linear regression!
LM_symm = fitlme(T,modelString)

%% asymmetrical ramp trial
% algo_1 = 'UTT';
% algo_2 = 'UT';

% y = LEN; % PEN or LEN

close all
modelString = 'Acc ~ (Quality*Algo) + (1|Subject)';

Acc = [y.ramp.asym.(algo_1) y.ramp.asym.(algo_2)];
Algo = repmat({algo_1, algo_1,algo_2,algo_2},8,1);
Subject = repmat((1:8)',1,size(Acc,2));
Quality = repmat({'low','high'}, 8,size(Acc,2)/2);

N = numel(Quality);

T = table(reshape(Acc,N,1), reshape(Algo,N,1), reshape(Subject,N,1), reshape(Quality,N,1));
T.Properties.VariableNames = {'Acc' 'Algo' 'Subject','Quality'};


% do linear regression!
LM_asym = fitlme(T,modelString)

%% passive trials
close all
modelString = 'Acc ~ (Speed*Algo) + (1|Subject)';

% algo_1 = 'UTT';
% algo_2 = 'UT';

% y = LEN; % PEN or LEN

Acc = [y.pas.slow.(algo_1) y.pas.slow.(algo_2) y.pas.med.(algo_1) y.pas.med.(algo_2) y.pas.fast.(algo_1) y.pas.fast.(algo_2)];
Algo = repmat({algo_1,algo_2},8,3);
Subject = repmat((1:8)',1,size(Acc,2));
Speed = [repmat(5,8,2) repmat(30,8,2) repmat(120,8,2)];
N = numel(Speed);

T = table(reshape(Acc,N,1), reshape(Algo,N,1), reshape(Subject,N,1), reshape(Speed, N, 1));
T.Properties.VariableNames = {'Acc' 'Algo' 'Subject','Speed'};


% do linear regression!
LM_pass = fitlme(T,modelString)

%%

% cd('C:\Users\timvd\Documents\RUB_ultrasound_study\figures\data')
% save('data_for_statistics.mat','PEN','LEN')
