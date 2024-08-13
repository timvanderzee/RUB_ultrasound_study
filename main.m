clear all; close all; clc

%% get paths
datafolder = uigetdir('','Find the folder where the two figshare subfolders were unzipped and select it')
codefolder = uigetdir('','Find the RUB_ultrasound_study folder and select it')
addpath(genpath(codefolder));

%% create missing folder
if ~exist([codefolder '\data\emg'], 'dir')
       mkdir([codefolder '\data\emg'])
end

%% torque
disp('Processing torque ...')
process_torque
create_target

%% EMG
disp('Processing EMG ...')
process_emg
process_emg_mvc

%% ultrasound
disp('Processing ultrasound ...')
process_passive
process_ramps
process_sinusoids

%% figures
disp('Creating figures ...')
fig3A
fig3B
fig4
fig5
fig6
fig7
figS1