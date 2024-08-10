clear all; close all; clc

datafolder = % should end with "figshare"
codefolder = % shouldn end with "RUB_ultrasound_study"
addpath(genpath(codefolder));

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
