clear all; close all; clc

datafolder = 'C:\Users\u0167448\Desktop\data'% should contain folders downloaded from figshare
codefolder = 'C:\Users\u0167448\Documents\GitHub\RUB_ultrasound_study'% should end with "RUB_ultrasound_study"
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
