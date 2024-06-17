clear all; close all; clc

load('cycle_averages_ramps.mat','msdphis')
rmsdphis = msdphis;

load('cycle_averages_sines_v2.mat','msdphis')
load('variability_passive.mat','mspen')

%% save some summary data to do statistics
clear sine ramp pas
sine.low.UT = squeeze(msdphis([10, 119], :, 2, 1))';
sine.low.TT = squeeze(msdphis([10, 119], :, 2, 9))';
sine.low.UTT = squeeze(msdphis([10, 119], :, 2, 5))';

sine.high.UT = squeeze(msdphis([10, 119], :, 1, 1))';
sine.high.TT = squeeze(msdphis([10, 119], :, 1, 9))';
sine.high.UTT = squeeze(msdphis([10, 119], :, 1, 5))';



ramp.slow.UT = [squeeze(rmsdphis(10, :, 1, 1,1))' squeeze(rmsdphis(10, :, 1, 1,2))'];
ramp.slow.TT = [squeeze(rmsdphis(10, :, 1, 9,1))' squeeze(rmsdphis(10, :, 1, 9,2))'];
ramp.slow.UTT = [squeeze(rmsdphis(10, :, 1, 5,1))' squeeze(rmsdphis(10, :, 1, 5,2))'];

ramp.med.UT = [squeeze(rmsdphis(10, :, 2, 1,1))' squeeze(rmsdphis(10, :, 2, 1,2))'];
ramp.med.TT = [squeeze(rmsdphis(10, :, 2, 9,1))' squeeze(rmsdphis(10, :, 2, 9,2))'];
ramp.med.UTT = [squeeze(rmsdphis(10, :, 2, 5,1))' squeeze(rmsdphis(10, :, 2, 5,2))'];

ramp.fast.UT = [squeeze(rmsdphis(10, :, 3, 1,1))' squeeze(rmsdphis(10, :, 3, 1,2))'];
ramp.fast.TT = [squeeze(rmsdphis(10, :, 3, 9,1))' squeeze(rmsdphis(10, :, 3, 9,2))'];
ramp.fast.UTT = [squeeze(rmsdphis(10, :, 3, 5,1))' squeeze(rmsdphis(10, :, 3, 5,2))'];

ramp.asym.UT = [squeeze(rmsdphis(10, :, 4, 1,1))' squeeze(rmsdphis(10, :, 3, 4,2))'];
ramp.asym.TT = [squeeze(rmsdphis(10, :, 4, 9,1))' squeeze(rmsdphis(10, :, 3, 4,2))'];
ramp.asym.UTT = [squeeze(rmsdphis(10, :, 4, 5,1))' squeeze(rmsdphis(10, :, 3, 4,2))'];


pas.slow.UT = squeeze(mspen(1,:,1))';
pas.slow.TT = squeeze(mspen(2,:,1))';
pas.slow.UTT = squeeze(mspen(2,:,1))';

pas.med.UT = squeeze(mspen(1,:,2))';
pas.med.TT = squeeze(mspen(2,:,2))';
pas.med.UTT = squeeze(mspen(2,:,2))';

pas.fast.UT = squeeze(mspen(1,:,3))';
pas.fast.TT = squeeze(mspen(2,:,3))';
pas.fast.UTT = squeeze(mspen(2,:,3))';

cd('C:\Users\timvd\Documents\RUB_ultrasound_study\figures\data')
save('data_for_statistics.mat','sine','pas','ramp')
