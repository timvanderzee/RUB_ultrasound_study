clear all; close all; clc

load('cycle_averages_ramps.mat','msdphis','msdlens')
rmsdphis = msdphis;
rmsdlens = msdlens;

load('cycle_averages_sines_v2.mat','msdphis','msdlens')
load('variability_passive.mat','mspen','mslen')

%% pennation
PEN.sine.low.UT = squeeze(msdphis([10, 119], :, 2, 1))';
PEN.sine.low.TT = squeeze(msdphis([10, 119], :, 2, 9))';
PEN.sine.low.UTT = squeeze(msdphis([10, 119], :, 2, 5))';

PEN.sine.high.UT = squeeze(msdphis([10, 119], :, 1, 1))';
PEN.sine.high.TT = squeeze(msdphis([10, 119], :, 1, 9))';
PEN.sine.high.UTT = squeeze(msdphis([10, 119], :, 1, 5))';



PEN.ramp.slow.UT = [squeeze(rmsdphis(10, :, 1, 1,1))' squeeze(rmsdphis(10, :, 1, 1,2))'];
PEN.ramp.slow.TT = [squeeze(rmsdphis(10, :, 1, 9,1))' squeeze(rmsdphis(10, :, 1, 9,2))'];
PEN.ramp.slow.UTT = [squeeze(rmsdphis(10, :, 1, 5,1))' squeeze(rmsdphis(10, :, 1, 5,2))'];

PEN.ramp.med.UT = [squeeze(rmsdphis(10, :, 2, 1,1))' squeeze(rmsdphis(10, :, 2, 1,2))'];
PEN.ramp.med.TT = [squeeze(rmsdphis(10, :, 2, 9,1))' squeeze(rmsdphis(10, :, 2, 9,2))'];
PEN.ramp.med.UTT = [squeeze(rmsdphis(10, :, 2, 5,1))' squeeze(rmsdphis(10, :, 2, 5,2))'];

PEN.ramp.fast.UT = [squeeze(rmsdphis(10, :, 3, 1,1))' squeeze(rmsdphis(10, :, 3, 1,2))'];
PEN.ramp.fast.TT = [squeeze(rmsdphis(10, :, 3, 9,1))' squeeze(rmsdphis(10, :, 3, 9,2))'];
PEN.ramp.fast.UTT = [squeeze(rmsdphis(10, :, 3, 5,1))' squeeze(rmsdphis(10, :, 3, 5,2))'];

PEN.ramp.asym.UT = [squeeze(rmsdphis(10, :, 4, 1,1))' squeeze(rmsdphis(10, :, 4, 1,2))'];
PEN.ramp.asym.TT = [squeeze(rmsdphis(10, :, 4, 9,1))' squeeze(rmsdphis(10, :, 4, 9,2))'];
PEN.ramp.asym.UTT = [squeeze(rmsdphis(10, :, 4, 5,1))' squeeze(rmsdphis(10, :, 4, 5,2))'];


PEN.pas.slow.UT = squeeze(mspen(1,:,1))';
PEN.pas.slow.TT = squeeze(mspen(2,:,1))';
PEN.pas.slow.UTT = squeeze(mspen(3,:,1))';

PEN.pas.med.UT = squeeze(mspen(1,:,2))';
PEN.pas.med.TT = squeeze(mspen(2,:,2))';
PEN.pas.med.UTT = squeeze(mspen(3,:,2))';

PEN.pas.fast.UT = squeeze(mspen(1,:,3))';
PEN.pas.fast.TT = squeeze(mspen(2,:,3))';
PEN.pas.fast.UTT = squeeze(mspen(3,:,3))';

%% length
LEN.sine.low.UT = squeeze(msdlens([10, 119], :, 2, 1))';
LEN.sine.low.TT = squeeze(msdlens([10, 119], :, 2, 9))';
LEN.sine.low.UTT = squeeze(msdlens([10, 119], :, 2, 5))';

LEN.sine.high.UT = squeeze(msdlens([10, 119], :, 1, 1))';
LEN.sine.high.TT = squeeze(msdlens([10, 119], :, 1, 9))';
LEN.sine.high.UTT = squeeze(msdlens([10, 119], :, 1, 5))';



LEN.ramp.slow.UT = [squeeze(rmsdlens(10, :, 1, 1,1))' squeeze(rmsdlens(10, :, 1, 1,2))'];
LEN.ramp.slow.TT = [squeeze(rmsdlens(10, :, 1, 9,1))' squeeze(rmsdlens(10, :, 1, 9,2))'];
LEN.ramp.slow.UTT = [squeeze(rmsdlens(10, :, 1, 5,1))' squeeze(rmsdlens(10, :, 1, 5,2))'];

LEN.ramp.med.UT = [squeeze(rmsdlens(10, :, 2, 1,1))' squeeze(rmsdlens(10, :, 2, 1,2))'];
LEN.ramp.med.TT = [squeeze(rmsdlens(10, :, 2, 9,1))' squeeze(rmsdlens(10, :, 2, 9,2))'];
LEN.ramp.med.UTT = [squeeze(rmsdlens(10, :, 2, 5,1))' squeeze(rmsdlens(10, :, 2, 5,2))'];

LEN.ramp.fast.UT = [squeeze(rmsdlens(10, :, 3, 1,1))' squeeze(rmsdlens(10, :, 3, 1,2))'];
LEN.ramp.fast.TT = [squeeze(rmsdlens(10, :, 3, 9,1))' squeeze(rmsdlens(10, :, 3, 9,2))'];
LEN.ramp.fast.UTT = [squeeze(rmsdlens(10, :, 3, 5,1))' squeeze(rmsdlens(10, :, 3, 5,2))'];

LEN.ramp.asym.UT = [squeeze(rmsdlens(10, :, 4, 1,1))' squeeze(rmsdlens(10, :, 4, 1,2))'];
LEN.ramp.asym.TT = [squeeze(rmsdlens(10, :, 4, 9,1))' squeeze(rmsdlens(10, :, 4, 9,2))'];
LEN.ramp.asym.UTT = [squeeze(rmsdlens(10, :, 4, 5,1))' squeeze(rmsdlens(10, :, 4, 5,2))'];


LEN.pas.slow.UT = squeeze(mslen(1,:,1))';
LEN.pas.slow.TT = squeeze(mslen(2,:,1))';
LEN.pas.slow.UTT = squeeze(mslen(3,:,1))';

LEN.pas.med.UT = squeeze(mslen(1,:,2))';
LEN.pas.med.TT = squeeze(mslen(2,:,2))';
LEN.pas.med.UTT = squeeze(mslen(3,:,2))';

LEN.pas.fast.UT = squeeze(mslen(1,:,3))';
LEN.pas.fast.TT = squeeze(mslen(2,:,3))';
LEN.pas.fast.UTT = squeeze(mslen(3,:,3))';

cd('C:\Users\timvd\Documents\RUB_ultrasound_study\figures\data')
save('data_for_statistics.mat','PEN','LEN')
