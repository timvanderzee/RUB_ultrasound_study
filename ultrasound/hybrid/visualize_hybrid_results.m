clear all; close all; clc

filename = 'sine_020_hybrid_Test';

cd('C:\Users\timvd\Documents\RUB_ultrasound_study\ultrasound\hybrid\individual_subjects')
dates = dir(cd);

alphas = [0 .01 .1];

for j = 1:3
for i = 1:(length(dates)-2)
    load([filename,dates(i+2).name,'_alpha=',num2str(alphas(j)*100),'.mat'])
    
    phi = -atan2d(diff(ypos), diff(xpos));
    faslen = sqrt(diff(xpos).^2 + diff(ypos).^2);
    
    figure(1)
    subplot(4,2,i);
    plot(phi); hold on; box off
    xlabel('Sample #'); ylabel('\alpha (deg)')
    xlim([0 length(phi)]); title(dates(i+2).name)
    
    figure(2)
    subplot(4,2,i);
    plot(faslen); hold on; box off
    xlabel('Sample #'); ylabel('L_{fas} (pixels)')
     xlim([0 length(phi)]);title(dates(i+2).name)
end
end