clear all; close all; clc

% cd('C:\Users\timvd\Documents')
% v = VideoReader('p6_isokinetic030_01_push_01.mp4');

cd('C:\Users\timvd\Documents\RUB\Test_0211\videos\high line density')
v = VideoReader('high_high01_short.mp4');

Im = read(v);
do_flip = 0;

%%
[Imc ROI] = imcrop(Im(:,:,1,1));

%%

% Imc = Im(130:330,215:915,1,:);
[n,m] = size(Imc);
load('parms.mat')
parms.apo.apox = round(linspace(parms.apo.apomargin, m-parms.apo.apomargin, parms.apo.napo));

[fas_thres, super_obj, deep_obj] = filter_usimage(Imc,parms);

%% elipsoid
% Create ellipse
r(1) = size(fas_thres,1)/2;
r(2) = size(fas_thres,2)/2;
th = linspace(0,2*pi);
% re = [r(1), (r(2)*parms.w_ellipse_rel)];

% Ellipse boundary
xc = r(2) + r(2)*cos(th) ; 
yc = r(1) + r(1)*sin(th); 

[nx,ny] = size(fas_thres);
[X,Y] = meshgrid(1:ny,1:nx) ;
idx = inpolygon(X(:),Y(:),xc',yc);

% Mask
parms.fas.Emask = zeros(size(fas_thres));
parms.fas.Emask(idx) = 1;

parms.fas.Emask_radius = r;

%%

if ishandle(2), close(2); end; figure(2)
subplot(211)
h = line('xdata',[],'ydata',[]);
g = line('xdata',[],'ydata',[]);
l = line('xdata',[],'ydata',[]);

subplot(212);
k = line('xdata',[],'ydata',[]);
m = line('xdata',[],'ydata',[]);

parms.fas.show = 0;
parms.fas.thetares = .5;

parms.fas.range = [1 89];

for i = 1:size(Im,4)
    
    
    Imcrop = imcrop(Im(:,:,1,i), ROI);
    
    if do_flip 
        Imcrop = flip(Imcrop,2);
    end
    
    fas_filt = FrangiFilter2D(double(Imcrop), parms.fas.frangi);  % filtered fascicle image

    % Fascicle image
    fas_thres = imbinarize(fas_filt,'adaptive','sensitivity', parms.fas.th); % threshold
    
    [alphas, w,hmat,gamma] = dohough(fas_thres,parms.fas);
    hmat_mean = movmean(max(hmat),10);
    
    if i == 1
        alpha(i) = weightedMedian(alphas,w);
%         alpha(i) = 2;
    end
    
    alpha2(i) = weightedMedian(alphas,w);
    
        
    gm = gmdistribution(alpha(i)+.1,2);
    p = pdf(gm, gamma(:));
    
    hslope = grad5(hmat_mean(:), mean(diff(gamma)));
    
    hslopei = interp1(gamma, hslope, alpha(i));
    
    beta = atan2d(hslopei,1);
    
    title(num2str((sind(beta))))
   

    
    alpha(i+1) = alpha(i) + sind(beta);
    % slope
   
%     alpha = weightedMedian(alphas,w);
    
    subplot(211)
    ylim([200 500])
    set(h, 'xdata', gamma, 'ydata',hmat_mean ,'linewidth',2);
    set(g, 'xdata', gamma, 'ydata', max(hmat),'linestyle','--'); 
    
    set(l,'xdata', alpha(i), 'ydata', interp1(gamma, hmat_mean, alpha(i)), 'marker','o','color','red')
    
    
    subplot(212);
    set(k,'xdata', gamma, 'ydata', hslope)
    set(m,'xdata', alpha(i), 'ydata', hslopei, 'marker','o','color','red')
    
    drawnow
    
%     pause
end
    
%%
if ishandle(3), close(3); end
figure(3);
plot(alpha); hold on
plot(alpha2)

plot(movmean(alpha,5))