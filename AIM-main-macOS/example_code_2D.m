%% AIM Drift Correction
clc
clear
try, close all; catch, end
warning('off')
addpath(genpath('./AIM'))
addpath(genpath('./Data'))

%% Load experimental data
fname = 'Origami_PAINT.mat';
load(fname);
imSize = 2048;
render_zoom = 20;

%% Organize data
Localizations(:,1) = F;
Localizations(:,2) = X;
Localizations(:,3) = Y;

%% AIM drift correction
trackInterval = 50;
t_start = tic;
[LocAIM, AIM_Drift] = AIM(Localizations, trackInterval);
AIM_time = toc(t_start)

%% Save results
save([fname(1:end-4) '_AIM.mat'], 'LocAIM', 'AIM_Drift', 'AIM_time');
save_imSR(X, Y, F, AIM_Drift*0, [fname(1:end-4) '_RAW'], imSize, render_zoom);
save_imSR(X, Y, F, AIM_Drift, [fname(1:end-4) '_AIM'], imSize, render_zoom);

%% Plot drift
figure(1)
plot(100*AIM_Drift(:,1), 'r')
xlabel('Frame (100fps)')
ylabel('X drift (nm)')
title(sprintf('AIM Time: %.2f s', AIM_time))
grid on; box on

figure(2)
plot(100*AIM_Drift(:,2), 'r')
xlabel('Frame (100fps)')
ylabel('Y drift (nm)')
title(sprintf('AIM Time: %.2f s', AIM_time))
grid on; box on