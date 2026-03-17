%% AIM 3D Drift Correction
clc
clear
try, close all; catch, end
warning('off')
addpath(genpath('./AIM'))
addpath(genpath('./Data'))

%% Load experimental data
fname = 'Microtublue_3d.mat';
load(fname);

%% Organize data (3D)
Localizations(:,1) = F;
Localizations(:,2) = X;
Localizations(:,3) = Y;
Localizations(:,4) = Z;

%% AIM drift correction
trackInterval = 20;
t_start = tic;
[LocAIM, AIM_Drift] = AIM(Localizations, trackInterval);
AIM_time = toc(t_start)

%% Save results
save([fname(1:end-4) '_AIM.mat'], 'LocAIM', 'AIM_Drift', 'AIM_time');

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

figure(3)
plot(100*AIM_Drift(:,3), 'r')
xlabel('Frame (100fps)')
ylabel('Z drift (nm)')
title(sprintf('AIM Time: %.2f s', AIM_time))
grid on; box on