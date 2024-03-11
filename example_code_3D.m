% Example code
% This code compares the performance of drift correction for AIM, RCC and
% DME using 3D localization coordinates of simulated data or experimental data of microtubules.
% This code was tested using MATLAB 2020 and 2021.

clc
clear
close all
warning('off')
addpath(genpath('./AIM'))
addpath(genpath('./DME_RCC'))
addpath(genpath('./Data'))

%% load experimental data
% fname = 'Microtublue_3d.mat';
% load(fname)

%% load simulation data
driftRMS = 0.002; % pixels
frameNUM = 20000; % frames
imSize = 2048; % pixels
density = 0.03; % number of localized emitters per um^2
precision = 0.02; % pixels

fname = 'simulationSMLM.mat';
[F,X,Y,Z,driftXT,driftYT,driftZT] = simulationSMLM(driftRMS,frameNUM,imSize,density,precision);
% save(fname,'F','X','Y','Z','driftXT','driftYT','driftZT');

%% data orgnization
dimensions = 3;

Localizations(:,1) = F;
Localizations(:,2) = X;
Localizations(:,3) = Y;
Localizations(:,4) = Z;

%% AIM drift correction
trackInterval = 20; % time interval for drift tracking, Unit: frames 
t_start = tic;
[LocAIM, AIM_Drift] = AIM(Localizations, trackInterval);
AIM_time = toc(t_start)


%% RCC computation
sigma = 1;
timebins = 10; % total number of drift estimation  
zoom = 5;
t_start = tic;
RCC_Drift = rcc3D(Localizations(:,2:4), F, timebins, zoom, sigma, 0);
RCC_time = toc(t_start)


%% DME computation
coarse_est = true;                      % Coarse drift estimation (bool)
precision_est = true;                  % Precision estimation (bool)
coarse_frames_per_bin = int32(10);      % Number of bins for coarse est. (int32)
framesperbin = int32(2); % int32((segSize)) % Number of frames per bin (int32)
maxneighbors_coarse = int32(1000);   % Max neighbors for coarse and precision est. (int32)
maxneighbors_regular = int32(1000);     % Max neighbors for regular est. (int32)
coarseSigma= single([0.02,0.02,0.02]);     % Localization precision for coarse estimation (single/float)                     
max_iter_coarse = int32(1000);          % Max iterations coarse est. (int32)
max_iter = int32(10000);                % Max iterations (int32)
gradientstep = single(1e-6);            % Gradient (single/float)
crlb = repmat([0.02 0.02 0.02], length(Localizations), 1);
t_start = tic;
RCC_Drift = rcc3D(Localizations(:,2:4), F, timebins, zoom, sigma, 0);
drift = dme_estimate(Localizations(:,2:4), F, crlb, RCC_Drift, 0, coarse_frames_per_bin, ...
    framesperbin, maxneighbors_coarse, maxneighbors_regular, coarseSigma, max_iter_coarse,max_iter, gradientstep, precision_est );
DME_Drift(:,1) = drift(:,1) - drift(1,1);
DME_Drift(:,2) = drift(:,2) - drift(1,2);
DME_Drift(:,3) = drift(:,3) - drift(1,3);
DME_time = toc(t_start)


%% save all data
save([fname(1:end-4) '_compare_results.mat'],'F','X','Y','Z', 'AIM_Drift', 'AIM_time', 'RCC_Drift', 'RCC_time', 'DME_Drift', 'DME_time')

%% precision
AIM_X_precision = std(driftXT-AIM_Drift(:,1)');
AIM_Y_precision = std(driftYT-AIM_Drift(:,2)');
AIM_Z_precision = std(driftZT-AIM_Drift(:,3)');
RCC_X_precision = std(driftXT-RCC_Drift(:,1)');
RCC_Y_precision = std(driftYT-RCC_Drift(:,2)');
RCC_Z_precision = std(driftZT-RCC_Drift(:,3)');
DME_X_precision = std(driftXT-DME_Drift(:,1)');
DME_Y_precision = std(driftYT-DME_Drift(:,2)');
DME_Z_precision = std(driftZT-DME_Drift(:,3)');

figure(1)
hold on
plot(100*driftXT,'k')
plot(100*AIM_Drift(:,1),'r')
plot(100*DME_Drift(:,1),'b')
plot(100*RCC_Drift(:,1),'g')
xlabel('Frame (100fps)')
ylabel('X drift (nm)')
grid
box

figure(2)
hold on
plot(100*driftYT,'k')
plot(100*AIM_Drift(:,2),'r')
plot(100*DME_Drift(:,2),'b')
plot(100*RCC_Drift(:,2),'g')
xlabel('Frame (100fps)')
ylabel('Y drift (nm)')
grid
box


figure(3)
hold on
plot(100*driftZT,'k')
plot(100*AIM_Drift(:,3),'r')
plot(100*DME_Drift(:,3),'b')
plot(100*RCC_Drift(:,3),'g')
xlabel('Frame (100fps)')
ylabel('Z drift (nm)')
grid
box

