% Example code
    % This code compares the performance of drift correction for AIM, RCC and
    % DME using 2D localization coordinates of experimental data of DNA origami or simulated data.
	% This code was tested using MATLAB 2020 and 2021.
clc
clear
close all
warning('off')
addpath(genpath('./AIM'))
addpath(genpath('./DME_RCC'))
addpath(genpath('./Data'))

%% load experimental data
% fname = 'Origami_PAINT.mat'; 
% load(fname)

%% Simulation data
driftRMS = 0.002;
frameNUM = 20000;
imSize = 2048;
render_zoom = 20;
density = 0.005; % number of localized emitters per um^2
precision = 0.02;
[F,X,Y,Z,driftXT,driftYT,driftZT] = simulationSMLM(driftRMS,frameNUM,imSize,density,precision);
save('simulationSMLM.mat','F','X','Y','Z','driftXT','driftYT','driftZT');

%% data orgnization
dimensions = 2;

Localizations(:,1) = F;
Localizations(:,2) = X;
Localizations(:,3) = Y;

%% AIM drift correction
trackInterval = 50; % time interval for drift tracking, Unit: frames 
t_start = tic;
[LocAIM, AIM_Drift] = AIM(Localizations, trackInterval);
AIM_time = toc(t_start)


%% RCC computation
sigma = 1;
timebins = 10; % total number of drift estimation  
zoom = 5;
t_start = tic;
RCC_Drift = rcc(Localizations(:,2:3), F, timebins, zoom, sigma, 0);
RCC_time = toc(t_start)


%% DME computation
coarse_est = true;                      % Coarse drift estimation (bool)
precision_est = true;                  % Precision estimation (bool)
coarse_frames_per_bin = int32(10);      % Number of bins for coarse est. (int32)
framesperbin = int32(2); % int32((segSize)) % Number of frames per bin (int32)
maxneighbors_coarse = int32(1000);   % Max neighbors for coarse and precision est. (int32)
maxneighbors_regular = int32(1000);     % Max neighbors for regular est. (int32)
coarseSigma= single([0.1,0.1]);     % Localization precision for coarse estimation (single/float)                     
max_iter_coarse = int32(1000);          % Max iterations coarse est. (int32)
max_iter = int32(10000);                % Max iterations (int32)
gradientstep = single(1e-6);            % Gradient (single/float)
crlb = repmat([0.1 0.1], length(Localizations), 1);
t_start = tic;
RCC_Drift = rcc(Localizations(:,2:3), F, timebins, zoom, sigma, 0);
drift = dme_estimate(Localizations(:,2:3), F, crlb, RCC_Drift, 0, coarse_frames_per_bin, ...
    framesperbin, maxneighbors_coarse, maxneighbors_regular, coarseSigma, max_iter_coarse,max_iter, gradientstep, precision_est);
DME_Drift(:,1) = drift(:,1) - drift(1,1);
DME_Drift(:,2) = drift(:,2) - drift(1,2);
DME_time = toc(t_start)

%% save all data
save([fname(1:end-4) '_compare_results.mat'],'F','X','Y', 'AIM_Drift', 'AIM_time', 'RCC_Drift', 'RCC_time', 'DME_Drift', 'DME_time')

%% Image Saving
save_imSR(X,Y,F,AIM_Drift*0,[fname(1:end-4) '_RAW'],imSize,render_zoom);
save_imSR(X,Y,F,AIM_Drift,[fname(1:end-4) '_AIM'],imSize,render_zoom);
save_imSR(X,Y,F,RCC_Drift,[fname(1:end-4) '_RCC'],imSize,render_zoom);
save_imSR(X,Y,F,DME_Drift,[fname(1:end-4) '_DME'],imSize,render_zoom);

figure(1)
hold on
plot(100*RCC_Drift(:,1),'g')
plot(100*DME_Drift(:,1),'b')
plot(100*AIM_Drift(:,1),'r')
xlabel('Frame (100fps)')
ylabel('X drift (nm)')
grid
box
figure(2)
hold on
plot(100*RCC_Drift(:,2),'g')
plot(100*DME_Drift(:,2),'b')
plot(100*AIM_Drift(:,2),'r')
xlabel('Frame (100fps)')
ylabel('Y drift (nm)')
grid
box

