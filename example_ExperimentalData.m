%% Example code
    % This code performs drift correction with AIM on 2D or 3D localization coordinates of experimental data.
	% This code was tested using MATLAB 2020 and 2021.

clear;
clc
clear
close all
warning('off')
addpath(genpath('./AIM'))
addpath(genpath('./DME_RCC'))
addpath(genpath('./Data'))

%% Load experimental data and set parameters
fname = 'Origami_PAINT.mat'; % Position lists of localization coordinates with three variabls: F, X and Y
load(fname);
imSize = 2048; % specify the image size in number of pixels, e.g., 2048 x 2048
render_zoom = 20; % magnification in the final rendered image
%% Asign the experimental frame number and localization coordinates (x and y) to localizations variable
Localizations(:,1) = F; %frame_id
Localizations(:,2) = X; % x position of the localization coordinates
Localizations(:,3) = Y; % y position of the localization coordinates

%% AIM drift correction and parameter settings
trackInterval = 50; % time interval for drift tracking, Unit: frames 
t_start = tic;
[LocAIM, AIM_Drift] = AIM(Localizations, trackInterval);
AIM_time = toc(t_start)

%% Save all data
save([fname(1:end-4) '_AIM.mat'],'LocAIM', 'AIM_Drift');
save_imSR(X,Y,F,AIM_Drift*0,[fname(1:end-4) '_RAW'],imSize,render_zoom);
save_imSR(X,Y,F,AIM_Drift,[fname(1:end-4) '_AIM'],imSize,render_zoom);
