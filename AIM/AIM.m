% Adaptive Intersection Maximization (AIM) based drift correction 
% Develped by Hongqiang Ma, University of Pittsburgh, January 2023
% Updated by Hongqiang Ma, University of Illinois at Urbana Champaign, February 2024
% function: [LocAIM, Drift] = AIM(Localizations, trackInterval)
% Input: 
%   Localizations: Localization list contains frame index (:,1), x position (:,2), y position (:,3), and z position (:,4)
%   trackInterval: the tracking interval, unit: frames
% Output:
%   LocAIM: the drift corrected localization with the same format to the input localization
%   Drift: the drift estimated by AIM 
%
% Contact information: mhq@illinois.edu 
%

function [LocAIM, Drift] = AIM(Localizations, trackInterval)

Dim = size(Localizations);
IntersectD = 0.2;
% IntersectD(IntersectD<0.1) = 0.1;

% check the format of dataset, 2D or 3D
if(Dim(2)==3) % 2D dataset
    F = Localizations(:,1);
    X = Localizations(:,2);
    Y = Localizations(:,3);
elseif(Dim(2)==4) % 3D dataset
    F = Localizations(:,1);
    X = Localizations(:,2);
    Y = Localizations(:,3);
    Z = Localizations(:,4);
else
    disp('localiation format error ')
end


F = F-min(F)+1;
frameNUM = max(F); 
trackNUM = floor(frameNUM/trackInterval);
frameNUM = trackNUM*trackInterval;
F(F>frameNUM) = frameNUM;
imageWidth = round(max(X))+ 10;

disp('AIM drift correction . . .  ')


refXList = X(F<=trackInterval);
refYList = Y(F<=trackInterval);
% the 1st round of AIM, using the 1st subset (X,Y) as reference
[Xpdc, Ypdc, driftX1, driftY1] = IntersectionMax(X,Y,refXList,refYList,F,trackNUM,trackInterval,imageWidth,IntersectD);
% the 2nd round of AIM, using the entire subset (Xpdc,Ypdc) as reference
[Xpdc, Ypdc, driftX2, driftY2] = IntersectionMax(Xpdc,Ypdc,Xpdc,Ypdc,F,trackNUM,trackInterval,imageWidth,IntersectD);

LocAIM(:,1) = F;
LocAIM(:,2) = Xpdc;
LocAIM(:,3) = Ypdc;

Drift_x = driftX1 + driftX2;
Drift_y = driftY1 + driftY2;
Drift_x = Drift_x - Drift_x(1);
Drift_y = Drift_y - Drift_y(1);

Drift(:,1) = Drift_x;
Drift(:,2) = Drift_y;

% for 3D dataset, perform drift tracking on Z diemsnion
if(Dim(2)==4)
    refXList = Xpdc(F<=trackInterval);
    refYList = Ypdc(F<=trackInterval);
    refZList = Z(F<=trackInterval);
    [Zpdc, driftZ1] = IntersectionMaxZ(Xpdc,Ypdc,Z,refXList,refYList,refZList,F,trackNUM,trackInterval,imageWidth,IntersectD);
    [Zpdc, driftZ2] = IntersectionMaxZ(Xpdc,Ypdc,Zpdc,Xpdc,Ypdc,Zpdc,F,trackNUM,trackInterval,imageWidth,IntersectD);
    
    LocAIM(:,4) = Zpdc;
    
    Drift_Z = driftZ1 + driftZ2;
    Drift_Z = Drift_Z - Drift_Z(1);
    
    Drift(:,3) = Drift_Z;
end

disp('AIM complete ')
end