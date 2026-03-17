% Adaptive Intersection Maximization (AIM) based drift correction 
% Developed by Hongqiang Ma, University of Pittsburgh, January 2023
% Updated by Hongqiang Ma, University of Illinois Urbana-Champaign, February 2024
% Modified to use compiled C++ core for speed

function [Xpdc, Ypdc, driftX, driftY] = IntersectionMax(XList,YList,refXList,refYList,fID,trackNUM,trackInterval,imW,IntersectD)

% Call the compiled C++ core for the main tracking loop
[driftX, driftY] = IntersectionMax_core(XList, YList, refXList, refYList, fID, trackNUM, trackInterval, imW, IntersectD);

% Spline interpolation (kept in MATLAB - fast and simple)
drift_X = [2*driftX(1)-driftX(2) driftX 2*driftX(end)-driftX(end-1)]*IntersectD; 
drift_Y = [2*driftY(1)-driftY(2) driftY 2*driftY(end)-driftY(end-1)]*IntersectD;

driftX = interp1((-0.5:(trackNUM+0.5))*trackInterval,drift_X,1:trackNUM*trackInterval,'spline'); 
driftY = interp1((-0.5:(trackNUM+0.5))*trackInterval,drift_Y,1:trackNUM*trackInterval,'spline');

Xpdc = XList - driftX(fID)';
Ypdc = YList - driftY(fID)';

end
