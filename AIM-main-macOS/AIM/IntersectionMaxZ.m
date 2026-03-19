% Adaptive Intersection Maximization (AIM) based drift correction 
% Developed by Hongqiang Ma, University of Pittsburgh, January 2023
% Updated by Hongqiang Ma, University of Illinois Urbana-Champaign, February 2024
% Modified to use compiled C++ core for speed

function [Zpdc, driftZ] = IntersectionMaxZ(XList,YList,ZList,refXList,refYList,refZList,fID,trackNUM,trackInterval,imW,IntersectD)

% Call the compiled C++ core for the main tracking loop
driftZ = IntersectionMaxZ_core(XList, YList, ZList, refXList, refYList, refZList, fID, trackNUM, trackInterval, imW, IntersectD);

% Spline interpolation
drift_Z = [2*driftZ(1)-driftZ(2) driftZ 2*driftZ(end)-driftZ(end-1)]*IntersectD;  
driftZ = interp1((-0.5:(trackNUM+0.5))*trackInterval,drift_Z,1:trackNUM*trackInterval,'spline'); 

Zpdc = ZList - driftZ(fID)';

end
