% Adaptive Intersection Maximization (AIM) based drift correction 
% Develped by Hongqiang Ma, University of Pittsburgh, January 2023
% Updated by Hongqiang Ma, University of Illinois at Urbana Champaign, February 2024
% function: [Zpdc, driftZ] = IntersectionMaxZ(XList,YList,ZList,refXList,refYList,refZList,fID,trackNUM,trackInterval,imW,IntersectD)
% Input: 
%   XList: original full x position list 
%   YList: original full y position list
%   ZList: original full z position list
%   refXList: x position list in refrence subset
%   refYList: y position list in refrence subset
%   refZList: z position list in refrence subset
%   fID: frame index
%   trackNUM: segmentation number with specific tracking interval
%   trackInterval: tracking interval, unit: frames
%   imW: image width
%   IntersectD: intersection distance, 20nm
% Output:
%   Zpdc: drift corrected z position list 
%   driftZ: estimated drift in z dimension
%
% Contact information: mhq@illinois.edu 
%
function [Zpdc, driftZ] = IntersectionMaxZ(XList,YList,ZList,refXList,refYList,refZList,fID,trackNUM,trackInterval,imW,IntersectD)

pList = round(ZList/IntersectD)*imW/IntersectD*imW/IntersectD + round(YList/IntersectD)*imW/IntersectD + round(XList/IntersectD);
refList = round(refZList/IntersectD)*imW/IntersectD*imW/IntersectD + round(refYList/IntersectD)*imW/IntersectD + round(refXList/IntersectD);
roiR = 3;
[Vxyz0,Pxyz0] = groupcounts(refList);
refz = 0;
driftZ = zeros(1,trackNUM);
ROI_size = 2*roiR + 1;
for s = 1:trackNUM
    pList1 = pList((fID>(s-1)*trackInterval)&(fID<=s*trackInterval));
    [Vxyz1,Pxyz1] = groupcounts(pList1);
    
	sft = round(refz) * imW / IntersectD * imW / IntersectD;	
	Pxyz1 = Pxyz1 + sft;
    % count the number of intersected localizations
	ROIcc = PointIntersect1D(Pxyz0, Vxyz0, length(Vxyz0), Pxyz1, Vxyz1, length(Vxyz1), imW / IntersectD, roiR);

    fft_values = fft(ROIcc);
    angZ = angle(fft_values(2));
    angZ=angZ-2*pi*(angZ>0);
    PZ = (abs(angZ)/(2*pi/ROI_size) + 1) - (ROI_size+1)/2;
    
    refz = round(refz) + (PZ);
    driftZ(s) = -(refz);    
end

drift_Z = [2*driftZ(1)-driftZ(2) driftZ 2*driftZ(end)-driftZ(end-1)]*IntersectD;

driftZ = interp1((-0.5:(trackNUM+0.5))*trackInterval,drift_Z,1:trackNUM*trackInterval,'spline');

Zpdc = ZList - driftZ(fID)';

end
