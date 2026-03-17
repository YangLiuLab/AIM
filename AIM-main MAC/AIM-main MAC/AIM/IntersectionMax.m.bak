% Adaptive Intersection Maximization (AIM) based drift correction 
% Developed by Hongqiang Ma, University of Pittsburgh, January 2023
% Updated by Hongqiang Ma, University of Illinois Urbana-Champaign, February 2024

% function: [Xpdc, Ypdc, driftX, driftY] = IntersectionMax(XList,YList,refXList,refYList,fID,trackNUM,trackInterval,imW,IntersectD)

% Input: 
%   XList: original full x position list 
%   YList: original full y position list
%   refXList: x position list in the reference subset
%   refYList: y position list in the reference subset
%   fID: image frame index
%   trackNUM: segmentation number with specific tracking interval
%   trackInterval: tracking interval, unit: frames
%   imW: image width
%   IntersectD: intersection distance, 20nm

% Output:
%   Xpdc: drift-corrected x position list 
%   Ypdc: drift-corrected y position list 
%   driftX: estimated drift in x dimension
%   driftY: estimated drift in y dimension
%
% Contact information: mhq@illinois.edu 
%

function [Xpdc, Ypdc, driftX, driftY] = IntersectionMax(XList,YList,refXList,refYList,fID,trackNUM,trackInterval,imW,IntersectD)

% transfer the 2D position to 1D (x+y*imW)
pList = round(YList/IntersectD)*imW/IntersectD + round(XList/IntersectD);
refList = round(refYList/IntersectD)*imW/IntersectD + round(refXList/IntersectD);
roiR = 3; % radius of the local search region
[Vxy0,Pxy0] = groupcounts(refList);  % sorting of the 1D localization position
refx = 0;
refy = 0;
ROI_size = 2*roiR+1;
num = 0;
array_sft = zeros(ROI_size*ROI_size,1);
array_idx = zeros(ROI_size*ROI_size,1);

% transfer 2D local search region to 1D
for r = -roiR:roiR
    for c = -roiR:roiR
        num = num + 1;
        array_sft(num) = r * imW / IntersectD + c;
		array_idx(num) = (r + roiR) * (roiR * 2 + 1) + (c + roiR);
    end
end

array_len = ROI_size*ROI_size;
driftX = zeros(1,trackNUM);
driftY = zeros(1,trackNUM);
for s = 2:trackNUM
    pList1 = pList((fID>(s-1)*trackInterval)&(fID<=s*trackInterval)); % extract localization list in each segmented subset 
    [Vxy1,Pxy1] = groupcounts(pList1); 
	sft = round(refy) * imW / IntersectD + round(refx);
	Pxy1 = Pxy1 + sft;
    % count the number of intersected localizations
	img_crr = PointIntersect2D(Pxy0, Vxy0, length(Vxy0), Pxy1, Vxy1, length(Vxy1), array_sft, array_idx, array_len, roiR);
	ROIcc = reshape(img_crr, 2*roiR+1, 2*roiR+1);
	
	% estimate the precise sub-pixel position of the peak with a fast FFT based single-molecule localization algorithm 
    fft_values = fft2(ROIcc');
    angX = angle(fft_values(1,2));
    angX=angX-2*pi*(angX>0);
    PX = (abs(angX)/(2*pi/ROI_size) + 1) - (ROI_size+1)/2;
    angY = angle(fft_values(2,1));
    angY=angY-2*pi*(angY>0);
    PY = (abs(angY)/(2*pi/ROI_size) + 1) - (ROI_size+1)/2;
   
	% update the relative drift reference for subsequent segmented subset
    refx = round(refx) + (PX);
    refy = round(refy) + (PY);
    driftX(s) = -(refx);
    driftY(s) = -(refy);    
end


% spline interpolation 
drift_X = [2*driftX(1)-driftX(2) driftX 2*driftX(end)-driftX(end-1)]*IntersectD; 
drift_Y = [2*driftY(1)-driftY(2) driftY 2*driftY(end)-driftY(end-1)]*IntersectD;

driftX = interp1((-0.5:(trackNUM+0.5))*trackInterval,drift_X,1:trackNUM*trackInterval,'spline'); 
driftY = interp1((-0.5:(trackNUM+0.5))*trackInterval,drift_Y,1:trackNUM*trackInterval,'spline');

Xpdc = XList - driftX(fID)';
Ypdc = YList - driftY(fID)';



end
