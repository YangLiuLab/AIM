function save_imSR(X,Y,F,drift,fname,imSize,zoom)

Xd = X;
Yd = Y;

for n=1:round(length(X)*0.997)
    Xd(n) = X(n) - drift(F(n),1);
    Yd(n) = Y(n) - drift(F(n),2);
end

Xd(Xd>imSize-1) = imSize-1;
Xd(Xd<1) = 1;
Yd(Yd>imSize-1) = imSize-1;
Yd(Yd<1) = 1;

imSR = zeros(imSize*zoom,imSize*zoom);
for n=1:length(Xd)
    imSR(round(Xd(n)*zoom),round(Yd(n)*zoom)) = imSR(round(Xd(n)*zoom),round(Yd(n)*zoom)) + 1;
end
imwrite(uint16(imSR),[fname '_SR.tif'])

% roi = round(imSize*zoom/2)-4100:round(imSize*zoom/2)+4099;
%% Images needed for calculating FRC resolutions
% imSR1 = zeros(imSize*zoom,imSize*zoom);
% for n=1:round(length(Xd)/2)
%     imSR1(round(Xd(n)*zoom),round(Yd(n)*zoom)) = imSR1(round(Xd(n)*zoom),round(Yd(n)*zoom)) + 1;
% end
% imwrite(uint16(imSR1(roi,roi)),[fname '_FRC1.tif'])
% 
% imSR2 = zeros(imSize*zoom,imSize*zoom);
% for n=round(length(Xd)/2)+1:length(Xd)
%     imSR2(round(Xd(n)*zoom),round(Yd(n)*zoom)) = imSR2(round(Xd(n)*zoom),round(Yd(n)*zoom)) + 1;
% end
% imwrite(uint16(imSR2(roi,roi)),[fname '_FRC2.tif'])

%%
% xrange = round(imSize*zoom/2)-6000:round(imSize*zoom/2)+6000;
% yrange = round(imSize*zoom/2)-6000:round(imSize*zoom/2)+6000;
% imSR1_ = imSR1(xrange,yrange);
% imSR2_ = imSR2(xrange,yrange);
% [Ffrc, Vfrc] = FRC(imSR1_,imSR2_);
% figure
% plot(Ffrc,Vfrc)

end