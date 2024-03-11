function [F,X,Y,Z,driftX,driftY,driftZ] = simulationSMLM(driftRMS,frameNUM,imSize,density,precision)


driftSTD = driftRMS*frameNUM/400;
drift_xyz = normrnd(0, driftSTD,400,3);
drift_xyz = cumsum(drift_xyz,1);
drift_xyz = drift_xyz - drift_xyz(1,:);
ratio = (frameNUM-1)/(400-1);

driftX = interp1(1:ratio:frameNUM,drift_xyz(:,1) - drift_xyz(1,1),1:frameNUM,'spline');
driftY = interp1(1:ratio:frameNUM,drift_xyz(:,2) - drift_xyz(1,2),1:frameNUM,'spline');
driftZ = interp1(1:ratio:frameNUM,drift_xyz(:,3) - drift_xyz(1,3),1:frameNUM,'spline');

d = 5*precision;
dotN = round(0.5/d);
mNUM = round(imSize*imSize*density/100);
patNUM = round(200000*density*imSize*imSize/2048/2048);

x = zeros(patNUM,dotN*dotN);
y = zeros(patNUM,dotN*dotN);
z = zeros(patNUM,dotN*dotN);

for n=1:patNUM
    x0 = rand()*(imSize-50)+25;
    y0 = rand()*(imSize-50)+25;
    z0 = rand()*8 - 4;
    a0 = rand()*pi*2;
    
    for r=0:(dotN-1)
        x1 = x0 + r*d*sin(a0);
        y1 = y0 + r*d*cos(a0);
        for c=1:dotN
            x(n,r*dotN+c) = x1 + (c-1)*d*cos(2*pi-a0);
            y(n,r*dotN+c) = y1 + (c-1)*d*sin(2*pi-a0);
            z(n,r*dotN+c) = z0;
        end
    end
end

F = zeros(mNUM*frameNUM,1);
X = zeros(mNUM*frameNUM,1);
Y = zeros(mNUM*frameNUM,1);
Z = zeros(mNUM*frameNUM,1);

for f = 1:frameNUM
    dx = driftX(f);
    dy = driftY(f);
    dz = driftZ(f);
    
    ID = ((f-1)*mNUM+1):((f-1)*mNUM + mNUM);
    
    F(ID) = f;
    id1 = round(rand(mNUM,1)*(patNUM-0.1) + 0.501)-1;
    id2 = round(rand(mNUM,1)*(dotN*dotN-0.1) + 0.501);
    id = id1*dotN*dotN+id2;
    X(ID) = x(id) + dx  + normrnd(0,precision,mNUM,1);
    Y(ID) = y(id) + dy  + normrnd(0,precision,mNUM,1);
    Z(ID) = z(id) + dz  + normrnd(0,precision,mNUM,1);  
end


end