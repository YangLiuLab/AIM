function save_imSR(X, Y, F, Drift, fname, imSize, zoom)
% SAVE_IMSR  Render and save a super-resolution image as 16-bit TIFF
%
% Uses sparse matrix + strip writing + LZW compression.

    nDrift = size(Drift, 1);
    Fc = F(:);
    Fc(Fc > nDrift) = nDrift;
    Fc(Fc < 1) = 1;

    Xc = X(:) - Drift(Fc, 1);
    Yc = Y(:) - Drift(Fc, 2);

    imgSz = imSize * zoom;
    xi = round(Xc * zoom);
    yi = round(Yc * zoom);

    valid = (xi >= 1) & (xi <= imgSz) & (yi >= 1) & (yi <= imgSz);
    xi = xi(valid);
    yi = yi(valid);

    imSR = sparse(double(yi), double(xi), 1, imgSz, imgSz);

    outFile = [fname '.tif'];
    stripRows = 256;
    tobj = Tiff(outFile, 'w');

    tobj.setTag('Photometric', Tiff.Photometric.MinIsBlack);
    tobj.setTag('ImageLength', imgSz);
    tobj.setTag('ImageWidth', imgSz);
    tobj.setTag('BitsPerSample', 16);
    tobj.setTag('SamplesPerPixel', 1);
    tobj.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    tobj.setTag('RowsPerStrip', stripRows);
    tobj.setTag('Compression', Tiff.Compression.LZW);  % ~30x smaller for sparse data
    tobj.setTag('SampleFormat', Tiff.SampleFormat.UInt);

    nStrips = ceil(imgSz / stripRows);
    for s = 1:nStrips
        r1 = (s-1) * stripRows + 1;
        r2 = min(s * stripRows, imgSz);
        strip = uint16(full(imSR(r1:r2, :)));
        tobj.writeEncodedStrip(s, strip);
    end

    tobj.close();
    fprintf('Saved: %s (%dx%d)\n', outFile, imgSz, imgSz);
end