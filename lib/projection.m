function pdY = projection(pdX, param, bfig)

if nargin < 3
    bfig = false;
end

switch param.device
    case 'matlab'
        pdY     = projection_matlab(pdX, param);
    case 'clang'
        pdY     = projection_clang(pdX, param.dView, param.nView, param.DSO, param.DSD, ...
                                param.dDctY, param.dDctX, param.nDctY, param.nDctX, param.dOffsetDctY, param.dOffsetDctX, ...
                                param.dImgY, param.dImgX, param.dImgZ, param.nImgY, param.nImgX, param.nImgZ, param.dOffsetImgY, param.dOffsetImgX, param.dOffsetImgZ);
    case 'gpu'
        disp('=========================================');
        pdY     = [];
end

if bfig
    figure('name', 'projection' ); colormap gray;
    imagesc(pdY);   title('Projection');
    xlabel('angle'); ylable('detector');
    drawnow();
end

end
