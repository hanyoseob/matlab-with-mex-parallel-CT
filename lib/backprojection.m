function out = backprojection(in, param, bfig)

if nargin < 3
    bfig = false;
end

switch param.device
    case 'matlab'
        out     = backprojection_matlab(in, param);
    case 'clang'
        out     = backprojection_clang(in, param.dView, param.nView, param.DSO, param.DSD, ...
                                       param.dDctY, param.dDctX, param.nDctY, param.nDctX, param.dOffsetDctY, param.dOffsetDctX, ...
                                       param.dImgY, param.dImgX, param.dImgZ, param.nImgY, param.nImgX, param.nImgZ, param.dOffsetImgY, param.dOffsetImgX, param.dOffsetImgZ);
    case 'gpu'
        disp('=========================================');
        out     = zeros(param.nImgY, param.nImgX, 'like', in);
end
                    
if bfig
    figure('name', 'backprojection'); colormap gray;
    imagesc(out);   title('BackProjection');
    xlabel('x-axis'); ylabel('y-axis');
    drawnow();
end

end
