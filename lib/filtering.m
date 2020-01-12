function out = filtering(in, param)

switch param.device
    case 'matlab'
        out	= filtering_matlab(in, param);
    case 'clang'
        switch param.compute_filtering
            case 'conv'
                out = filtering_with_conv1d_clang(in, param.dView, param.nView, param.DSO, param.DSD, ...
                                param.dDctY, param.dDctX, param.nDctY, param.nDctX, param.dOffsetDctY, param.dOffsetDctX, ...
                                param.dImgY, param.dImgX, param.dImgZ, param.nImgY, param.nImgX, param.nImgZ, param.dOffsetImgY, param.dOffsetImgX, param.dOffsetImgZ);
            case 'fft'
                out = filtering_with_fft1d_clang(in, param.dView, param.nView, param.DSO, param.DSD, ...
                                param.dDctY, param.dDctX, param.nDctY, param.nDctX, param.dOffsetDctY, param.dOffsetDctX, ...
                                param.dImgY, param.dImgX, param.dImgZ, param.nImgY, param.nImgX, param.nImgZ, param.dOffsetImgY, param.dOffsetImgX, param.dOffsetImgZ);
        end
    case 'gpu'
        out = [];
end

end

