%% COMPILE the MEX FILES for projection, backprojection, and filtering using conv1d and fft1d
% If you are want to compile with fft1d,
% then fft3w library was downloaded from "http://www.fftw.org/download.html"
%
% Some details are described in
% "https://stackoverflow.com/questions/39675436/how-to-get-fftw-working-on-windows-for-dummies"

mex -outdir .\mex                                                   .\src\projection_clang.cpp
mex -outdir .\mex                                                   .\src\backprojection_clang.cpp
mex -outdir .\mex                                                   .\src\filtering_with_conv1d_clang.cpp
mex -outdir .\mex	-L.\src\lib	-L.\src\include	-llibfftw3f-3.lib	.\src\filtering_with_fft1d_clang.cpp
copyfile('.\src\include\libfftw3-3.dll', '.\mex');
copyfile('.\src\include\libfftw3f-3.dll', '.\mex');
copyfile('.\src\include\libfftw3l-3.dll', '.\mex');