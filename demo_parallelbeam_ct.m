clear ;
close all;

%%
restoredefaultpath();
addpath('./lib');
addpath('./util');
addpath('./matlab');
addpath('./mex');

%% Select the device to execute the simulation
% device: computed language to execute the simulation, choise=['matlab', 'clang'; (string)]
param.device    = 'matlab';

%% X-ray CT System parameter
% dAngle: Measure from 0 until the angle [degree; (float)]
% nView	: # of the views [element; (uint)]
% dView	: Gap between view_(k) - view_(k-1) [degree; (float)]
% DSO   : Distance from the Source to the Object    [mm; (float)]
% DSD 	: Distance from the Source to the Detector  [mm; (float)]

param.dAngle        = 360;
param.nView         = 360;
param.dView         = param.dAngle/param.nView;
param.DSO           = 400;
param.DSD           = 800;

%% X-ray detector parameter
% dDctY, dDctX              : Detector pitch [mm; (float)]
% nDctY, nDctX            	: Number of detector [element; (uint)]
% dOffsetDctY, dOffsetDctX  : Index of shifted detector [element; (float; +, -)]
% compute_filtering         : Filtering method, choise=['conv', 'fft'; (string)]
% ---------------------------------------------------
% '*DctY' parameters are only used when 3D CT system.
% ---------------------------------------------------

param.dDctY         = 1;
param.dDctX         = 0.7;

param.nDctY         = 1;
param.nDctX         = 400;

param.dOffsetDctY   = 0;
param.dOffsetDctX   = 0;

param.compute_filtering = 'conv';

%% Object parameter
% dImgY, dImgX, dImgZ	: Pixel resolution [mm; (float)]
% nImgY, nImgX, nImgZ	: Matrix size of image [element; (uint)]s
% dOffsetImgY, dOffsetImgX, dOffsetImgZ	: Index of shifted image [element; (float; +, -)]
% ----------------------------------------------------
% '*ImgZ' parameters are only used when 3D CT system. 
% ----------------------------------------------------

param.dImgY         = 1;
param.dImgX         = 1;
param.dImgZ         = 1;

param.nImgY         = 256;
param.nImgX         = 256;
param.nImgZ         = 1;

param.dOffsetImgY	= 0;
param.dOffsetImgX   = 0;
param.dOffsetImgZ   = 0;

%% Load image
fid                 = fopen(num2str(param.nImgY, 'input%d.raw'), 'rb');
input               = single(fread(fid, [param.nImgY, param.nImgX], 'single'));
fclose(fid);

%% Run implemented function
disp ('implementation ver.')

disp ('projection - implementation of Ch.3 Equation (3.5) & (3.6)');
tic;
prj         = projection(input, param);
toc;

disp ('filtering - implementation of Ch.3 Equation (3.21) & (3.29) & (3.30)');
tic;
prj_flt 	= filtering(prj, param);
toc;

disp ('backprojection - implementation of Ch.3 Equation (3.22)');
tic;
output   	= backprojection(prj_flt, param);
toc;

%% Display reconstruction images
wndImg      = [0, max(input(:))];

figure('name', 'parallel-beam CT'); colormap gray;
subplot(2,3,1);     imagesc(input, wndImg);     	
                axis image;     xlabel('X-axis'); 	ylabel('Y-axis');	title('ground truth');
subplot(2,3,2);     imagesc(output, wndImg);
                axis image;     xlabel('X-axis'); 	ylabel('Y-axis');	title(['reconstruction_{implementation, # of view = ' num2str(param.nView) '}']);
subplot(2,3,3);     imagesc(input - output);
                axis image;     xlabel('X-axis'); 	ylabel('Y-axis');	title({'difference', 'ground truch - reconstruction'});
subplot(2,3,[4,6]); imagesc(prj);  
                    xlabel(['Angle ( \Delta\theta : ' num2str(param.dView) ' \circ )']);	ylabel('Detector');	title(['sinogram_{implementation, # of view = ' num2str(param.nView) '}']);
                    ax              = gca;
                    ax.XTick        = linspace(0, param.nView, 10);
                    ax.XTickLabel	= linspace(0, param.dAngle, 10);