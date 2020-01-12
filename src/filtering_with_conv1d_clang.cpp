#include "mex.h"
#include "./include/util.h"

#include <math.h>

using namespace std;


void    runConvolution1d(float *pout, float *pin,
                            float DSO, float DSD, float dView, int nView,
                            float *pdDct, int *pnDct, float *pdOffsetDct,
                            float *pdImg, int *pnImg, float *pdOffsetImg);


void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float   *pin    = (float *) mxGetData(prhs[0]);
     
    /*
     *  X-ray CT System parameters
     *
     * dView	: Gap between view_(k) - view_(k-1) [degree (float)]
     * nView	: # of the views [element (uint)]
     * DSO      : Distance from the Source to the Object    [mm (float)]
     * DSD      : Distance from the Source to the Detector  [mm (float)]
     */
    
    float   dView           = (float)   mxGetScalar(prhs[1]);
    int     nView           = (int)     mxGetScalar(prhs[2]);
    float   DSO             = (float)   mxGetScalar(prhs[3]);
    float   DSD             = (float)   mxGetScalar(prhs[4]);
    
    /*
     *  X-ray CT Detector parameters
     *
     * pdDct[2]         : Detector pitch [mm (float)]
     * pnDct[2]         : Number of detector [element (uint)]
     * pdOffsetDct[2]	: Index of shifted detector [element (float; +, -)]
     * -----------------------------------------------------
     * '*Dct[Y]' parameters are only used when 3D CT system. 
     * -----------------------------------------------------
     */
    
    float   pdDct[2];
    pdDct[Y]        = (float)   mxGetScalar(prhs[5]);
    pdDct[X]        = (float)   mxGetScalar(prhs[6]);
    
    int     pnDct[2];
    pnDct[Y]        = (int)     mxGetScalar(prhs[7]);
    pnDct[X]        = (int)     mxGetScalar(prhs[8]);
    
    float   pdOffsetDct[2];
    pdOffsetDct[Y]	= (float)   mxGetScalar(prhs[9]);
    pdOffsetDct[X]	= (float)   mxGetScalar(prhs[10]);
    
    /*
     *  Image Object parameters
     *
     *  pdImg[3]        : The resolution of the element of the image voxel
     *  pnImg[3]        : The number of the element of the the image voxel
     *  pdOffsetImg[3]  : The offset (shift) from the centor of the image
     * ----------------------------------------------------
     * '*ImgZ' parameters are only used when 3D CT system. 
     * ----------------------------------------------------
     */
    
    float   pdImg[3];
    pdImg[Y]        = (float)   mxGetScalar(prhs[11]);
    pdImg[X]        = (float)   mxGetScalar(prhs[12]);
    pdImg[Z]        = (float)   mxGetScalar(prhs[13]);
    
    int     pnImg[3];
    pnImg[Y]        = (int)     mxGetScalar(prhs[14]);
    pnImg[X]        = (int)     mxGetScalar(prhs[15]);
    pnImg[Z]        = (int)     mxGetScalar(prhs[16]);
    
    float   pdOffsetImg[3];
    pdOffsetImg[Y]	= (float)   mxGetScalar(prhs[17]);
    pdOffsetImg[X]	= (float)   mxGetScalar(prhs[18]);
    pdOffsetImg[Z]	= (float)   mxGetScalar(prhs[19]);
    
    /*
     *  Sinogram data
     */
    plhs[0]                 = mxCreateNumericMatrix(pnDct[X], nView, mxSINGLE_CLASS, mxREAL);
    float   *pout           = (float *) mxGetData(plhs[0]);
    
    /*
     *  Run filtering operator
     */
    runConvolution1d(pout, pin,
                        DSO, DSD, dView, nView,
                        pdDct, pnDct, pdOffsetDct,
                        pdImg, pnImg, pdOffsetImg);
    
    return ;
}

void    runConvolution1d(float *pout, float *pin,
                        float DSO, float DSD, float dView, int nView,
                        float *pdDct, int *pnDct, float *pdOffsetDct,
                        float *pdImg, int *pnImg, float *pdOffsetImg) {
    
    float   pnDctExt[2] = {pnDct[Y], 2*pnDct[X] - 1};
    
    float   *pker       = (float*)   malloc(pnDctExt[X] * sizeof(float));
    memset(pker, 0, pnDctExt[X] * sizeof(float));
    
    generate_filter(pker, pdDct[X], pnDct[X]);
    
    for (int iView = 0; iView < nView; iView++) {
        convolution1d(&pout[pnDct[X]*iView], &pin[pnDct[X]*iView], pker, pnDct[X]);
    }
    
    free(pker);
    pker = 0;
        
    return ;
}