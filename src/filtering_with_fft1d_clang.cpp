/*
 * IF, datatype is double
 * mex -L.\lib -L.\include -llibfftw3-3.lib  filtering_with_fft1d_cpu.cpp
 *
 * ELSE IF, datatype is single
 * mex -L.\lib -L.\include -llibfftw3f-3.lib  filtering_with_fft1d_cpu.cpp
 *
 * ELSE IF, datatype is long
 * mex -L.\lib -L.\include -llibfftw3l-3.lib  filtering_with_fft1d_cpu.cpp
 *
 */

#include "mex.h"
#include "./include/util.h"
#include "./include/fftw3.h"
#include <math.h>

using namespace std;

void    runFourierTransform1d(float *pout, float *pin,
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
    runFourierTransform1d(pout, pin,
                        DSO, DSD, dView, nView,
                        pdDct, pnDct, pdOffsetDct,
                        pdImg, pnImg, pdOffsetImg);
    
    return ;
}

void    runFourierTransform1d(float *pout, float *pin,
                        float DSO, float DSD, float dView, int nView,
                        float *pdDct, int *pnDct, float *pdOffsetDct,
                        float *pdImg, int *pnImg, float *pdOffsetImg) {
    
    
    int     nDctExtX	= pow(2, ceil(log2(2.0f*pnDct[X])));
    float   pnDctExt[2] = {pnDct[Y], nDctExtX};
    
    float   *pker       = (float*)   malloc(pnDctExt[X] * sizeof(float));
    memset(pker, 0, pnDctExt[X] * sizeof(float));
    
    float   *pview      = (float*)   malloc(pnDctExt[X] * sizeof(float));
    memset(pview, 0, pnDctExt[X] * sizeof(float));
    
    generate_filter(pker, pdDct[X], pnDct[X]);
    
    /*
     *
     */  
    
//     fftwf_complex	*pker_ft	= (fftwf_complex*)  fftwf_malloc(pnDct[X] * sizeof(fftwf_complex));
    fftwf_complex	*pker_ft	= (fftwf_complex*)  fftwf_malloc(pnDctExt[X] * sizeof(fftwf_complex));
    memset(pker_ft, 0, pnDctExt[X] * sizeof(fftwf_complex));
    
//     fftwf_complex	*pview_ft	= (fftwf_complex*)  fftwf_malloc(pnDct[X] * sizeof(fftwf_complex));
    fftwf_complex	*pview_ft	= (fftwf_complex*)  fftwf_malloc(pnDctExt[X] * sizeof(fftwf_complex));
    memset(pview_ft, 0, pnDctExt[X] * sizeof(fftwf_complex));

    
	fftwf_plan p, pi;
	p = fftwf_plan_dft_r2c_1d(pnDctExt[X], pker, pker_ft, FFTW_ESTIMATE);    
    
    fftwf_execute(p); 
    
    /*
     *
     */
    
    memset(pker, 0, pnDctExt[X] * sizeof(float));
    
	p = fftwf_plan_dft_r2c_1d(pnDctExt[X], pker, pview_ft, FFTW_ESTIMATE);
    pi= fftwf_plan_dft_c2r_1d(pnDctExt[X], pview_ft, pview, FFTW_ESTIMATE);
    
    fftwf_complex pview_ft_;
    
    for (int iView = 0; iView < nView; iView++) {
        
        memcpy(pker, &pin[pnDct[X]*iView], pnDct[X] * sizeof(float));
        fftwf_execute(p);  
        
        for (int iDctX = 0; iDctX < pnDctExt[X]; iDctX++) {
            pview_ft_[REAL]         = pview_ft[iDctX][REAL];
            pview_ft_[IMAG]         = pview_ft[iDctX][IMAG];
            pview_ft[iDctX][REAL]   = 1.0f/pnDctExt[X]*(pview_ft_[REAL]*pker_ft[iDctX][REAL] - pview_ft_[IMAG]*pker_ft[iDctX][IMAG]);
            pview_ft[iDctX][IMAG]   = 1.0f/pnDctExt[X]*(pview_ft_[REAL]*pker_ft[iDctX][IMAG] + pview_ft_[IMAG]*pker_ft[iDctX][REAL]);
        }

        fftwf_execute(pi); 
        memcpy(&pout[pnDct[X]*iView], &pview[pnDct[X] - 1], pnDct[X] * sizeof(float));
    }
    
	fftwf_destroy_plan(p);
    fftwf_destroy_plan(pi);
    fftwf_free(pker_ft); 
    fftwf_free(pview_ft); 
    free(pker);
    free(pview);
        
    return ;
}