#include "mex.h"
#include "./include/util.h"

#include <math.h>

using namespace std;

void    runProjection(float *pout, float *pin,
                        float DSO, float DSD, float dView, int nView,
                        float *pdDct, int *pnDct, float *pdOffsetDct,
                        float *pdImg, int *pnImg, float *pdOffsetImg);

/*
 * Please read the Ch.3 Image Reconstruction
 *
 * Implementation for projection operator based on Ch.3 Equation (3.5) & (3.6)
 * Projection operator is implemented using ray-driven method
*/

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
    
    plhs[0]                 = mxCreateNumericMatrix(pnDct[Y]*pnDct[X], nView, mxSINGLE_CLASS, mxREAL);
    float   *pout           = (float *) mxGetData(plhs[0]);
    
    /*
     *  Run projection operator
     */
    
    runProjection(pout, pin,
                    DSO, DSD, dView, nView,
                    pdDct, pnDct, pdOffsetDct,
                    pdImg, pnImg, pdOffsetImg);
    
    return ;
}

void    runProjection(float *pout, float *pin,
                        float DSO, float DSD, float dView, int nView,
                        float *pdDct, int *pnDct, float *pdOffsetDct,
                        float *pdImg, int *pnImg, float *pdOffsetImg) {
    
    /*
     *  Initialize projection parameters
     */
    float   pdSizeImg[2]	= {pnImg[Y]*pdImg[Y], pnImg[X]*pdImg[X]};
    
    float   dDiameter       = sqrt(pdSizeImg[Y]*pdSizeImg[Y] + pdSizeImg[X]*pdSizeImg[X]);
    float   dRadius         = dDiameter/2.0f;
    
    float   dSample         = min(pdImg[X], pdImg[Y]);
    int     nSample         = int(ceil(dDiameter/dSample));

    float	dBeta, dGamma;
    
    float   pdOriImg[2]     = {0.0f, 0.0f};
    float   pdOriNorDir[2]  = {0.0f, 0.0f};
    
    float   pdPosImg[2]     = {0.0f, 0.0f};
    float   pdNorDir[2]     = {0.0f, 0.0f};
    
    float   pnCurIdImg[2]   = {0.0f, 0.0f};
    
    float   out;
    
    // Ch.3 Equation (3.6)
    // Projection operator
    for (int iView = 0; iView < nView; iView++) {
        dBeta   = iView*dView * PI / 180.0f;
        
        // Rotation angle for geometry
        for (int iDctX = 0; iDctX < pnDct[X]; iDctX++) {    
            pdOriImg[Y]     = -dRadius;
            pdOriImg[X]     = id2pos(iDctX + pdOffsetDct[X], pdDct[X], pnDct[X]);
        
            dGamma          = 0.0f;
            
            // Normal vector of incident X-ray [mm, mm] 
            pdOriNorDir[Y]	= dSample*cos(dGamma);
            pdOriNorDir[X]	= dSample*sin(dGamma);
            
            // Ch.3 Equation (3.5)
            // Rotated X-ray source position [mm, mm]
            pdPosImg[Y] = -sin(dBeta)*pdOriImg[X] + cos(dBeta)*pdOriImg[Y];
            pdPosImg[X] = cos(dBeta)*pdOriImg[X] + sin(dBeta)*pdOriImg[Y];
            
            // Ch.3 Equation (3.5)
            // Rotated normal vector of incident X-ray [mm, mm]
            pdNorDir[Y] = -sin(dBeta)*pdOriNorDir[X] + cos(dBeta)*pdOriNorDir[Y];
            pdNorDir[X] = cos(dBeta)*pdOriNorDir[X] + sin(dBeta)*pdOriNorDir[Y];
            
            out             = 0.0f;
            
            for (int iSample = 0; iSample < nSample; iSample++) {
                pnCurIdImg[Y]	= pos2id(pdPosImg[Y] + iSample*pdNorDir[Y], pdImg[Y], pnImg[Y]) - pdOffsetImg[Y];
                pnCurIdImg[X]	= pos2id(pdPosImg[X] + iSample*pdNorDir[X], pdImg[X], pnImg[X]) - pdOffsetImg[X];
                
                out             += interpolation2d(pin, pnCurIdImg, pnImg);
            }
            
            pout[pnDct[X]*iView + iDctX]     = dSample*out;
        }
    }
    
    return ;
}