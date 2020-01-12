#include "mex.h"
#include "./include/util.h"

#include <math.h>

using namespace std;

void    runBackprojection(float *pout, float *pin,
                        float DSO, float DSD, float dView, int nView,
                        float *pdDct, int *pnDct, float *pdOffsetDct,
                        float *pdImg, int *pnImg, float *pdOffsetImg);

/*
 * Please read the Ch.3 Image Reconstruction
 *
 * Implementation for backprojection operator based on Ch.3 Equation (3.22)
 * Backprojection operator is implemented using pixel-driven method
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
     *  Image data
     */
    plhs[0]                 = mxCreateNumericMatrix(pnImg[Y], pnImg[X], mxSINGLE_CLASS, mxREAL);
    float   *pout           = (float *) mxGetData(plhs[0]);
    
    /*
     *  Run backprojection operator
     */
    runBackprojection(pout, pin,
                        DSO, DSD, dView, nView,
                        pdDct, pnDct, pdOffsetDct,
                        pdImg, pnImg, pdOffsetImg);
    
    return ;
}

void    runBackprojection(float *pout, float *pin,
                        float DSO, float DSD, float dView, int nView,
                        float *pdDct, int *pnDct, float *pdOffsetDct,
                        float *pdImg, int *pnImg, float *pdOffsetImg) {
    
    /*
     *  Initialize backprojection parameters
     */
    
    float   dRadius;
    
    float	dBeta, dGamma, dPhi;
    
    float   pdOriImg[2]     = {0.0f, 0.0f};
    
    float   pdPosImg[2]     = {0.0f, 0.0f}; 
    float   pdDist[2]       = {0.0f, 0.0f};
    
    float   pnCurIdDct[2]   = {0.0f, 0.0f};
    
    float   out;
    
    // Ch.3 Equation (3.22)
    // Backprojection operator
    for (int iView = 0; iView < nView; iView++) {
        
        // Rotation angle for geometry
        dBeta   = -iView*dView * PI / 180.0f;
        
        for (int iImgX = 0; iImgX < pnImg[X]; iImgX++) {
            
            // X position of object
            pdPosImg[X]     = id2pos(iImgX + pdOffsetImg[X], pdImg[X], pnImg[X]);
            
            for (int iImgY = 0; iImgY < pnImg[Y]; iImgY++) {
                
                // Y position of object
                pdPosImg[Y]     = id2pos(iImgY + pdOffsetImg[Y], pdImg[Y], pnImg[Y]);
                
                // Calculate a detector position based on Figure 3.34
                dRadius         = sqrt(pdPosImg[Y]*pdPosImg[Y] + pdPosImg[X]*pdPosImg[X]);
                dPhi            = atan2(pdPosImg[Y], pdPosImg[X]);
                
                // Detector position related with (X, Y) position of object
                pdDist[X]       = dRadius*cos(dBeta - dPhi);
                
                pnCurIdDct[X]   = pos2id(pdDist[X], pdDct[X], pnDct[X]) - pdOffsetDct[X];
                
                out                             = interpolation1d(&pin[pnDct[X]*iView], pnCurIdDct[X], pnDct[X]);
                pout[pnImg[Y]*iImgX + iImgY]	+= PI/nView * out;
            }
        }
    }
    
    return ;
}