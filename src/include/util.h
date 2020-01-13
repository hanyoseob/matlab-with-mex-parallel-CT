#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;

/*
 *
 */
#define PI          3.14159265359f
#define Y           0
#define X           1
#define Z           2

#define REAL        0
#define IMAG        1

#define min(x, y)   (x < y ? x : y)
#define max(x, y)   (x > y ? x : y)

/*
 *
 */
float   id2pos(float id, float d, int n);
float   pos2id(float pos, float d, int n);
void    rot2d(float *rot, float *pos, float theta);

void	generate_filter(float *pker, float d, int n);
void	convolution1d(float *pout, float *pin, float *pker, int n);

float   interpolation1d(float *pin, float curid, int n);
float   interpolation2d(float *pin, float *pcurid, int *pn);


/*
 *
 */
float id2pos(float id, float d, int n) {
    float pos = (id - (n - 1.0f)/2.0f)*d;
    return pos ;
}

float pos2id(float pos, float d, int n) {
    float id = pos/d + (n - 1.0f)/2.0f;
    return id ;
}

void rot2d(float *rot, float *pos, float rad) {
    float COS_	= cos(rad);
    float SIN_	= sin(rad);

    rot[X]	= COS_*pos[X] - SIN_*pos[Y];
    rot[Y]	= SIN_*pos[X] + COS_*pos[Y];
    
    return ;
}


/*
 *
 */
void generate_filter(float *pker, float d, int n) {
    for (int i = -(n-1); i <= (n-1); i++) {
        if (i == 0) {
            pker[(n-1) + i] = 1.0f/(4.0f*d*d);
        } else {
            if (i%2) {
                pker[(n-1) + i] = -1.0f/((i*PI*d)*(i*PI*d));
            } else  {
                pker[(n-1) + i] = 0;
            }
        }
        
        pker[(n-1) + i] *= d;
    }
    return ;
}

void convolution1d(float *pout, float *pin, float *pker, int n) {
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < n; j++) 
            pout[i]    += pker[n-1 + i - j]*pin[j];
    
    return ;
}

/*
 *
 */
float interpolation1d(float *pin, float curid, int n) {
    
    if (curid < 0 || curid >= n-1) {
        return 0;
    }
    
    int preid  	= int(curid);
    int postid 	= preid + 1;
    
    float out   = (postid - curid)*pin[preid] + (curid - preid)*pin[postid];
    return out ;
}

float interpolation2d(float *pin, float *pcurid, int *pn) {

    float curidy    = pcurid[Y];
    float curidx    = pcurid[X];

    int ny          = pn[Y];
    int nx          = pn[X];
    
    if (curidx < 0 || curidy < 0 || curidx >= nx-1 || curidy >= ny-1) {
        return 0;
    }
    
    int preidy      = int(curidy);
    int postidy     = preidy + 1;
    
    int preidx      = int(curidx);
    int postidx     = preidx + 1;
            
    float prewgty   = postidy - curidy;
    float postwgty	= curidy - preidy;
    
    float prewgtx   = postidx - curidx;
    float postwgtx	= curidx - preidx;

    float out = prewgtx*(prewgty*pin[preidy + ny*preidx] + postwgty*pin[postidy + ny*preidx]) + 
                postwgtx*(prewgty*pin[preidy + ny*postidx] + postwgty*pin[postidy + ny*postidx]);
    return out;
}