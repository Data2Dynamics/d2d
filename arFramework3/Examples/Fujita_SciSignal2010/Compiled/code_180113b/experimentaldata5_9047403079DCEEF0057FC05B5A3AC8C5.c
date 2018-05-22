#include "experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5.h"
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_sparse.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_sparse.h>
#include <cvodes/cvodes_klu.h>
#include <udata.h>
#include <math.h>
#include <mex.h>
#include <arInputFunctionsC.h>





 void fy_experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z){
  y[ny*nt*iruns+it+nt*0] = z[nz*ntlink*iruns+itlink+ntlink*2];
  y[ny*nt*iruns+it+nt*1] = z[nz*ntlink*iruns+itlink+ntlink*3];
  y[ny*nt*iruns+it+nt*2] = z[nz*ntlink*iruns+itlink+ntlink*4];

  return;
}


 void fystd_experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z){
  ystd[it+nt*0] = p[20];
  ystd[it+nt*1] = p[19];
  ystd[it+nt*2] = p[21];

  return;
}


 void fsy_experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz){
  sy[it+nt*0] = sz[itlink+ntlink*2];
  sy[it+nt*1] = sz[itlink+ntlink*3];
  sy[it+nt*2] = sz[itlink+ntlink*4];
  sy[it+nt*3] = sz[itlink+ntlink*7];
  sy[it+nt*4] = sz[itlink+ntlink*8];
  sy[it+nt*5] = sz[itlink+ntlink*9];
  sy[it+nt*6] = sz[itlink+ntlink*12];
  sy[it+nt*7] = sz[itlink+ntlink*13];
  sy[it+nt*8] = sz[itlink+ntlink*14];
  sy[it+nt*9] = sz[itlink+ntlink*17];
  sy[it+nt*10] = sz[itlink+ntlink*18];
  sy[it+nt*11] = sz[itlink+ntlink*19];
  sy[it+nt*12] = sz[itlink+ntlink*22];
  sy[it+nt*13] = sz[itlink+ntlink*23];
  sy[it+nt*14] = sz[itlink+ntlink*24];
  sy[it+nt*15] = sz[itlink+ntlink*27];
  sy[it+nt*16] = sz[itlink+ntlink*28];
  sy[it+nt*17] = sz[itlink+ntlink*29];
  sy[it+nt*18] = sz[itlink+ntlink*32];
  sy[it+nt*19] = sz[itlink+ntlink*33];
  sy[it+nt*20] = sz[itlink+ntlink*34];
  sy[it+nt*21] = sz[itlink+ntlink*37];
  sy[it+nt*22] = sz[itlink+ntlink*38];
  sy[it+nt*23] = sz[itlink+ntlink*39];
  sy[it+nt*24] = sz[itlink+ntlink*42];
  sy[it+nt*25] = sz[itlink+ntlink*43];
  sy[it+nt*26] = sz[itlink+ntlink*44];
  sy[it+nt*27] = sz[itlink+ntlink*47];
  sy[it+nt*28] = sz[itlink+ntlink*48];
  sy[it+nt*29] = sz[itlink+ntlink*49];
  sy[it+nt*30] = sz[itlink+ntlink*52];
  sy[it+nt*31] = sz[itlink+ntlink*53];
  sy[it+nt*32] = sz[itlink+ntlink*54];
  sy[it+nt*33] = sz[itlink+ntlink*57];
  sy[it+nt*34] = sz[itlink+ntlink*58];
  sy[it+nt*35] = sz[itlink+ntlink*59];
  sy[it+nt*36] = sz[itlink+ntlink*62];
  sy[it+nt*37] = sz[itlink+ntlink*63];
  sy[it+nt*38] = sz[itlink+ntlink*64];
  sy[it+nt*39] = sz[itlink+ntlink*67];
  sy[it+nt*40] = sz[itlink+ntlink*68];
  sy[it+nt*41] = sz[itlink+ntlink*69];
  sy[it+nt*42] = sz[itlink+ntlink*72];
  sy[it+nt*43] = sz[itlink+ntlink*73];
  sy[it+nt*44] = sz[itlink+ntlink*74];
  sy[it+nt*45] = sz[itlink+ntlink*77];
  sy[it+nt*46] = sz[itlink+ntlink*78];
  sy[it+nt*47] = sz[itlink+ntlink*79];
  sy[it+nt*48] = sz[itlink+ntlink*82];
  sy[it+nt*49] = sz[itlink+ntlink*83];
  sy[it+nt*50] = sz[itlink+ntlink*84];
  sy[it+nt*51] = sz[itlink+ntlink*87];
  sy[it+nt*52] = sz[itlink+ntlink*88];
  sy[it+nt*53] = sz[itlink+ntlink*89];
  sy[it+nt*54] = sz[itlink+ntlink*92];
  sy[it+nt*55] = sz[itlink+ntlink*93];
  sy[it+nt*56] = sz[itlink+ntlink*94];

  return;
}


 void fsystd_experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz){
  systd[it+nt*58] = 1.0;
  systd[it+nt*60] = 1.0;
  systd[it+nt*65] = 1.0;

  return;
}


 void fy_scale_experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx){
  y_scale[ny*nt*iruns+it+nt*0] = dfzdx[nx*ntlink*iruns+itlink+ntlink*2];
  y_scale[ny*nt*iruns+it+nt*1] = dfzdx[nx*ntlink*iruns+itlink+ntlink*3];
  y_scale[ny*nt*iruns+it+nt*2] = dfzdx[nx*ntlink*iruns+itlink+ntlink*4];
  y_scale[ny*nt*iruns+it+nt*3] = dfzdx[nx*ntlink*iruns+itlink+ntlink*7];
  y_scale[ny*nt*iruns+it+nt*4] = dfzdx[nx*ntlink*iruns+itlink+ntlink*8];
  y_scale[ny*nt*iruns+it+nt*5] = dfzdx[nx*ntlink*iruns+itlink+ntlink*9];
  y_scale[ny*nt*iruns+it+nt*6] = dfzdx[nx*ntlink*iruns+itlink+ntlink*12];
  y_scale[ny*nt*iruns+it+nt*7] = dfzdx[nx*ntlink*iruns+itlink+ntlink*13];
  y_scale[ny*nt*iruns+it+nt*8] = dfzdx[nx*ntlink*iruns+itlink+ntlink*14];
  y_scale[ny*nt*iruns+it+nt*9] = dfzdx[nx*ntlink*iruns+itlink+ntlink*17];
  y_scale[ny*nt*iruns+it+nt*10] = dfzdx[nx*ntlink*iruns+itlink+ntlink*18];
  y_scale[ny*nt*iruns+it+nt*11] = dfzdx[nx*ntlink*iruns+itlink+ntlink*19];
  y_scale[ny*nt*iruns+it+nt*12] = dfzdx[nx*ntlink*iruns+itlink+ntlink*22];
  y_scale[ny*nt*iruns+it+nt*13] = dfzdx[nx*ntlink*iruns+itlink+ntlink*23];
  y_scale[ny*nt*iruns+it+nt*14] = dfzdx[nx*ntlink*iruns+itlink+ntlink*24];
  y_scale[ny*nt*iruns+it+nt*15] = dfzdx[nx*ntlink*iruns+itlink+ntlink*27];
  y_scale[ny*nt*iruns+it+nt*16] = dfzdx[nx*ntlink*iruns+itlink+ntlink*28];
  y_scale[ny*nt*iruns+it+nt*17] = dfzdx[nx*ntlink*iruns+itlink+ntlink*29];
  y_scale[ny*nt*iruns+it+nt*18] = dfzdx[nx*ntlink*iruns+itlink+ntlink*32];
  y_scale[ny*nt*iruns+it+nt*19] = dfzdx[nx*ntlink*iruns+itlink+ntlink*33];
  y_scale[ny*nt*iruns+it+nt*20] = dfzdx[nx*ntlink*iruns+itlink+ntlink*34];
  y_scale[ny*nt*iruns+it+nt*21] = dfzdx[nx*ntlink*iruns+itlink+ntlink*37];
  y_scale[ny*nt*iruns+it+nt*22] = dfzdx[nx*ntlink*iruns+itlink+ntlink*38];
  y_scale[ny*nt*iruns+it+nt*23] = dfzdx[nx*ntlink*iruns+itlink+ntlink*39];
  y_scale[ny*nt*iruns+it+nt*24] = dfzdx[nx*ntlink*iruns+itlink+ntlink*42];
  y_scale[ny*nt*iruns+it+nt*25] = dfzdx[nx*ntlink*iruns+itlink+ntlink*43];
  y_scale[ny*nt*iruns+it+nt*26] = dfzdx[nx*ntlink*iruns+itlink+ntlink*44];

  return;
}

