#include "ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400.h"
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





 void fy_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z){
  y[ny*nt*iruns+it+nt*0] = x[nx*ntlink*iruns+itlink+ntlink*2];
  y[ny*nt*iruns+it+nt*1] = x[nx*ntlink*iruns+itlink+ntlink*1];
  y[ny*nt*iruns+it+nt*2] = x[nx*ntlink*iruns+itlink+ntlink*0];

  return;
}


 void fystd_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z){
  ystd[it+nt*0] = p[5];
  ystd[it+nt*1] = p[4];
  ystd[it+nt*2] = p[3];

  return;
}


 void fsy_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz){
  sy[it+nt*0] = sx[itlink+ntlink*2];
  sy[it+nt*1] = sx[itlink+ntlink*1];
  sy[it+nt*2] = sx[itlink+ntlink*0];
  sy[it+nt*3] = sx[itlink+ntlink*5];
  sy[it+nt*4] = sx[itlink+ntlink*4];
  sy[it+nt*5] = sx[itlink+ntlink*3];
  sy[it+nt*6] = sx[itlink+ntlink*8];
  sy[it+nt*7] = sx[itlink+ntlink*7];
  sy[it+nt*8] = sx[itlink+ntlink*6];

  return;
}


 void fsystd_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz){
  systd[it+nt*11] = 1.0;
  systd[it+nt*13] = 1.0;
  systd[it+nt*15] = 1.0;

  return;
}


 void fy_scale_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx){
  y_scale[ny*nt*iruns+it+nt*2] = 1.0;
  y_scale[ny*nt*iruns+it+nt*4] = 1.0;
  y_scale[ny*nt*iruns+it+nt*6] = 1.0;

  return;
}


