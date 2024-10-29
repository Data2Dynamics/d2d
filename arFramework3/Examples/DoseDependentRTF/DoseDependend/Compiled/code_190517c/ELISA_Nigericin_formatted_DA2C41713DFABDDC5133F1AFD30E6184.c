#include "ELISA_Nigericin_formatted_DA2C41713DFABDDC5133F1AFD30E6184.h"
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





 void fy_ELISA_Nigericin_formatted_DA2C41713DFABDDC5133F1AFD30E6184(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z){
  y[ny*nt*iruns+it+nt*0] = p[0];

  return;
}


 void fystd_ELISA_Nigericin_formatted_DA2C41713DFABDDC5133F1AFD30E6184(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z){
  ystd[it+nt*0] = p[1];

  return;
}


 void fsy_ELISA_Nigericin_formatted_DA2C41713DFABDDC5133F1AFD30E6184(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz){
  sy[it+nt*0] = 1.0;
  sy[it+nt*1] = 0.0;

  return;
}


 void fsystd_ELISA_Nigericin_formatted_DA2C41713DFABDDC5133F1AFD30E6184(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz){
  systd[it+nt*0] = 0.0;
  systd[it+nt*1] = 1.0;

  return;
}


 void fy_scale_ELISA_Nigericin_formatted_DA2C41713DFABDDC5133F1AFD30E6184(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx){


  return;
}


