#ifndef _MY_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400
#define _MY_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400

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



 void fy_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z);
 void fystd_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z);
 void fsy_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz);
 void fsystd_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz);

 void fy_scale_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx);
#endif /* _MY_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400 */



