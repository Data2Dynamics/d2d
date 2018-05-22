#ifndef _MY_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A
#define _MY_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_sparse.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <cvodes/cvodes_klu.h>
#include <sundials/sundials_sparse.h>
#include <udata.h>
#include <math.h>
#include <mex.h>
#include <arInputFunctionsC.h>



 void fu_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(void *user_data, double t);
 void fsu_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(void *user_data, double t);
 void fv_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x, void *user_data);
 void dvdx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x, void *user_data);
 void dvdu_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x, void *user_data);
 void dvdp_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x, void *user_data);
 int fx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x, N_Vector xdot, void *user_data);
 void fxdouble_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x, double *xdot_tmp, void *user_data);
 void fx0_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(N_Vector x0, void *user_data);
 int dfxdx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(long int N, realtype t, N_Vector x,N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 int dfxdx_out_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x, realtype* J, void *user_data); int dfxdx_sparse_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x,N_Vector fx, SlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 int fsx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip, N_Vector sx, N_Vector sxdot, void *user_data,N_Vector tmp1, N_Vector tmp2);
 int subfsx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip, N_Vector sx, N_Vector sxdot, void *user_data,N_Vector tmp1, N_Vector tmp2);
 void fsx0_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(int ip, N_Vector sx0, void *user_data);
 void subfsx0_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(int ip, N_Vector sx0, void *user_data);
 void csv_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data);
 void dfxdp0_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x, double *dfxdp0, void *user_data);

 void dfxdp_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(realtype t, N_Vector x, double *dfxdp, void *user_data);

 void fz_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x);
 void fsz_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx);

 void dfzdx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x);
#endif /* _MY_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A */



