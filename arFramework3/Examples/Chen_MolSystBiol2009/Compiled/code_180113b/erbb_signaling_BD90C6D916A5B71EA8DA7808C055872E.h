#ifndef _MY_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E
#define _MY_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E

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



 void fu_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(void *user_data, double t);
 void fsu_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(void *user_data, double t);
 void fv_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x, void *user_data);
 void dvdx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x, void *user_data);
 void dvdu_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x, void *user_data);
 void dvdp_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x, void *user_data);
 int fx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x, N_Vector xdot, void *user_data);
 void fxdouble_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x, double *xdot_tmp, void *user_data);
 void fx0_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(N_Vector x0, void *user_data);
 int dfxdx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(long int N, realtype t, N_Vector x,N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 int dfxdx_out_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x, realtype* J, void *user_data); int dfxdx_sparse_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x,N_Vector fx, SlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 int fsx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip, N_Vector sx, N_Vector sxdot, void *user_data,N_Vector tmp1, N_Vector tmp2);
 int subfsx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip, N_Vector sx, N_Vector sxdot, void *user_data,N_Vector tmp1, N_Vector tmp2);
 void fsx0_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(int ip, N_Vector sx0, void *user_data);
 void subfsx0_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(int ip, N_Vector sx0, void *user_data);
 void csv_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data);
 void dfxdp0_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x, double *dfxdp0, void *user_data);

 void dfxdp_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(realtype t, N_Vector x, double *dfxdp, void *user_data);

 void fz_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x);
 void fsz_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx);

 void dfzdx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x);
#endif /* _MY_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E */



