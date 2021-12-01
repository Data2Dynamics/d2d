#include "ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3.h"
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





 void fu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(void *user_data, double t)
{

  return;
}


 void fsu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(void *user_data, double t)
{

  return;
}


 void fv_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->v[0] = p[1]*x_tmp[0];
  data->v[1] = p[2]*x_tmp[1];

  return;
}


 void dvdx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdx[0] = p[1];
  data->dvdx[3] = p[2];

  return;
}


 void dvdu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, void *user_data)
{

  return;
}


 void dvdp_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdp[2] = x_tmp[0];
  data->dvdp[5] = x_tmp[1];

  return;
}


 int fx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  UserData data = (UserData) user_data;
  int is;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  fu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(data, t);
  fv_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
  xdot_tmp[0] = -v[0];
  xdot_tmp[1] = v[0]-v[1];
  xdot_tmp[2] = v[1];
  for (is=0; is<3; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
    if(qpositivex[is]>0.5 && x_tmp[is]<0.0 && xdot_tmp[is]<0.0) xdot_tmp[is] = -xdot_tmp[is];
  }

  return(*(data->abort));
}


 void fxdouble_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, double *xdot_tmp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  fu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(data, t);
  fv_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
  xdot_tmp[0] = -v[0];
  xdot_tmp[1] = v[0]-v[1];
  xdot_tmp[2] = v[1];
  for (is=0; is<3; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
  }

  return;
}


 void fx0_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  x0_tmp[0] = p[0];

  return;
}


 int dfxdx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(long int N, realtype t, N_Vector x, 
  	N_Vector fx, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
  for (is=0; is<9; is++) {
    J->data[is] = 0.0;
  }
  J->data[0] = -dvdx[0];
  J->data[1] = dvdx[0];
  J->data[4] = -dvdx[3];
  J->data[5] = dvdx[3];
  for (is=0; is<9; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = 0.0;
  }

  return(0);
}


 int dfxdx_out_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, realtype* J, void *user_data) {
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
  for (is=0; is<9; is++) {
    J[is] = 0.0;
  }
  J[0] = -dvdx[0];
  J[1] = dvdx[0];
  J[4] = -dvdx[3];
  J[5] = dvdx[3];
  for (is=0; is<9; is++) {
    if(mxIsNaN(J[is])) J[is] = 0.0;
  }

  return(0);
}


 int dfxdx_sparse_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, 
  	N_Vector fx, SlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
  SlsSetToZero(J);
    J->rowvals[0] = 0;
    J->rowvals[1] = 1;
    J->rowvals[2] = 1;
    J->rowvals[3] = 2;
    J->rowvals[4] = 2;

    J->colptrs[0] = 0;
    J->colptrs[1] = 2;
    J->colptrs[2] = 4;
    J->colptrs[3] = 5;

  J->data[0] = -dvdx[0];
  J->data[1] = dvdx[0];
  J->data[2] = -dvdx[3];
  J->data[3] = dvdx[3];
  J->data[4] = RCONST(0.0);
  for (is=0; is<5; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = RCONST(0.0);
  }

  return(0);
}


 int fsx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(int Ns, realtype t, N_Vector x, N_Vector xdot, 
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *sv = data->sv;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *dvdp = data->dvdp;
  double *x_tmp = N_VGetArrayPointer(x);
  double *su = data->su;
  double *sx_tmp = N_VGetArrayPointer(sx);
  double *sxdot_tmp = N_VGetArrayPointer(sxdot);
  fsu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(data, t);
  dvdx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
  dvdu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
  dvdp_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
  for (is=0; is<2; is++) {
    sv[is] = 0.0;
  }
  sv[0] = dvdx[0]*sx_tmp[0];
  sv[1] = dvdx[3]*sx_tmp[1];
  switch (ip) {
    case 0: {

    } break;
    case 1: {
      sv[0] += dvdp[2];
    } break;
    case 2: {
      sv[1] += dvdp[5];
    } break;
  }
  sxdot_tmp[0] = -sv[0];
  sxdot_tmp[1] = sv[0]-sv[1];
  sxdot_tmp[2] = sv[1];
  for (is=0; is<3; is++) {
    if(mxIsNaN(sxdot_tmp[is])) sxdot_tmp[is] = 0.0;
  }

  return(0);
}


 int subfsx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(int Ns, realtype t, N_Vector x, N_Vector xdot, 
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
 {
   UserData data = (UserData) user_data;
   return fsx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(Ns, t, x, xdot, data->sensIndices[ip], sx, sxdot, user_data, tmp1, tmp2);
 };

 void csv_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *sv = data->sv;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *dvdp = data->dvdp;
  double *x_tmp = N_VGetArrayPointer(x);
  double *su = data->su;
  double *sx_tmp = N_VGetArrayPointer(sx);
  for (is=0; is<2; is++) {
    sv[is] = 0.0;
  }
  sv[0] = dvdx[0]*sx_tmp[0];
  sv[1] = dvdx[3]*sx_tmp[1];
  switch (ip) {
    case 0: {

    } break;
    case 1: {
      sv[0] += dvdp[2];
    } break;
    case 2: {
      sv[1] += dvdp[5];
    } break;
  }

  return;
}


 #pragma optimize("", off)
 void fsx0_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  switch (ip) {
    case 0: {
      sx0_tmp[0] = 1.0;
    } break;
  }

  return;
}


 void subfsx0_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(int ip, N_Vector sx0, void *user_data)
 {
   UserData data = (UserData) user_data;
   fsx0_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(data->sensIndices[ip], sx0, user_data);
 };

 #pragma optimize("", on)
 void dfxdp0_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, double *dfxdp0, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdp = data->dvdp;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *x_tmp = N_VGetArrayPointer(x);
  dfxdp0[0] = -dvdx[0];
  dfxdp0[1] = dvdx[0];
  dfxdp0[3] = -dvdp[2];
  dfxdp0[4] = dvdp[2];
  dfxdp0[7] = -dvdp[5];
  dfxdp0[8] = dvdp[5];
  for (is=0; is<9; is++) {
    if(mxIsNaN(dfxdp0[is])) dfxdp0[is] = 0.0;
  }

  return;
}


 void dfxdp_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(realtype t, N_Vector x, double *dfxdp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdp = data->dvdp;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *x_tmp = N_VGetArrayPointer(x);
  dfxdp[3] = -dvdp[2];
  dfxdp[4] = dvdp[2];
  dfxdp[7] = -dvdp[5];
  dfxdp[8] = dvdp[5];
  for (is=0; is<9; is++) {
    if(mxIsNaN(dfxdp[is])) dfxdp[is] = 0.0;
  }

  return;
}


 void fz_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x){

  return;
}


 void fsz_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx){

  return;
}


 void dfzdx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x){

  return;
}


