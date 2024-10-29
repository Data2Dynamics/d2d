#include "RTF_DoseDependent_F298A94C5A6864FEB754E83604898892.h"
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





 void fu_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(void *user_data, double t)
{

  return;
}


 void fsu_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(void *user_data, double t)
{

  return;
}


 void fv_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);


  return;
}


 void dvdx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);


  return;
}


 void dvdu_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, void *user_data)
{

  return;
}


 void dvdp_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);

  return;
}


 int fx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  UserData data = (UserData) user_data;
  int is;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  fu_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(data, t);
  fv_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(t, x, data);

  xdot_tmp[0] = 0.0;
  for (is=0; is<1; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
    if(qpositivex[is]>0.5 && x_tmp[is]<0.0 && xdot_tmp[is]<0.0) xdot_tmp[is] = -xdot_tmp[is];
  }

  return(*(data->abort));
}


 void fxdouble_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, double *xdot_tmp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  fu_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(data, t);
  fv_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(t, x, data);

  xdot_tmp[0] = 0.0;
  for (is=0; is<1; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
  }

  return;
}


 void fx0_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);


  return;
}


 int dfxdx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(long int N, realtype t, N_Vector x, 
  	N_Vector fx, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(t, x, data);
  for (is=0; is<1; is++) {
    J->data[is] = 0.0;
  }

  for (is=0; is<1; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = 0.0;
  }

  return(0);
}


 int dfxdx_out_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, realtype* J, void *user_data) {
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(t, x, data);
  for (is=0; is<1; is++) {
    J[is] = 0.0;
  }

  for (is=0; is<1; is++) {
    if(mxIsNaN(J[is])) J[is] = 0.0;
  }

  return(0);
}


 int dfxdx_sparse_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, 
  	N_Vector fx, SlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(t, x, data);
  SlsSetToZero(J);

    J->colptrs[0] = 0;


  for (is=0; is<0; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = RCONST(0.0);
  }

  return(0);
}


 int fsx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(int Ns, realtype t, N_Vector x, N_Vector xdot, 
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
  fsu_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(data, t);
  dvdx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(t, x, data);
  dvdu_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(t, x, data);
  dvdp_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(t, x, data);
  for (is=0; is<1; is++) {
    sv[is] = 0.0;
  }

  switch (ip) {
  }

  sxdot_tmp[0] = 0.0;
  for (is=0; is<1; is++) {
    if(mxIsNaN(sxdot_tmp[is])) sxdot_tmp[is] = 0.0;
  }

  return(0);
}


 int subfsx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(int Ns, realtype t, N_Vector x, N_Vector xdot, 
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
 {
   UserData data = (UserData) user_data;
   return fsx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(Ns, t, x, xdot, data->sensIndices[ip], sx, sxdot, user_data, tmp1, tmp2);
 };

 void csv_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data)
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
  for (is=0; is<1; is++) {
    sv[is] = 0.0;
  }

  switch (ip) {
  }

  return;
}


 #pragma optimize("", off)
 void fsx0_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  switch (ip) {
  }

  return;
}


 void subfsx0_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(int ip, N_Vector sx0, void *user_data)
 {
   UserData data = (UserData) user_data;
   fsx0_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(data->sensIndices[ip], sx0, user_data);
 };

 #pragma optimize("", on)
 void dfxdp0_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, double *dfxdp0, void *user_data)
{

  return;
}


 void dfxdp_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(realtype t, N_Vector x, double *dfxdp, void *user_data)
{

  return;
}


 void fz_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x){

  return;
}


 void fsz_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx){

  return;
}


 void dfzdx_RTF_DoseDependent_F298A94C5A6864FEB754E83604898892(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x){

  return;
}


