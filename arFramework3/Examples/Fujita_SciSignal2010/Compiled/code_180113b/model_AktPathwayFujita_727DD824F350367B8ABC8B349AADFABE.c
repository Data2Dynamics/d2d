#include "model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE.h"
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





 void fu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(void *user_data, double t)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  data->u[0] = 3.0E-1;
  data->u[3] = 6.0E1;
  data->u[4] = 3.6E3;

  return;
}


 void fsu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(void *user_data, double t)
{

  return;
}


 void fv_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->v[0] = p[4]*x_tmp[0]*3.0E-1-p[5]*x_tmp[8];
  data->v[1] = -p[7]*x_tmp[2]+p[6]*x_tmp[1]*x_tmp[3];
  data->v[2] = p[8]*x_tmp[2];
  data->v[3] = p[9]*x_tmp[1];
  data->v[4] = -p[11]*x_tmp[6]+p[10]*x_tmp[4]*x_tmp[5];
  data->v[5] = p[12]*x_tmp[6];
  data->v[6] = p[13]*x_tmp[4];
  data->v[7] = p[14]*x_tmp[7];
  data->v[8] = p[15]*x_tmp[8];
  data->v[9] = p[0]*x_tmp[0];
  data->v[10] = p[0]*6.819E4;

  return;
}


 void dvdx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdx[0] = p[4]*3.0E-1;
  data->dvdx[9] = p[0];
  data->dvdx[12] = p[6]*x_tmp[3];
  data->dvdx[14] = p[9];
  data->dvdx[23] = -p[7];
  data->dvdx[24] = p[8];
  data->dvdx[34] = p[6]*x_tmp[1];
  data->dvdx[48] = p[10]*x_tmp[5];
  data->dvdx[50] = p[13];
  data->dvdx[59] = p[10]*x_tmp[4];
  data->dvdx[70] = -p[11];
  data->dvdx[71] = p[12];
  data->dvdx[84] = p[14];
  data->dvdx[88] = -p[5];
  data->dvdx[96] = p[15];

  return;
}


 void dvdu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);


  return;
}


 void dvdp_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdp[9] = x_tmp[0];
  data->dvdp[10] = 6.819E4;
  data->dvdp[44] = x_tmp[0]*3.0E-1;
  data->dvdp[55] = -x_tmp[8];
  data->dvdp[67] = x_tmp[1]*x_tmp[3];
  data->dvdp[78] = -x_tmp[2];
  data->dvdp[90] = x_tmp[2];
  data->dvdp[102] = x_tmp[1];
  data->dvdp[114] = x_tmp[4]*x_tmp[5];
  data->dvdp[125] = -x_tmp[6];
  data->dvdp[137] = x_tmp[6];
  data->dvdp[149] = x_tmp[4];
  data->dvdp[161] = x_tmp[7];
  data->dvdp[173] = x_tmp[8];

  return;
}


 int fx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  UserData data = (UserData) user_data;
  int is;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  fu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(data, t);
  fv_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
  xdot_tmp[0] = -v[0]-v[9]+v[10];
  xdot_tmp[1] = -v[1]+v[2]-v[3]+v[8];
  xdot_tmp[2] = v[1]-v[2];
  xdot_tmp[3] = -v[1]+v[6];
  xdot_tmp[4] = v[2]-v[4]+v[5]-v[6];
  xdot_tmp[5] = -v[4]+v[7];
  xdot_tmp[6] = v[4]-v[5];
  xdot_tmp[7] = v[5]-v[7];
  xdot_tmp[8] = v[0]-v[8];
  for (is=0; is<9; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
    if(qpositivex[is]>0.5 && x_tmp[is]<0.0 && xdot_tmp[is]<0.0) xdot_tmp[is] = -xdot_tmp[is];
  }

  return(*(data->abort));
}


 void fxdouble_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, double *xdot_tmp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  fu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(data, t);
  fv_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
  xdot_tmp[0] = -v[0]-v[9]+v[10];
  xdot_tmp[1] = -v[1]+v[2]-v[3]+v[8];
  xdot_tmp[2] = v[1]-v[2];
  xdot_tmp[3] = -v[1]+v[6];
  xdot_tmp[4] = v[2]-v[4]+v[5]-v[6];
  xdot_tmp[5] = -v[4]+v[7];
  xdot_tmp[6] = v[4]-v[5];
  xdot_tmp[7] = v[5]-v[7];
  xdot_tmp[8] = v[0]-v[8];
  for (is=0; is<9; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
  }

  return;
}


 void fx0_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  x0_tmp[0] = p[2];
  x0_tmp[3] = p[1];
  x0_tmp[5] = p[3];

  return;
}


 int dfxdx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(long int N, realtype t, N_Vector x, 
  	N_Vector fx, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
  for (is=0; is<81; is++) {
    J->data[is] = 0.0;
  }
  J->data[0] = -dvdx[0]-dvdx[9];
  J->data[8] = dvdx[0];
  J->data[10] = -dvdx[12]-dvdx[14];
  J->data[11] = dvdx[12];
  J->data[12] = -dvdx[12];
  J->data[19] = -dvdx[23]+dvdx[24];
  J->data[20] = dvdx[23]-dvdx[24];
  J->data[21] = -dvdx[23];
  J->data[22] = dvdx[24];
  J->data[28] = -dvdx[34];
  J->data[29] = dvdx[34];
  J->data[30] = -dvdx[34];
  J->data[39] = dvdx[50];
  J->data[40] = -dvdx[48]-dvdx[50];
  J->data[41] = -dvdx[48];
  J->data[42] = dvdx[48];
  J->data[49] = -dvdx[59];
  J->data[50] = -dvdx[59];
  J->data[51] = dvdx[59];
  J->data[58] = -dvdx[70]+dvdx[71];
  J->data[59] = -dvdx[70];
  J->data[60] = dvdx[70]-dvdx[71];
  J->data[61] = dvdx[71];
  J->data[68] = dvdx[84];
  J->data[70] = -dvdx[84];
  J->data[72] = -dvdx[88];
  J->data[73] = dvdx[96];
  J->data[80] = dvdx[88]-dvdx[96];
  for (is=0; is<81; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = 0.0;
  }

  return(0);
}


 int dfxdx_out_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, realtype* J, void *user_data) {
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
  for (is=0; is<81; is++) {
    J[is] = 0.0;
  }
  J[0] = -dvdx[0]-dvdx[9];
  J[8] = dvdx[0];
  J[10] = -dvdx[12]-dvdx[14];
  J[11] = dvdx[12];
  J[12] = -dvdx[12];
  J[19] = -dvdx[23]+dvdx[24];
  J[20] = dvdx[23]-dvdx[24];
  J[21] = -dvdx[23];
  J[22] = dvdx[24];
  J[28] = -dvdx[34];
  J[29] = dvdx[34];
  J[30] = -dvdx[34];
  J[39] = dvdx[50];
  J[40] = -dvdx[48]-dvdx[50];
  J[41] = -dvdx[48];
  J[42] = dvdx[48];
  J[49] = -dvdx[59];
  J[50] = -dvdx[59];
  J[51] = dvdx[59];
  J[58] = -dvdx[70]+dvdx[71];
  J[59] = -dvdx[70];
  J[60] = dvdx[70]-dvdx[71];
  J[61] = dvdx[71];
  J[68] = dvdx[84];
  J[70] = -dvdx[84];
  J[72] = -dvdx[88];
  J[73] = dvdx[96];
  J[80] = dvdx[88]-dvdx[96];
  for (is=0; is<81; is++) {
    if(mxIsNaN(J[is])) J[is] = 0.0;
  }

  return(0);
}


 int dfxdx_sparse_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, 
  	N_Vector fx, SlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
  SlsSetToZero(J);
    J->rowvals[0] = 0;
    J->rowvals[1] = 8;
    J->rowvals[2] = 1;
    J->rowvals[3] = 2;
    J->rowvals[4] = 3;
    J->rowvals[5] = 1;
    J->rowvals[6] = 2;
    J->rowvals[7] = 3;
    J->rowvals[8] = 4;
    J->rowvals[9] = 1;
    J->rowvals[10] = 2;
    J->rowvals[11] = 3;
    J->rowvals[12] = 3;
    J->rowvals[13] = 4;
    J->rowvals[14] = 5;
    J->rowvals[15] = 6;
    J->rowvals[16] = 4;
    J->rowvals[17] = 5;
    J->rowvals[18] = 6;
    J->rowvals[19] = 4;
    J->rowvals[20] = 5;
    J->rowvals[21] = 6;
    J->rowvals[22] = 7;
    J->rowvals[23] = 5;
    J->rowvals[24] = 7;
    J->rowvals[25] = 0;
    J->rowvals[26] = 1;
    J->rowvals[27] = 8;

    J->colptrs[0] = 0;
    J->colptrs[1] = 2;
    J->colptrs[2] = 5;
    J->colptrs[3] = 9;
    J->colptrs[4] = 12;
    J->colptrs[5] = 16;
    J->colptrs[6] = 19;
    J->colptrs[7] = 23;
    J->colptrs[8] = 25;
    J->colptrs[9] = 28;

  J->data[0] = -dvdx[0]-dvdx[9];
  J->data[1] = dvdx[0];
  J->data[2] = -dvdx[12]-dvdx[14];
  J->data[3] = dvdx[12];
  J->data[4] = -dvdx[12];
  J->data[5] = -dvdx[23]+dvdx[24];
  J->data[6] = dvdx[23]-dvdx[24];
  J->data[7] = -dvdx[23];
  J->data[8] = dvdx[24];
  J->data[9] = -dvdx[34];
  J->data[10] = dvdx[34];
  J->data[11] = -dvdx[34];
  J->data[12] = dvdx[50];
  J->data[13] = -dvdx[48]-dvdx[50];
  J->data[14] = -dvdx[48];
  J->data[15] = dvdx[48];
  J->data[16] = -dvdx[59];
  J->data[17] = -dvdx[59];
  J->data[18] = dvdx[59];
  J->data[19] = -dvdx[70]+dvdx[71];
  J->data[20] = -dvdx[70];
  J->data[21] = dvdx[70]-dvdx[71];
  J->data[22] = dvdx[71];
  J->data[23] = dvdx[84];
  J->data[24] = -dvdx[84];
  J->data[25] = -dvdx[88];
  J->data[26] = dvdx[96];
  J->data[27] = dvdx[88]-dvdx[96];
  for (is=0; is<28; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = RCONST(0.0);
  }

  return(0);
}


 int fsx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(int Ns, realtype t, N_Vector x, N_Vector xdot, 
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
  fsu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(data, t);
  dvdx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
  dvdu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
  dvdp_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
  for (is=0; is<11; is++) {
    sv[is] = 0.0;
  }
  sv[0] = dvdx[0]*sx_tmp[0]+dvdx[88]*sx_tmp[8];
  sv[1] = dvdx[12]*sx_tmp[1]+dvdx[23]*sx_tmp[2]+dvdx[34]*sx_tmp[3];
  sv[2] = dvdx[24]*sx_tmp[2];
  sv[3] = dvdx[14]*sx_tmp[1];
  sv[4] = dvdx[48]*sx_tmp[4]+dvdx[59]*sx_tmp[5]+dvdx[70]*sx_tmp[6];
  sv[5] = dvdx[71]*sx_tmp[6];
  sv[6] = dvdx[50]*sx_tmp[4];
  sv[7] = dvdx[84]*sx_tmp[7];
  sv[8] = dvdx[96]*sx_tmp[8];
  sv[9] = dvdx[9]*sx_tmp[0];
  switch (ip) {
    case 0: {
      sv[9] += dvdp[9];
      sv[10] += dvdp[10];
    } break;
    case 1: {

    } break;
    case 2: {

    } break;
    case 3: {

    } break;
    case 4: {
      sv[0] += dvdp[44];
    } break;
    case 5: {
      sv[0] += dvdp[55];
    } break;
    case 6: {
      sv[1] += dvdp[67];
    } break;
    case 7: {
      sv[1] += dvdp[78];
    } break;
    case 8: {
      sv[2] += dvdp[90];
    } break;
    case 9: {
      sv[3] += dvdp[102];
    } break;
    case 10: {
      sv[4] += dvdp[114];
    } break;
    case 11: {
      sv[4] += dvdp[125];
    } break;
    case 12: {
      sv[5] += dvdp[137];
    } break;
    case 13: {
      sv[6] += dvdp[149];
    } break;
    case 14: {
      sv[7] += dvdp[161];
    } break;
    case 15: {
      sv[8] += dvdp[173];
    } break;
    case 16: {

    } break;
    case 17: {

    } break;
    case 18: {

    } break;
  }
  sxdot_tmp[0] = -sv[0]-sv[9]+sv[10];
  sxdot_tmp[1] = -sv[1]+sv[2]-sv[3]+sv[8];
  sxdot_tmp[2] = sv[1]-sv[2];
  sxdot_tmp[3] = -sv[1]+sv[6];
  sxdot_tmp[4] = sv[2]-sv[4]+sv[5]-sv[6];
  sxdot_tmp[5] = -sv[4]+sv[7];
  sxdot_tmp[6] = sv[4]-sv[5];
  sxdot_tmp[7] = sv[5]-sv[7];
  sxdot_tmp[8] = sv[0]-sv[8];
  for (is=0; is<9; is++) {
    if(mxIsNaN(sxdot_tmp[is])) sxdot_tmp[is] = 0.0;
  }

  return(0);
}


 int subfsx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(int Ns, realtype t, N_Vector x, N_Vector xdot, 
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
 {
   UserData data = (UserData) user_data;
   return fsx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(Ns, t, x, xdot, data->sensIndices[ip], sx, sxdot, user_data, tmp1, tmp2);
 };

 void csv_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data)
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
  for (is=0; is<11; is++) {
    sv[is] = 0.0;
  }
  sv[0] = dvdx[0]*sx_tmp[0]+dvdx[88]*sx_tmp[8];
  sv[1] = dvdx[12]*sx_tmp[1]+dvdx[23]*sx_tmp[2]+dvdx[34]*sx_tmp[3];
  sv[2] = dvdx[24]*sx_tmp[2];
  sv[3] = dvdx[14]*sx_tmp[1];
  sv[4] = dvdx[48]*sx_tmp[4]+dvdx[59]*sx_tmp[5]+dvdx[70]*sx_tmp[6];
  sv[5] = dvdx[71]*sx_tmp[6];
  sv[6] = dvdx[50]*sx_tmp[4];
  sv[7] = dvdx[84]*sx_tmp[7];
  sv[8] = dvdx[96]*sx_tmp[8];
  sv[9] = dvdx[9]*sx_tmp[0];
  switch (ip) {
    case 0: {
      sv[9] += dvdp[9];
      sv[10] += dvdp[10];
    } break;
    case 1: {

    } break;
    case 2: {

    } break;
    case 3: {

    } break;
    case 4: {
      sv[0] += dvdp[44];
    } break;
    case 5: {
      sv[0] += dvdp[55];
    } break;
    case 6: {
      sv[1] += dvdp[67];
    } break;
    case 7: {
      sv[1] += dvdp[78];
    } break;
    case 8: {
      sv[2] += dvdp[90];
    } break;
    case 9: {
      sv[3] += dvdp[102];
    } break;
    case 10: {
      sv[4] += dvdp[114];
    } break;
    case 11: {
      sv[4] += dvdp[125];
    } break;
    case 12: {
      sv[5] += dvdp[137];
    } break;
    case 13: {
      sv[6] += dvdp[149];
    } break;
    case 14: {
      sv[7] += dvdp[161];
    } break;
    case 15: {
      sv[8] += dvdp[173];
    } break;
    case 16: {

    } break;
    case 17: {

    } break;
    case 18: {

    } break;
  }

  return;
}


 void fsx0_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  switch (ip) {
    case 1: {
      sx0_tmp[3] = 1.0;
    } break;
    case 2: {
      sx0_tmp[0] = 1.0;
    } break;
    case 3: {
      sx0_tmp[5] = 1.0;
    } break;
  }

  return;
}


 void subfsx0_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(int ip, N_Vector sx0, void *user_data)
 {
   UserData data = (UserData) user_data;
   fsx0_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(data->sensIndices[ip], sx0, user_data);
 };

 void dfxdp0_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, double *dfxdp0, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdp = data->dvdp;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *x_tmp = N_VGetArrayPointer(x);
  dfxdp0[0] = -dvdp[9]+dvdp[10];
  dfxdp0[10] = -dvdx[34];
  dfxdp0[11] = dvdx[34];
  dfxdp0[12] = -dvdx[34];
  dfxdp0[18] = -dvdx[0]-dvdx[9];
  dfxdp0[26] = dvdx[0];
  dfxdp0[31] = -dvdx[59];
  dfxdp0[32] = -dvdx[59];
  dfxdp0[33] = dvdx[59];
  dfxdp0[36] = -dvdp[44];
  dfxdp0[44] = dvdp[44];
  dfxdp0[45] = -dvdp[55];
  dfxdp0[53] = dvdp[55];
  dfxdp0[55] = -dvdp[67];
  dfxdp0[56] = dvdp[67];
  dfxdp0[57] = -dvdp[67];
  dfxdp0[64] = -dvdp[78];
  dfxdp0[65] = dvdp[78];
  dfxdp0[66] = -dvdp[78];
  dfxdp0[73] = dvdp[90];
  dfxdp0[74] = -dvdp[90];
  dfxdp0[76] = dvdp[90];
  dfxdp0[82] = -dvdp[102];
  dfxdp0[94] = -dvdp[114];
  dfxdp0[95] = -dvdp[114];
  dfxdp0[96] = dvdp[114];
  dfxdp0[103] = -dvdp[125];
  dfxdp0[104] = -dvdp[125];
  dfxdp0[105] = dvdp[125];
  dfxdp0[112] = dvdp[137];
  dfxdp0[114] = -dvdp[137];
  dfxdp0[115] = dvdp[137];
  dfxdp0[120] = dvdp[149];
  dfxdp0[121] = -dvdp[149];
  dfxdp0[131] = dvdp[161];
  dfxdp0[133] = -dvdp[161];
  dfxdp0[136] = dvdp[173];
  dfxdp0[143] = -dvdp[173];
  for (is=0; is<171; is++) {
    if(mxIsNaN(dfxdp0[is])) dfxdp0[is] = 0.0;
  }

  return;
}


 void dfxdp_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(realtype t, N_Vector x, double *dfxdp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdp = data->dvdp;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *x_tmp = N_VGetArrayPointer(x);
  dfxdp[0] = -dvdp[9]+dvdp[10];
  dfxdp[36] = -dvdp[44];
  dfxdp[44] = dvdp[44];
  dfxdp[45] = -dvdp[55];
  dfxdp[53] = dvdp[55];
  dfxdp[55] = -dvdp[67];
  dfxdp[56] = dvdp[67];
  dfxdp[57] = -dvdp[67];
  dfxdp[64] = -dvdp[78];
  dfxdp[65] = dvdp[78];
  dfxdp[66] = -dvdp[78];
  dfxdp[73] = dvdp[90];
  dfxdp[74] = -dvdp[90];
  dfxdp[76] = dvdp[90];
  dfxdp[82] = -dvdp[102];
  dfxdp[94] = -dvdp[114];
  dfxdp[95] = -dvdp[114];
  dfxdp[96] = dvdp[114];
  dfxdp[103] = -dvdp[125];
  dfxdp[104] = -dvdp[125];
  dfxdp[105] = dvdp[125];
  dfxdp[112] = dvdp[137];
  dfxdp[114] = -dvdp[137];
  dfxdp[115] = dvdp[137];
  dfxdp[120] = dvdp[149];
  dfxdp[121] = -dvdp[149];
  dfxdp[131] = dvdp[161];
  dfxdp[133] = -dvdp[161];
  dfxdp[136] = dvdp[173];
  dfxdp[143] = -dvdp[173];
  for (is=0; is<171; is++) {
    if(mxIsNaN(dfxdp[is])) dfxdp[is] = 0.0;
  }

  return;
}


 void fz_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x){
  z[nz*nt*iruns+it+nt*0] = 6.819E4;
  z[nz*nt*iruns+it+nt*1] = 3.0E-1;
  z[nz*nt*iruns+it+nt*2] = (x[nx*nt*iruns+it+nt*1]+x[nx*nt*iruns+it+nt*2])*p[17];
  z[nz*nt*iruns+it+nt*3] = (x[nx*nt*iruns+it+nt*4]+x[nx*nt*iruns+it+nt*6])*p[16];
  z[nz*nt*iruns+it+nt*4] = p[18]*x[nx*nt*iruns+it+nt*7];

  return;
}


 void fsz_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx){
  int jp;
  for (jp=0; jp<np; jp++) {
      sz[it + nt*5*jp + nt*2] = p[17]*sx[it + nt*9*jp + nt*1]+p[17]*sx[it + nt*9*jp + nt*2];
      sz[it + nt*5*jp + nt*3] = p[16]*sx[it + nt*9*jp + nt*4]+p[16]*sx[it + nt*9*jp + nt*6];
      sz[it + nt*5*jp + nt*4] = p[18]*sx[it + nt*9*jp + nt*7];
      sz[it + nt*5*jp + nt*0] = 0.0;
      sz[it + nt*5*jp + nt*1] = 0.0;
  };

  sz[it+nt*83] += x[it+nt*4]+x[it+nt*6];
  sz[it+nt*87] += x[it+nt*1]+x[it+nt*2];
  sz[it+nt*94] += x[it+nt*7];

  return;
}


 void dfzdx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x){
      dfzdxs[nx*nt*iruns+it+nt*7] = p[17];
      dfzdxs[nx*nt*iruns+it+nt*12] = p[17];
      dfzdxs[nx*nt*iruns+it+nt*23] = p[16];
      dfzdxs[nx*nt*iruns+it+nt*33] = p[16];
      dfzdxs[nx*nt*iruns+it+nt*39] = p[18];


  return;
}


