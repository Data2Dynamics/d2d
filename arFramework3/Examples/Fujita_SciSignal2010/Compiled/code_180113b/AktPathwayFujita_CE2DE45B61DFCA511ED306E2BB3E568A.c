#include "AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A.h"
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





 void fu_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(void *user_data, double t)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  data->u[0] = 3.0E1;
  data->u[3] = 6.0E1;
  data->u[4] = 3.6E3;

  return;
}


 void fsu_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(void *user_data, double t)
{

  return;
}


 void fv_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->v[0] = p[9]*x_tmp[1]+(p[7]-x_tmp[0])*p[0]-p[8]*u[0]*x_tmp[0];
  data->v[1] = -p[9]*x_tmp[1]-p[10]*x_tmp[1]+p[8]*u[0]*x_tmp[0];
  data->v[2] = p[10]*x_tmp[1]-p[11]*x_tmp[2]+p[13]*x_tmp[4]+p[14]*x_tmp[4]-p[12]*x_tmp[2]*x_tmp[3];
  data->v[3] = p[13]*x_tmp[4]+p[15]*x_tmp[5]-p[12]*x_tmp[2]*x_tmp[3];
  data->v[4] = -p[13]*x_tmp[4]-p[14]*x_tmp[4]+p[12]*x_tmp[2]*x_tmp[3];
  data->v[5] = p[14]*x_tmp[4]-p[15]*x_tmp[5]+p[17]*x_tmp[7]+p[18]*x_tmp[7]-p[16]*x_tmp[5]*x_tmp[6];
  data->v[6] = p[17]*x_tmp[7]+p[19]*x_tmp[8]-p[16]*x_tmp[5]*x_tmp[6];
  data->v[7] = -p[17]*x_tmp[7]-p[18]*x_tmp[7]+p[16]*x_tmp[5]*x_tmp[6];
  data->v[8] = p[18]*x_tmp[7]-p[19]*x_tmp[8];

  return;
}


 void dvdx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdx[0] = -p[8]*u[0]-p[0];
  data->dvdx[1] = p[8]*u[0];
  data->dvdx[9] = p[9];
  data->dvdx[10] = -p[9]-p[10];
  data->dvdx[11] = p[10];
  data->dvdx[20] = -p[12]*x_tmp[3]-p[11];
  data->dvdx[21] = -p[12]*x_tmp[3];
  data->dvdx[22] = p[12]*x_tmp[3];
  data->dvdx[29] = -p[12]*x_tmp[2];
  data->dvdx[30] = -p[12]*x_tmp[2];
  data->dvdx[31] = p[12]*x_tmp[2];
  data->dvdx[38] = p[13]+p[14];
  data->dvdx[39] = p[13];
  data->dvdx[40] = -p[13]-p[14];
  data->dvdx[41] = p[14];
  data->dvdx[48] = p[15];
  data->dvdx[50] = -p[16]*x_tmp[6]-p[15];
  data->dvdx[51] = -p[16]*x_tmp[6];
  data->dvdx[52] = p[16]*x_tmp[6];
  data->dvdx[59] = -p[16]*x_tmp[5];
  data->dvdx[60] = -p[16]*x_tmp[5];
  data->dvdx[61] = p[16]*x_tmp[5];
  data->dvdx[68] = p[17]+p[18];
  data->dvdx[69] = p[17];
  data->dvdx[70] = -p[17]-p[18];
  data->dvdx[71] = p[18];
  data->dvdx[78] = p[19];
  data->dvdx[80] = -p[19];

  return;
}


 void dvdu_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdu[0] = -p[8]*x_tmp[0];
  data->dvdu[1] = p[8]*x_tmp[0];

  return;
}


 void dvdp_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdp[0] = p[7]-x_tmp[0];
  data->dvdp[63] = p[0];
  data->dvdp[72] = -u[0]*x_tmp[0];
  data->dvdp[73] = u[0]*x_tmp[0];
  data->dvdp[81] = x_tmp[1];
  data->dvdp[82] = -x_tmp[1];
  data->dvdp[91] = -x_tmp[1];
  data->dvdp[92] = x_tmp[1];
  data->dvdp[101] = -x_tmp[2];
  data->dvdp[110] = -x_tmp[2]*x_tmp[3];
  data->dvdp[111] = -x_tmp[2]*x_tmp[3];
  data->dvdp[112] = x_tmp[2]*x_tmp[3];
  data->dvdp[119] = x_tmp[4];
  data->dvdp[120] = x_tmp[4];
  data->dvdp[121] = -x_tmp[4];
  data->dvdp[128] = x_tmp[4];
  data->dvdp[130] = -x_tmp[4];
  data->dvdp[131] = x_tmp[4];
  data->dvdp[138] = x_tmp[5];
  data->dvdp[140] = -x_tmp[5];
  data->dvdp[149] = -x_tmp[5]*x_tmp[6];
  data->dvdp[150] = -x_tmp[5]*x_tmp[6];
  data->dvdp[151] = x_tmp[5]*x_tmp[6];
  data->dvdp[158] = x_tmp[7];
  data->dvdp[159] = x_tmp[7];
  data->dvdp[160] = -x_tmp[7];
  data->dvdp[167] = x_tmp[7];
  data->dvdp[169] = -x_tmp[7];
  data->dvdp[170] = x_tmp[7];
  data->dvdp[177] = x_tmp[8];
  data->dvdp[179] = -x_tmp[8];

  return;
}


 int fx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  UserData data = (UserData) user_data;
  int is;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  fu_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(data, t);
  fv_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(t, x, data);
  xdot_tmp[0] = v[0];
  xdot_tmp[1] = v[1];
  xdot_tmp[2] = v[2];
  xdot_tmp[3] = v[3];
  xdot_tmp[4] = v[4];
  xdot_tmp[5] = v[5];
  xdot_tmp[6] = v[6];
  xdot_tmp[7] = v[7];
  xdot_tmp[8] = v[8];
  for (is=0; is<9; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
    if(qpositivex[is]>0.5 && x_tmp[is]<0.0 && xdot_tmp[is]<0.0) xdot_tmp[is] = -xdot_tmp[is];
  }

  return(*(data->abort));
}


 void fxdouble_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, double *xdot_tmp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  fu_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(data, t);
  fv_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(t, x, data);
  xdot_tmp[0] = v[0];
  xdot_tmp[1] = v[1];
  xdot_tmp[2] = v[2];
  xdot_tmp[3] = v[3];
  xdot_tmp[4] = v[4];
  xdot_tmp[5] = v[5];
  xdot_tmp[6] = v[6];
  xdot_tmp[7] = v[7];
  xdot_tmp[8] = v[8];
  for (is=0; is<9; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
  }

  return;
}


 void fx0_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  x0_tmp[0] = p[2];
  x0_tmp[3] = p[1];
  x0_tmp[6] = p[3];

  return;
}


 int dfxdx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(long int N, realtype t, N_Vector x, 
  	N_Vector fx, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(t, x, data);
  for (is=0; is<81; is++) {
    J->data[is] = 0.0;
  }
  J->data[0] = dvdx[0];
  J->data[1] = dvdx[1];
  J->data[9] = dvdx[9];
  J->data[10] = dvdx[10];
  J->data[11] = dvdx[11];
  J->data[20] = dvdx[20];
  J->data[21] = dvdx[21];
  J->data[22] = dvdx[22];
  J->data[29] = dvdx[29];
  J->data[30] = dvdx[30];
  J->data[31] = dvdx[31];
  J->data[38] = dvdx[38];
  J->data[39] = dvdx[39];
  J->data[40] = dvdx[40];
  J->data[41] = dvdx[41];
  J->data[48] = dvdx[48];
  J->data[50] = dvdx[50];
  J->data[51] = dvdx[51];
  J->data[52] = dvdx[52];
  J->data[59] = dvdx[59];
  J->data[60] = dvdx[60];
  J->data[61] = dvdx[61];
  J->data[68] = dvdx[68];
  J->data[69] = dvdx[69];
  J->data[70] = dvdx[70];
  J->data[71] = dvdx[71];
  J->data[78] = dvdx[78];
  J->data[80] = dvdx[80];
  for (is=0; is<81; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = 0.0;
  }

  return(0);
}


 int dfxdx_out_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, realtype* J, void *user_data) {
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(t, x, data);
  for (is=0; is<81; is++) {
    J[is] = 0.0;
  }
  J[0] = dvdx[0];
  J[1] = dvdx[1];
  J[9] = dvdx[9];
  J[10] = dvdx[10];
  J[11] = dvdx[11];
  J[20] = dvdx[20];
  J[21] = dvdx[21];
  J[22] = dvdx[22];
  J[29] = dvdx[29];
  J[30] = dvdx[30];
  J[31] = dvdx[31];
  J[38] = dvdx[38];
  J[39] = dvdx[39];
  J[40] = dvdx[40];
  J[41] = dvdx[41];
  J[48] = dvdx[48];
  J[50] = dvdx[50];
  J[51] = dvdx[51];
  J[52] = dvdx[52];
  J[59] = dvdx[59];
  J[60] = dvdx[60];
  J[61] = dvdx[61];
  J[68] = dvdx[68];
  J[69] = dvdx[69];
  J[70] = dvdx[70];
  J[71] = dvdx[71];
  J[78] = dvdx[78];
  J[80] = dvdx[80];
  for (is=0; is<81; is++) {
    if(mxIsNaN(J[is])) J[is] = 0.0;
  }

  return(0);
}


 int dfxdx_sparse_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, 
  	N_Vector fx, SlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(t, x, data);
  SlsSetToZero(J);
    J->rowvals[0] = 0;
    J->rowvals[1] = 1;
    J->rowvals[2] = 0;
    J->rowvals[3] = 1;
    J->rowvals[4] = 2;
    J->rowvals[5] = 2;
    J->rowvals[6] = 3;
    J->rowvals[7] = 4;
    J->rowvals[8] = 2;
    J->rowvals[9] = 3;
    J->rowvals[10] = 4;
    J->rowvals[11] = 2;
    J->rowvals[12] = 3;
    J->rowvals[13] = 4;
    J->rowvals[14] = 5;
    J->rowvals[15] = 3;
    J->rowvals[16] = 5;
    J->rowvals[17] = 6;
    J->rowvals[18] = 7;
    J->rowvals[19] = 5;
    J->rowvals[20] = 6;
    J->rowvals[21] = 7;
    J->rowvals[22] = 5;
    J->rowvals[23] = 6;
    J->rowvals[24] = 7;
    J->rowvals[25] = 8;
    J->rowvals[26] = 6;
    J->rowvals[27] = 8;

    J->colptrs[0] = 0;
    J->colptrs[1] = 2;
    J->colptrs[2] = 5;
    J->colptrs[3] = 8;
    J->colptrs[4] = 11;
    J->colptrs[5] = 15;
    J->colptrs[6] = 19;
    J->colptrs[7] = 22;
    J->colptrs[8] = 26;
    J->colptrs[9] = 28;

  J->data[0] = dvdx[0];
  J->data[1] = dvdx[1];
  J->data[2] = dvdx[9];
  J->data[3] = dvdx[10];
  J->data[4] = dvdx[11];
  J->data[5] = dvdx[20];
  J->data[6] = dvdx[21];
  J->data[7] = dvdx[22];
  J->data[8] = dvdx[29];
  J->data[9] = dvdx[30];
  J->data[10] = dvdx[31];
  J->data[11] = dvdx[38];
  J->data[12] = dvdx[39];
  J->data[13] = dvdx[40];
  J->data[14] = dvdx[41];
  J->data[15] = dvdx[48];
  J->data[16] = dvdx[50];
  J->data[17] = dvdx[51];
  J->data[18] = dvdx[52];
  J->data[19] = dvdx[59];
  J->data[20] = dvdx[60];
  J->data[21] = dvdx[61];
  J->data[22] = dvdx[68];
  J->data[23] = dvdx[69];
  J->data[24] = dvdx[70];
  J->data[25] = dvdx[71];
  J->data[26] = dvdx[78];
  J->data[27] = dvdx[80];
  for (is=0; is<28; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = RCONST(0.0);
  }

  return(0);
}


 int fsx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(int Ns, realtype t, N_Vector x, N_Vector xdot, 
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
  fsu_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(data, t);
  dvdx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(t, x, data);
  dvdu_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(t, x, data);
  dvdp_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(t, x, data);
  for (is=0; is<9; is++) {
    sv[is] = 0.0;
  }
  sv[0] = dvdu[0]*su[(ip*5)+0]+dvdx[0]*sx_tmp[0]+dvdx[9]*sx_tmp[1];
  sv[1] = dvdu[1]*su[(ip*5)+0]+dvdx[1]*sx_tmp[0]+dvdx[10]*sx_tmp[1];
  sv[2] = dvdx[11]*sx_tmp[1]+dvdx[20]*sx_tmp[2]+dvdx[29]*sx_tmp[3]+dvdx[38]*sx_tmp[4];
  sv[3] = dvdx[21]*sx_tmp[2]+dvdx[30]*sx_tmp[3]+dvdx[39]*sx_tmp[4]+dvdx[48]*sx_tmp[5];
  sv[4] = dvdx[22]*sx_tmp[2]+dvdx[31]*sx_tmp[3]+dvdx[40]*sx_tmp[4];
  sv[5] = dvdx[41]*sx_tmp[4]+dvdx[50]*sx_tmp[5]+dvdx[59]*sx_tmp[6]+dvdx[68]*sx_tmp[7];
  sv[6] = dvdx[51]*sx_tmp[5]+dvdx[60]*sx_tmp[6]+dvdx[69]*sx_tmp[7]+dvdx[78]*sx_tmp[8];
  sv[7] = dvdx[52]*sx_tmp[5]+dvdx[61]*sx_tmp[6]+dvdx[70]*sx_tmp[7];
  sv[8] = dvdx[71]*sx_tmp[7]+dvdx[80]*sx_tmp[8];
  switch (ip) {
    case 0: {
      sv[0] += dvdp[0];
    } break;
    case 1: {

    } break;
    case 2: {

    } break;
    case 3: {

    } break;
    case 4: {

    } break;
    case 5: {

    } break;
    case 6: {

    } break;
    case 7: {
      sv[0] += dvdp[63];
    } break;
    case 8: {
      sv[0] += dvdp[72];
      sv[1] += dvdp[73];
    } break;
    case 9: {
      sv[0] += dvdp[81];
      sv[1] += dvdp[82];
    } break;
    case 10: {
      sv[1] += dvdp[91];
      sv[2] += dvdp[92];
    } break;
    case 11: {
      sv[2] += dvdp[101];
    } break;
    case 12: {
      sv[2] += dvdp[110];
      sv[3] += dvdp[111];
      sv[4] += dvdp[112];
    } break;
    case 13: {
      sv[2] += dvdp[119];
      sv[3] += dvdp[120];
      sv[4] += dvdp[121];
    } break;
    case 14: {
      sv[2] += dvdp[128];
      sv[4] += dvdp[130];
      sv[5] += dvdp[131];
    } break;
    case 15: {
      sv[3] += dvdp[138];
      sv[5] += dvdp[140];
    } break;
    case 16: {
      sv[5] += dvdp[149];
      sv[6] += dvdp[150];
      sv[7] += dvdp[151];
    } break;
    case 17: {
      sv[5] += dvdp[158];
      sv[6] += dvdp[159];
      sv[7] += dvdp[160];
    } break;
    case 18: {
      sv[5] += dvdp[167];
      sv[7] += dvdp[169];
      sv[8] += dvdp[170];
    } break;
    case 19: {
      sv[6] += dvdp[177];
      sv[8] += dvdp[179];
    } break;
  }
  sxdot_tmp[0] = sv[0];
  sxdot_tmp[1] = sv[1];
  sxdot_tmp[2] = sv[2];
  sxdot_tmp[3] = sv[3];
  sxdot_tmp[4] = sv[4];
  sxdot_tmp[5] = sv[5];
  sxdot_tmp[6] = sv[6];
  sxdot_tmp[7] = sv[7];
  sxdot_tmp[8] = sv[8];
  for (is=0; is<9; is++) {
    if(mxIsNaN(sxdot_tmp[is])) sxdot_tmp[is] = 0.0;
  }

  return(0);
}


 int subfsx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(int Ns, realtype t, N_Vector x, N_Vector xdot, 
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
 {
   UserData data = (UserData) user_data;
   return fsx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(Ns, t, x, xdot, data->sensIndices[ip], sx, sxdot, user_data, tmp1, tmp2);
 };

 void csv_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data)
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
  for (is=0; is<9; is++) {
    sv[is] = 0.0;
  }
  sv[0] = dvdu[0]*su[(ip*5)+0]+dvdx[0]*sx_tmp[0]+dvdx[9]*sx_tmp[1];
  sv[1] = dvdu[1]*su[(ip*5)+0]+dvdx[1]*sx_tmp[0]+dvdx[10]*sx_tmp[1];
  sv[2] = dvdx[11]*sx_tmp[1]+dvdx[20]*sx_tmp[2]+dvdx[29]*sx_tmp[3]+dvdx[38]*sx_tmp[4];
  sv[3] = dvdx[21]*sx_tmp[2]+dvdx[30]*sx_tmp[3]+dvdx[39]*sx_tmp[4]+dvdx[48]*sx_tmp[5];
  sv[4] = dvdx[22]*sx_tmp[2]+dvdx[31]*sx_tmp[3]+dvdx[40]*sx_tmp[4];
  sv[5] = dvdx[41]*sx_tmp[4]+dvdx[50]*sx_tmp[5]+dvdx[59]*sx_tmp[6]+dvdx[68]*sx_tmp[7];
  sv[6] = dvdx[51]*sx_tmp[5]+dvdx[60]*sx_tmp[6]+dvdx[69]*sx_tmp[7]+dvdx[78]*sx_tmp[8];
  sv[7] = dvdx[52]*sx_tmp[5]+dvdx[61]*sx_tmp[6]+dvdx[70]*sx_tmp[7];
  sv[8] = dvdx[71]*sx_tmp[7]+dvdx[80]*sx_tmp[8];
  switch (ip) {
    case 0: {
      sv[0] += dvdp[0];
    } break;
    case 1: {

    } break;
    case 2: {

    } break;
    case 3: {

    } break;
    case 4: {

    } break;
    case 5: {

    } break;
    case 6: {

    } break;
    case 7: {
      sv[0] += dvdp[63];
    } break;
    case 8: {
      sv[0] += dvdp[72];
      sv[1] += dvdp[73];
    } break;
    case 9: {
      sv[0] += dvdp[81];
      sv[1] += dvdp[82];
    } break;
    case 10: {
      sv[1] += dvdp[91];
      sv[2] += dvdp[92];
    } break;
    case 11: {
      sv[2] += dvdp[101];
    } break;
    case 12: {
      sv[2] += dvdp[110];
      sv[3] += dvdp[111];
      sv[4] += dvdp[112];
    } break;
    case 13: {
      sv[2] += dvdp[119];
      sv[3] += dvdp[120];
      sv[4] += dvdp[121];
    } break;
    case 14: {
      sv[2] += dvdp[128];
      sv[4] += dvdp[130];
      sv[5] += dvdp[131];
    } break;
    case 15: {
      sv[3] += dvdp[138];
      sv[5] += dvdp[140];
    } break;
    case 16: {
      sv[5] += dvdp[149];
      sv[6] += dvdp[150];
      sv[7] += dvdp[151];
    } break;
    case 17: {
      sv[5] += dvdp[158];
      sv[6] += dvdp[159];
      sv[7] += dvdp[160];
    } break;
    case 18: {
      sv[5] += dvdp[167];
      sv[7] += dvdp[169];
      sv[8] += dvdp[170];
    } break;
    case 19: {
      sv[6] += dvdp[177];
      sv[8] += dvdp[179];
    } break;
  }

  return;
}


 void fsx0_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(int ip, N_Vector sx0, void *user_data)
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
      sx0_tmp[6] = 1.0;
    } break;
  }

  return;
}


 void subfsx0_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(int ip, N_Vector sx0, void *user_data)
 {
   UserData data = (UserData) user_data;
   fsx0_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(data->sensIndices[ip], sx0, user_data);
 };

 void dfxdp0_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, double *dfxdp0, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdp = data->dvdp;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *x_tmp = N_VGetArrayPointer(x);
  dfxdp0[0] = dvdp[0];
  dfxdp0[11] = dvdx[29];
  dfxdp0[12] = dvdx[30];
  dfxdp0[13] = dvdx[31];
  dfxdp0[18] = dvdx[0];
  dfxdp0[19] = dvdx[1];
  dfxdp0[32] = dvdx[59];
  dfxdp0[33] = dvdx[60];
  dfxdp0[34] = dvdx[61];
  dfxdp0[63] = dvdp[63];
  dfxdp0[72] = dvdp[72];
  dfxdp0[73] = dvdp[73];
  dfxdp0[81] = dvdp[81];
  dfxdp0[82] = dvdp[82];
  dfxdp0[91] = dvdp[91];
  dfxdp0[92] = dvdp[92];
  dfxdp0[101] = dvdp[101];
  dfxdp0[110] = dvdp[110];
  dfxdp0[111] = dvdp[111];
  dfxdp0[112] = dvdp[112];
  dfxdp0[119] = dvdp[119];
  dfxdp0[120] = dvdp[120];
  dfxdp0[121] = dvdp[121];
  dfxdp0[128] = dvdp[128];
  dfxdp0[130] = dvdp[130];
  dfxdp0[131] = dvdp[131];
  dfxdp0[138] = dvdp[138];
  dfxdp0[140] = dvdp[140];
  dfxdp0[149] = dvdp[149];
  dfxdp0[150] = dvdp[150];
  dfxdp0[151] = dvdp[151];
  dfxdp0[158] = dvdp[158];
  dfxdp0[159] = dvdp[159];
  dfxdp0[160] = dvdp[160];
  dfxdp0[167] = dvdp[167];
  dfxdp0[169] = dvdp[169];
  dfxdp0[170] = dvdp[170];
  dfxdp0[177] = dvdp[177];
  dfxdp0[179] = dvdp[179];
  for (is=0; is<180; is++) {
    if(mxIsNaN(dfxdp0[is])) dfxdp0[is] = 0.0;
  }

  return;
}


 void dfxdp_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(realtype t, N_Vector x, double *dfxdp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdp = data->dvdp;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *x_tmp = N_VGetArrayPointer(x);
  dfxdp[0] = dvdp[0];
  dfxdp[63] = dvdp[63];
  dfxdp[72] = dvdp[72];
  dfxdp[73] = dvdp[73];
  dfxdp[81] = dvdp[81];
  dfxdp[82] = dvdp[82];
  dfxdp[91] = dvdp[91];
  dfxdp[92] = dvdp[92];
  dfxdp[101] = dvdp[101];
  dfxdp[110] = dvdp[110];
  dfxdp[111] = dvdp[111];
  dfxdp[112] = dvdp[112];
  dfxdp[119] = dvdp[119];
  dfxdp[120] = dvdp[120];
  dfxdp[121] = dvdp[121];
  dfxdp[128] = dvdp[128];
  dfxdp[130] = dvdp[130];
  dfxdp[131] = dvdp[131];
  dfxdp[138] = dvdp[138];
  dfxdp[140] = dvdp[140];
  dfxdp[149] = dvdp[149];
  dfxdp[150] = dvdp[150];
  dfxdp[151] = dvdp[151];
  dfxdp[158] = dvdp[158];
  dfxdp[159] = dvdp[159];
  dfxdp[160] = dvdp[160];
  dfxdp[167] = dvdp[167];
  dfxdp[169] = dvdp[169];
  dfxdp[170] = dvdp[170];
  dfxdp[177] = dvdp[177];
  dfxdp[179] = dvdp[179];
  for (is=0; is<180; is++) {
    if(mxIsNaN(dfxdp[is])) dfxdp[is] = 0.0;
  }

  return;
}


 void fz_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x){
  z[nz*nt*iruns+it+nt*0] = (x[nx*nt*iruns+it+nt*2]+x[nx*nt*iruns+it+nt*4])*p[5];
  z[nz*nt*iruns+it+nt*1] = (x[nx*nt*iruns+it+nt*5]+x[nx*nt*iruns+it+nt*7])*p[4];
  z[nz*nt*iruns+it+nt*2] = p[6]*x[nx*nt*iruns+it+nt*8];

  return;
}


 void fsz_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx){
  int jp;
  for (jp=0; jp<np; jp++) {
      sz[it + nt*3*jp + nt*0] = p[5]*sx[it + nt*9*jp + nt*2]+p[5]*sx[it + nt*9*jp + nt*4];
      sz[it + nt*3*jp + nt*1] = p[4]*sx[it + nt*9*jp + nt*5]+p[4]*sx[it + nt*9*jp + nt*7];
      sz[it + nt*3*jp + nt*2] = p[6]*sx[it + nt*9*jp + nt*8];
  };

  sz[it+nt*13] += x[it+nt*5]+x[it+nt*7];
  sz[it+nt*15] += x[it+nt*2]+x[it+nt*4];
  sz[it+nt*20] += x[it+nt*8];

  return;
}


 void dfzdx_AktPathwayFujita_CE2DE45B61DFCA511ED306E2BB3E568A(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x){
      dfzdxs[nx*nt*iruns+it+nt*6] = p[5];
      dfzdxs[nx*nt*iruns+it+nt*12] = p[5];
      dfzdxs[nx*nt*iruns+it+nt*16] = p[4];
      dfzdxs[nx*nt*iruns+it+nt*22] = p[4];
      dfzdxs[nx*nt*iruns+it+nt*26] = p[6];


  return;
}


