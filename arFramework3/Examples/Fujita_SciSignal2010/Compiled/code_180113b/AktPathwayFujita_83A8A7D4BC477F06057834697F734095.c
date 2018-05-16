#include "AktPathwayFujita_83A8A7D4BC477F06057834697F734095.h"
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





 void fu_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(void *user_data, double t)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  data->u[1] = 3.0E1;
  data->u[3] = 6.0E1;
  data->u[4] = 3.6E3;

  return;
}


 void fsu_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(void *user_data, double t)
{

  return;
}


 void fv_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->v[1] = p[5]*x_tmp[2]+(x_tmp[0]-x_tmp[1])*p[0]-p[4]*x_tmp[1]*x_tmp[10];
  data->v[2] = -p[5]*x_tmp[2]-p[6]*x_tmp[2]+p[4]*x_tmp[1]*x_tmp[10];
  data->v[3] = p[6]*x_tmp[2]-p[7]*x_tmp[3]+p[9]*x_tmp[5]+p[10]*x_tmp[5]-p[8]*x_tmp[3]*x_tmp[4];
  data->v[4] = p[9]*x_tmp[5]+p[11]*x_tmp[6]-p[8]*x_tmp[3]*x_tmp[4];
  data->v[5] = -p[9]*x_tmp[5]-p[10]*x_tmp[5]+p[8]*x_tmp[3]*x_tmp[4];
  data->v[6] = p[10]*x_tmp[5]-p[11]*x_tmp[6]+p[13]*x_tmp[8]+p[14]*x_tmp[8]-p[12]*x_tmp[6]*x_tmp[7];
  data->v[7] = p[13]*x_tmp[8]+p[15]*x_tmp[9]-p[12]*x_tmp[6]*x_tmp[7];
  data->v[8] = -p[13]*x_tmp[8]-p[14]*x_tmp[8]+p[12]*x_tmp[6]*x_tmp[7];
  data->v[9] = p[14]*x_tmp[8]-p[15]*x_tmp[9];

  return;
}


 void dvdx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdx[1] = p[0];
  data->dvdx[12] = -p[4]*x_tmp[10]-p[0];
  data->dvdx[13] = p[4]*x_tmp[10];
  data->dvdx[23] = p[5];
  data->dvdx[24] = -p[5]-p[6];
  data->dvdx[25] = p[6];
  data->dvdx[36] = -p[8]*x_tmp[4]-p[7];
  data->dvdx[37] = -p[8]*x_tmp[4];
  data->dvdx[38] = p[8]*x_tmp[4];
  data->dvdx[47] = -p[8]*x_tmp[3];
  data->dvdx[48] = -p[8]*x_tmp[3];
  data->dvdx[49] = p[8]*x_tmp[3];
  data->dvdx[58] = p[9]+p[10];
  data->dvdx[59] = p[9];
  data->dvdx[60] = -p[9]-p[10];
  data->dvdx[61] = p[10];
  data->dvdx[70] = p[11];
  data->dvdx[72] = -p[12]*x_tmp[7]-p[11];
  data->dvdx[73] = -p[12]*x_tmp[7];
  data->dvdx[74] = p[12]*x_tmp[7];
  data->dvdx[83] = -p[12]*x_tmp[6];
  data->dvdx[84] = -p[12]*x_tmp[6];
  data->dvdx[85] = p[12]*x_tmp[6];
  data->dvdx[94] = p[13]+p[14];
  data->dvdx[95] = p[13];
  data->dvdx[96] = -p[13]-p[14];
  data->dvdx[97] = p[14];
  data->dvdx[106] = p[15];
  data->dvdx[108] = -p[15];
  data->dvdx[111] = -p[4]*x_tmp[1];
  data->dvdx[112] = p[4]*x_tmp[1];

  return;
}


 void dvdu_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);


  return;
}


 void dvdp_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdp[1] = x_tmp[0]-x_tmp[1];
  data->dvdp[45] = -x_tmp[1]*x_tmp[10];
  data->dvdp[46] = x_tmp[1]*x_tmp[10];
  data->dvdp[56] = x_tmp[2];
  data->dvdp[57] = -x_tmp[2];
  data->dvdp[68] = -x_tmp[2];
  data->dvdp[69] = x_tmp[2];
  data->dvdp[80] = -x_tmp[3];
  data->dvdp[91] = -x_tmp[3]*x_tmp[4];
  data->dvdp[92] = -x_tmp[3]*x_tmp[4];
  data->dvdp[93] = x_tmp[3]*x_tmp[4];
  data->dvdp[102] = x_tmp[5];
  data->dvdp[103] = x_tmp[5];
  data->dvdp[104] = -x_tmp[5];
  data->dvdp[113] = x_tmp[5];
  data->dvdp[115] = -x_tmp[5];
  data->dvdp[116] = x_tmp[5];
  data->dvdp[125] = x_tmp[6];
  data->dvdp[127] = -x_tmp[6];
  data->dvdp[138] = -x_tmp[6]*x_tmp[7];
  data->dvdp[139] = -x_tmp[6]*x_tmp[7];
  data->dvdp[140] = x_tmp[6]*x_tmp[7];
  data->dvdp[149] = x_tmp[8];
  data->dvdp[150] = x_tmp[8];
  data->dvdp[151] = -x_tmp[8];
  data->dvdp[160] = x_tmp[8];
  data->dvdp[162] = -x_tmp[8];
  data->dvdp[163] = x_tmp[8];
  data->dvdp[172] = x_tmp[9];
  data->dvdp[174] = -x_tmp[9];

  return;
}


 int fx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  UserData data = (UserData) user_data;
  int is;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  fu_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(data, t);
  fv_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(t, x, data);
  xdot_tmp[0] = v[0];
  xdot_tmp[1] = v[1];
  xdot_tmp[2] = v[2];
  xdot_tmp[3] = v[3];
  xdot_tmp[4] = v[4];
  xdot_tmp[5] = v[5];
  xdot_tmp[6] = v[6];
  xdot_tmp[7] = v[7];
  xdot_tmp[8] = v[8];
  xdot_tmp[9] = v[9];
  xdot_tmp[10] = v[10];
  for (is=0; is<11; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
    if(qpositivex[is]>0.5 && x_tmp[is]<0.0 && xdot_tmp[is]<0.0) xdot_tmp[is] = -xdot_tmp[is];
  }

  return(*(data->abort));
}


 void fxdouble_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, double *xdot_tmp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  fu_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(data, t);
  fv_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(t, x, data);
  xdot_tmp[0] = v[0];
  xdot_tmp[1] = v[1];
  xdot_tmp[2] = v[2];
  xdot_tmp[3] = v[3];
  xdot_tmp[4] = v[4];
  xdot_tmp[5] = v[5];
  xdot_tmp[6] = v[6];
  xdot_tmp[7] = v[7];
  xdot_tmp[8] = v[8];
  xdot_tmp[9] = v[9];
  xdot_tmp[10] = v[10];
  for (is=0; is<11; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
  }

  return;
}


 void fx0_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  x0_tmp[0] = 6.81901837333797E4;
  x0_tmp[1] = 6.81901837333797E4;
  x0_tmp[4] = 4.33090165709309E-2;
  x0_tmp[7] = 3.54316740542218;
  x0_tmp[10] = 3.0E1;

  return;
}


 int dfxdx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(long int N, realtype t, N_Vector x, 
  	N_Vector fx, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(t, x, data);
  for (is=0; is<121; is++) {
    J->data[is] = 0.0;
  }
  J->data[1] = dvdx[1];
  J->data[12] = dvdx[12];
  J->data[13] = dvdx[13];
  J->data[23] = dvdx[23];
  J->data[24] = dvdx[24];
  J->data[25] = dvdx[25];
  J->data[36] = dvdx[36];
  J->data[37] = dvdx[37];
  J->data[38] = dvdx[38];
  J->data[47] = dvdx[47];
  J->data[48] = dvdx[48];
  J->data[49] = dvdx[49];
  J->data[58] = dvdx[58];
  J->data[59] = dvdx[59];
  J->data[60] = dvdx[60];
  J->data[61] = dvdx[61];
  J->data[70] = dvdx[70];
  J->data[72] = dvdx[72];
  J->data[73] = dvdx[73];
  J->data[74] = dvdx[74];
  J->data[83] = dvdx[83];
  J->data[84] = dvdx[84];
  J->data[85] = dvdx[85];
  J->data[94] = dvdx[94];
  J->data[95] = dvdx[95];
  J->data[96] = dvdx[96];
  J->data[97] = dvdx[97];
  J->data[106] = dvdx[106];
  J->data[108] = dvdx[108];
  J->data[111] = dvdx[111];
  J->data[112] = dvdx[112];
  for (is=0; is<121; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = 0.0;
  }

  return(0);
}


 int dfxdx_out_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, realtype* J, void *user_data) {
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(t, x, data);
  for (is=0; is<121; is++) {
    J[is] = 0.0;
  }
  J[1] = dvdx[1];
  J[12] = dvdx[12];
  J[13] = dvdx[13];
  J[23] = dvdx[23];
  J[24] = dvdx[24];
  J[25] = dvdx[25];
  J[36] = dvdx[36];
  J[37] = dvdx[37];
  J[38] = dvdx[38];
  J[47] = dvdx[47];
  J[48] = dvdx[48];
  J[49] = dvdx[49];
  J[58] = dvdx[58];
  J[59] = dvdx[59];
  J[60] = dvdx[60];
  J[61] = dvdx[61];
  J[70] = dvdx[70];
  J[72] = dvdx[72];
  J[73] = dvdx[73];
  J[74] = dvdx[74];
  J[83] = dvdx[83];
  J[84] = dvdx[84];
  J[85] = dvdx[85];
  J[94] = dvdx[94];
  J[95] = dvdx[95];
  J[96] = dvdx[96];
  J[97] = dvdx[97];
  J[106] = dvdx[106];
  J[108] = dvdx[108];
  J[111] = dvdx[111];
  J[112] = dvdx[112];
  for (is=0; is<121; is++) {
    if(mxIsNaN(J[is])) J[is] = 0.0;
  }

  return(0);
}


 int dfxdx_sparse_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, 
  	N_Vector fx, SlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(t, x, data);
  SlsSetToZero(J);
    J->rowvals[0] = 1;
    J->rowvals[1] = 1;
    J->rowvals[2] = 2;
    J->rowvals[3] = 1;
    J->rowvals[4] = 2;
    J->rowvals[5] = 3;
    J->rowvals[6] = 3;
    J->rowvals[7] = 4;
    J->rowvals[8] = 5;
    J->rowvals[9] = 3;
    J->rowvals[10] = 4;
    J->rowvals[11] = 5;
    J->rowvals[12] = 3;
    J->rowvals[13] = 4;
    J->rowvals[14] = 5;
    J->rowvals[15] = 6;
    J->rowvals[16] = 4;
    J->rowvals[17] = 6;
    J->rowvals[18] = 7;
    J->rowvals[19] = 8;
    J->rowvals[20] = 6;
    J->rowvals[21] = 7;
    J->rowvals[22] = 8;
    J->rowvals[23] = 6;
    J->rowvals[24] = 7;
    J->rowvals[25] = 8;
    J->rowvals[26] = 9;
    J->rowvals[27] = 7;
    J->rowvals[28] = 9;
    J->rowvals[29] = 1;
    J->rowvals[30] = 2;

    J->colptrs[0] = 0;
    J->colptrs[1] = 1;
    J->colptrs[2] = 3;
    J->colptrs[3] = 6;
    J->colptrs[4] = 9;
    J->colptrs[5] = 12;
    J->colptrs[6] = 16;
    J->colptrs[7] = 20;
    J->colptrs[8] = 23;
    J->colptrs[9] = 27;
    J->colptrs[10] = 29;
    J->colptrs[11] = 31;

  J->data[0] = dvdx[1];
  J->data[1] = dvdx[12];
  J->data[2] = dvdx[13];
  J->data[3] = dvdx[23];
  J->data[4] = dvdx[24];
  J->data[5] = dvdx[25];
  J->data[6] = dvdx[36];
  J->data[7] = dvdx[37];
  J->data[8] = dvdx[38];
  J->data[9] = dvdx[47];
  J->data[10] = dvdx[48];
  J->data[11] = dvdx[49];
  J->data[12] = dvdx[58];
  J->data[13] = dvdx[59];
  J->data[14] = dvdx[60];
  J->data[15] = dvdx[61];
  J->data[16] = dvdx[70];
  J->data[17] = dvdx[72];
  J->data[18] = dvdx[73];
  J->data[19] = dvdx[74];
  J->data[20] = dvdx[83];
  J->data[21] = dvdx[84];
  J->data[22] = dvdx[85];
  J->data[23] = dvdx[94];
  J->data[24] = dvdx[95];
  J->data[25] = dvdx[96];
  J->data[26] = dvdx[97];
  J->data[27] = dvdx[106];
  J->data[28] = dvdx[108];
  J->data[29] = dvdx[111];
  J->data[30] = dvdx[112];
  for (is=0; is<31; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = RCONST(0.0);
  }

  return(0);
}


 int fsx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(int Ns, realtype t, N_Vector x, N_Vector xdot, 
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
  fsu_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(data, t);
  dvdx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(t, x, data);
  dvdu_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(t, x, data);
  dvdp_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(t, x, data);
  for (is=0; is<11; is++) {
    sv[is] = 0.0;
  }
  sv[1] = dvdx[1]*sx_tmp[0]+dvdx[12]*sx_tmp[1]+dvdx[23]*sx_tmp[2]+dvdx[111]*sx_tmp[10];
  sv[2] = dvdx[13]*sx_tmp[1]+dvdx[24]*sx_tmp[2]+dvdx[112]*sx_tmp[10];
  sv[3] = dvdx[25]*sx_tmp[2]+dvdx[36]*sx_tmp[3]+dvdx[47]*sx_tmp[4]+dvdx[58]*sx_tmp[5];
  sv[4] = dvdx[37]*sx_tmp[3]+dvdx[48]*sx_tmp[4]+dvdx[59]*sx_tmp[5]+dvdx[70]*sx_tmp[6];
  sv[5] = dvdx[38]*sx_tmp[3]+dvdx[49]*sx_tmp[4]+dvdx[60]*sx_tmp[5];
  sv[6] = dvdx[61]*sx_tmp[5]+dvdx[72]*sx_tmp[6]+dvdx[83]*sx_tmp[7]+dvdx[94]*sx_tmp[8];
  sv[7] = dvdx[73]*sx_tmp[6]+dvdx[84]*sx_tmp[7]+dvdx[95]*sx_tmp[8]+dvdx[106]*sx_tmp[9];
  sv[8] = dvdx[74]*sx_tmp[6]+dvdx[85]*sx_tmp[7]+dvdx[96]*sx_tmp[8];
  sv[9] = dvdx[97]*sx_tmp[8]+dvdx[108]*sx_tmp[9];
  switch (ip) {
    case 0: {
      sv[1] += dvdp[1];
    } break;
    case 1: {

    } break;
    case 2: {

    } break;
    case 3: {

    } break;
    case 4: {
      sv[1] += dvdp[45];
      sv[2] += dvdp[46];
    } break;
    case 5: {
      sv[1] += dvdp[56];
      sv[2] += dvdp[57];
    } break;
    case 6: {
      sv[2] += dvdp[68];
      sv[3] += dvdp[69];
    } break;
    case 7: {
      sv[3] += dvdp[80];
    } break;
    case 8: {
      sv[3] += dvdp[91];
      sv[4] += dvdp[92];
      sv[5] += dvdp[93];
    } break;
    case 9: {
      sv[3] += dvdp[102];
      sv[4] += dvdp[103];
      sv[5] += dvdp[104];
    } break;
    case 10: {
      sv[3] += dvdp[113];
      sv[5] += dvdp[115];
      sv[6] += dvdp[116];
    } break;
    case 11: {
      sv[4] += dvdp[125];
      sv[6] += dvdp[127];
    } break;
    case 12: {
      sv[6] += dvdp[138];
      sv[7] += dvdp[139];
      sv[8] += dvdp[140];
    } break;
    case 13: {
      sv[6] += dvdp[149];
      sv[7] += dvdp[150];
      sv[8] += dvdp[151];
    } break;
    case 14: {
      sv[6] += dvdp[160];
      sv[8] += dvdp[162];
      sv[9] += dvdp[163];
    } break;
    case 15: {
      sv[7] += dvdp[172];
      sv[9] += dvdp[174];
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
  sxdot_tmp[9] = sv[9];
  sxdot_tmp[10] = sv[10];
  for (is=0; is<11; is++) {
    if(mxIsNaN(sxdot_tmp[is])) sxdot_tmp[is] = 0.0;
  }

  return(0);
}


 int subfsx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(int Ns, realtype t, N_Vector x, N_Vector xdot, 
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
 {
   UserData data = (UserData) user_data;
   return fsx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(Ns, t, x, xdot, data->sensIndices[ip], sx, sxdot, user_data, tmp1, tmp2);
 };

 void csv_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data)
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
  sv[1] = dvdx[1]*sx_tmp[0]+dvdx[12]*sx_tmp[1]+dvdx[23]*sx_tmp[2]+dvdx[111]*sx_tmp[10];
  sv[2] = dvdx[13]*sx_tmp[1]+dvdx[24]*sx_tmp[2]+dvdx[112]*sx_tmp[10];
  sv[3] = dvdx[25]*sx_tmp[2]+dvdx[36]*sx_tmp[3]+dvdx[47]*sx_tmp[4]+dvdx[58]*sx_tmp[5];
  sv[4] = dvdx[37]*sx_tmp[3]+dvdx[48]*sx_tmp[4]+dvdx[59]*sx_tmp[5]+dvdx[70]*sx_tmp[6];
  sv[5] = dvdx[38]*sx_tmp[3]+dvdx[49]*sx_tmp[4]+dvdx[60]*sx_tmp[5];
  sv[6] = dvdx[61]*sx_tmp[5]+dvdx[72]*sx_tmp[6]+dvdx[83]*sx_tmp[7]+dvdx[94]*sx_tmp[8];
  sv[7] = dvdx[73]*sx_tmp[6]+dvdx[84]*sx_tmp[7]+dvdx[95]*sx_tmp[8]+dvdx[106]*sx_tmp[9];
  sv[8] = dvdx[74]*sx_tmp[6]+dvdx[85]*sx_tmp[7]+dvdx[96]*sx_tmp[8];
  sv[9] = dvdx[97]*sx_tmp[8]+dvdx[108]*sx_tmp[9];
  switch (ip) {
    case 0: {
      sv[1] += dvdp[1];
    } break;
    case 1: {

    } break;
    case 2: {

    } break;
    case 3: {

    } break;
    case 4: {
      sv[1] += dvdp[45];
      sv[2] += dvdp[46];
    } break;
    case 5: {
      sv[1] += dvdp[56];
      sv[2] += dvdp[57];
    } break;
    case 6: {
      sv[2] += dvdp[68];
      sv[3] += dvdp[69];
    } break;
    case 7: {
      sv[3] += dvdp[80];
    } break;
    case 8: {
      sv[3] += dvdp[91];
      sv[4] += dvdp[92];
      sv[5] += dvdp[93];
    } break;
    case 9: {
      sv[3] += dvdp[102];
      sv[4] += dvdp[103];
      sv[5] += dvdp[104];
    } break;
    case 10: {
      sv[3] += dvdp[113];
      sv[5] += dvdp[115];
      sv[6] += dvdp[116];
    } break;
    case 11: {
      sv[4] += dvdp[125];
      sv[6] += dvdp[127];
    } break;
    case 12: {
      sv[6] += dvdp[138];
      sv[7] += dvdp[139];
      sv[8] += dvdp[140];
    } break;
    case 13: {
      sv[6] += dvdp[149];
      sv[7] += dvdp[150];
      sv[8] += dvdp[151];
    } break;
    case 14: {
      sv[6] += dvdp[160];
      sv[8] += dvdp[162];
      sv[9] += dvdp[163];
    } break;
    case 15: {
      sv[7] += dvdp[172];
      sv[9] += dvdp[174];
    } break;
  }

  return;
}


 void fsx0_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  switch (ip) {
  }

  return;
}


 void subfsx0_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(int ip, N_Vector sx0, void *user_data)
 {
   UserData data = (UserData) user_data;
   fsx0_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(data->sensIndices[ip], sx0, user_data);
 };

 void dfxdp0_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, double *dfxdp0, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdp = data->dvdp;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *x_tmp = N_VGetArrayPointer(x);
  dfxdp0[1] = dvdp[1];
  dfxdp0[45] = dvdp[45];
  dfxdp0[46] = dvdp[46];
  dfxdp0[56] = dvdp[56];
  dfxdp0[57] = dvdp[57];
  dfxdp0[68] = dvdp[68];
  dfxdp0[69] = dvdp[69];
  dfxdp0[80] = dvdp[80];
  dfxdp0[91] = dvdp[91];
  dfxdp0[92] = dvdp[92];
  dfxdp0[93] = dvdp[93];
  dfxdp0[102] = dvdp[102];
  dfxdp0[103] = dvdp[103];
  dfxdp0[104] = dvdp[104];
  dfxdp0[113] = dvdp[113];
  dfxdp0[115] = dvdp[115];
  dfxdp0[116] = dvdp[116];
  dfxdp0[125] = dvdp[125];
  dfxdp0[127] = dvdp[127];
  dfxdp0[138] = dvdp[138];
  dfxdp0[139] = dvdp[139];
  dfxdp0[140] = dvdp[140];
  dfxdp0[149] = dvdp[149];
  dfxdp0[150] = dvdp[150];
  dfxdp0[151] = dvdp[151];
  dfxdp0[160] = dvdp[160];
  dfxdp0[162] = dvdp[162];
  dfxdp0[163] = dvdp[163];
  dfxdp0[172] = dvdp[172];
  dfxdp0[174] = dvdp[174];
  for (is=0; is<176; is++) {
    if(mxIsNaN(dfxdp0[is])) dfxdp0[is] = 0.0;
  }

  return;
}


 void dfxdp_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(realtype t, N_Vector x, double *dfxdp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdp = data->dvdp;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *x_tmp = N_VGetArrayPointer(x);
  dfxdp[1] = dvdp[1];
  dfxdp[45] = dvdp[45];
  dfxdp[46] = dvdp[46];
  dfxdp[56] = dvdp[56];
  dfxdp[57] = dvdp[57];
  dfxdp[68] = dvdp[68];
  dfxdp[69] = dvdp[69];
  dfxdp[80] = dvdp[80];
  dfxdp[91] = dvdp[91];
  dfxdp[92] = dvdp[92];
  dfxdp[93] = dvdp[93];
  dfxdp[102] = dvdp[102];
  dfxdp[103] = dvdp[103];
  dfxdp[104] = dvdp[104];
  dfxdp[113] = dvdp[113];
  dfxdp[115] = dvdp[115];
  dfxdp[116] = dvdp[116];
  dfxdp[125] = dvdp[125];
  dfxdp[127] = dvdp[127];
  dfxdp[138] = dvdp[138];
  dfxdp[139] = dvdp[139];
  dfxdp[140] = dvdp[140];
  dfxdp[149] = dvdp[149];
  dfxdp[150] = dvdp[150];
  dfxdp[151] = dvdp[151];
  dfxdp[160] = dvdp[160];
  dfxdp[162] = dvdp[162];
  dfxdp[163] = dvdp[163];
  dfxdp[172] = dvdp[172];
  dfxdp[174] = dvdp[174];
  for (is=0; is<176; is++) {
    if(mxIsNaN(dfxdp[is])) dfxdp[is] = 0.0;
  }

  return;
}


 void fz_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x){
  z[nz*nt*iruns+it+nt*0] = (x[nx*nt*iruns+it+nt*3]+x[nx*nt*iruns+it+nt*5])*p[2];
  z[nz*nt*iruns+it+nt*1] = (x[nx*nt*iruns+it+nt*6]+x[nx*nt*iruns+it+nt*8])*p[1];
  z[nz*nt*iruns+it+nt*2] = p[3]*x[nx*nt*iruns+it+nt*9];

  return;
}


 void fsz_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx){
  int jp;
  for (jp=0; jp<np; jp++) {
      sz[it + nt*3*jp + nt*0] = p[2]*sx[it + nt*11*jp + nt*3]+p[2]*sx[it + nt*11*jp + nt*5];
      sz[it + nt*3*jp + nt*1] = p[1]*sx[it + nt*11*jp + nt*6]+p[1]*sx[it + nt*11*jp + nt*8];
      sz[it + nt*3*jp + nt*2] = p[3]*sx[it + nt*11*jp + nt*9];
  };

  sz[it+nt*4] += x[it+nt*6]+x[it+nt*8];
  sz[it+nt*6] += x[it+nt*3]+x[it+nt*5];
  sz[it+nt*11] += x[it+nt*9];

  return;
}


 void dfzdx_AktPathwayFujita_83A8A7D4BC477F06057834697F734095(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x){
      dfzdxs[nx*nt*iruns+it+nt*9] = p[2];
      dfzdxs[nx*nt*iruns+it+nt*15] = p[2];
      dfzdxs[nx*nt*iruns+it+nt*19] = p[1];
      dfzdxs[nx*nt*iruns+it+nt*25] = p[1];
      dfzdxs[nx*nt*iruns+it+nt*29] = p[3];


  return;
}


