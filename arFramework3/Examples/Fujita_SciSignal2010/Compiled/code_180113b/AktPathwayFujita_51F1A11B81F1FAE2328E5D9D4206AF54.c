#include "AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54.h"
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





 void fu_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(void *user_data, double t)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  data->u[0] = 3.0;
  data->u[3] = 6.0E1;
  data->u[4] = 3.6E3;

  return;
}


 void fsu_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(void *user_data, double t)
{

  return;
}


 void fv_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->v[1] = p[6]*x_tmp[2]+(x_tmp[0]-x_tmp[1])*p[0]-p[5]*u[0]*x_tmp[1];
  data->v[2] = -p[6]*x_tmp[2]-p[7]*x_tmp[2]+p[5]*u[0]*x_tmp[1];
  data->v[3] = p[7]*x_tmp[2]-p[8]*x_tmp[3]+p[10]*x_tmp[5]+p[11]*x_tmp[5]-p[9]*x_tmp[3]*x_tmp[4];
  data->v[4] = p[10]*x_tmp[5]+p[12]*x_tmp[6]-p[9]*x_tmp[3]*x_tmp[4];
  data->v[5] = -p[10]*x_tmp[5]-p[11]*x_tmp[5]+p[9]*x_tmp[3]*x_tmp[4];
  data->v[6] = p[11]*x_tmp[5]-p[12]*x_tmp[6]+p[14]*x_tmp[8]+p[15]*x_tmp[8]-p[13]*x_tmp[6]*x_tmp[7];
  data->v[7] = p[14]*x_tmp[8]+p[16]*x_tmp[9]-p[13]*x_tmp[6]*x_tmp[7];
  data->v[8] = -p[14]*x_tmp[8]-p[15]*x_tmp[8]+p[13]*x_tmp[6]*x_tmp[7];
  data->v[9] = p[15]*x_tmp[8]-p[16]*x_tmp[9];

  return;
}


 void dvdx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdx[1] = p[0];
  data->dvdx[11] = -p[5]*u[0]-p[0];
  data->dvdx[12] = p[5]*u[0];
  data->dvdx[21] = p[6];
  data->dvdx[22] = -p[6]-p[7];
  data->dvdx[23] = p[7];
  data->dvdx[33] = -p[9]*x_tmp[4]-p[8];
  data->dvdx[34] = -p[9]*x_tmp[4];
  data->dvdx[35] = p[9]*x_tmp[4];
  data->dvdx[43] = -p[9]*x_tmp[3];
  data->dvdx[44] = -p[9]*x_tmp[3];
  data->dvdx[45] = p[9]*x_tmp[3];
  data->dvdx[53] = p[10]+p[11];
  data->dvdx[54] = p[10];
  data->dvdx[55] = -p[10]-p[11];
  data->dvdx[56] = p[11];
  data->dvdx[64] = p[12];
  data->dvdx[66] = -p[13]*x_tmp[7]-p[12];
  data->dvdx[67] = -p[13]*x_tmp[7];
  data->dvdx[68] = p[13]*x_tmp[7];
  data->dvdx[76] = -p[13]*x_tmp[6];
  data->dvdx[77] = -p[13]*x_tmp[6];
  data->dvdx[78] = p[13]*x_tmp[6];
  data->dvdx[86] = p[14]+p[15];
  data->dvdx[87] = p[14];
  data->dvdx[88] = -p[14]-p[15];
  data->dvdx[89] = p[15];
  data->dvdx[97] = p[16];
  data->dvdx[99] = -p[16];

  return;
}


 void dvdu_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdu[1] = -p[5]*x_tmp[1];
  data->dvdu[2] = p[5]*x_tmp[1];

  return;
}


 void dvdp_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdp[1] = x_tmp[0]-x_tmp[1];
  data->dvdp[51] = -u[0]*x_tmp[1];
  data->dvdp[52] = u[0]*x_tmp[1];
  data->dvdp[61] = x_tmp[2];
  data->dvdp[62] = -x_tmp[2];
  data->dvdp[72] = -x_tmp[2];
  data->dvdp[73] = x_tmp[2];
  data->dvdp[83] = -x_tmp[3];
  data->dvdp[93] = -x_tmp[3]*x_tmp[4];
  data->dvdp[94] = -x_tmp[3]*x_tmp[4];
  data->dvdp[95] = x_tmp[3]*x_tmp[4];
  data->dvdp[103] = x_tmp[5];
  data->dvdp[104] = x_tmp[5];
  data->dvdp[105] = -x_tmp[5];
  data->dvdp[113] = x_tmp[5];
  data->dvdp[115] = -x_tmp[5];
  data->dvdp[116] = x_tmp[5];
  data->dvdp[124] = x_tmp[6];
  data->dvdp[126] = -x_tmp[6];
  data->dvdp[136] = -x_tmp[6]*x_tmp[7];
  data->dvdp[137] = -x_tmp[6]*x_tmp[7];
  data->dvdp[138] = x_tmp[6]*x_tmp[7];
  data->dvdp[146] = x_tmp[8];
  data->dvdp[147] = x_tmp[8];
  data->dvdp[148] = -x_tmp[8];
  data->dvdp[156] = x_tmp[8];
  data->dvdp[158] = -x_tmp[8];
  data->dvdp[159] = x_tmp[8];
  data->dvdp[167] = x_tmp[9];
  data->dvdp[169] = -x_tmp[9];

  return;
}


 int fx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  UserData data = (UserData) user_data;
  int is;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  fu_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(data, t);
  fv_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(t, x, data);
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
  for (is=0; is<10; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
    if(qpositivex[is]>0.5 && x_tmp[is]<0.0 && xdot_tmp[is]<0.0) xdot_tmp[is] = -xdot_tmp[is];
  }

  return(*(data->abort));
}


 void fxdouble_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, double *xdot_tmp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  fu_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(data, t);
  fv_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(t, x, data);
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
  for (is=0; is<10; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
  }

  return;
}


 void fx0_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  x0_tmp[0] = 6.819E4;
  x0_tmp[1] = p[1];
  x0_tmp[4] = 4.33E-2;
  x0_tmp[7] = 3.5432;

  return;
}


 int dfxdx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(long int N, realtype t, N_Vector x, 
  	N_Vector fx, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(t, x, data);
  for (is=0; is<100; is++) {
    J->data[is] = 0.0;
  }
  J->data[1] = dvdx[1];
  J->data[11] = dvdx[11];
  J->data[12] = dvdx[12];
  J->data[21] = dvdx[21];
  J->data[22] = dvdx[22];
  J->data[23] = dvdx[23];
  J->data[33] = dvdx[33];
  J->data[34] = dvdx[34];
  J->data[35] = dvdx[35];
  J->data[43] = dvdx[43];
  J->data[44] = dvdx[44];
  J->data[45] = dvdx[45];
  J->data[53] = dvdx[53];
  J->data[54] = dvdx[54];
  J->data[55] = dvdx[55];
  J->data[56] = dvdx[56];
  J->data[64] = dvdx[64];
  J->data[66] = dvdx[66];
  J->data[67] = dvdx[67];
  J->data[68] = dvdx[68];
  J->data[76] = dvdx[76];
  J->data[77] = dvdx[77];
  J->data[78] = dvdx[78];
  J->data[86] = dvdx[86];
  J->data[87] = dvdx[87];
  J->data[88] = dvdx[88];
  J->data[89] = dvdx[89];
  J->data[97] = dvdx[97];
  J->data[99] = dvdx[99];
  for (is=0; is<100; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = 0.0;
  }

  return(0);
}


 int dfxdx_out_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, realtype* J, void *user_data) {
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(t, x, data);
  for (is=0; is<100; is++) {
    J[is] = 0.0;
  }
  J[1] = dvdx[1];
  J[11] = dvdx[11];
  J[12] = dvdx[12];
  J[21] = dvdx[21];
  J[22] = dvdx[22];
  J[23] = dvdx[23];
  J[33] = dvdx[33];
  J[34] = dvdx[34];
  J[35] = dvdx[35];
  J[43] = dvdx[43];
  J[44] = dvdx[44];
  J[45] = dvdx[45];
  J[53] = dvdx[53];
  J[54] = dvdx[54];
  J[55] = dvdx[55];
  J[56] = dvdx[56];
  J[64] = dvdx[64];
  J[66] = dvdx[66];
  J[67] = dvdx[67];
  J[68] = dvdx[68];
  J[76] = dvdx[76];
  J[77] = dvdx[77];
  J[78] = dvdx[78];
  J[86] = dvdx[86];
  J[87] = dvdx[87];
  J[88] = dvdx[88];
  J[89] = dvdx[89];
  J[97] = dvdx[97];
  J[99] = dvdx[99];
  for (is=0; is<100; is++) {
    if(mxIsNaN(J[is])) J[is] = 0.0;
  }

  return(0);
}


 int dfxdx_sparse_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, 
  	N_Vector fx, SlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(t, x, data);
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

  J->data[0] = dvdx[1];
  J->data[1] = dvdx[11];
  J->data[2] = dvdx[12];
  J->data[3] = dvdx[21];
  J->data[4] = dvdx[22];
  J->data[5] = dvdx[23];
  J->data[6] = dvdx[33];
  J->data[7] = dvdx[34];
  J->data[8] = dvdx[35];
  J->data[9] = dvdx[43];
  J->data[10] = dvdx[44];
  J->data[11] = dvdx[45];
  J->data[12] = dvdx[53];
  J->data[13] = dvdx[54];
  J->data[14] = dvdx[55];
  J->data[15] = dvdx[56];
  J->data[16] = dvdx[64];
  J->data[17] = dvdx[66];
  J->data[18] = dvdx[67];
  J->data[19] = dvdx[68];
  J->data[20] = dvdx[76];
  J->data[21] = dvdx[77];
  J->data[22] = dvdx[78];
  J->data[23] = dvdx[86];
  J->data[24] = dvdx[87];
  J->data[25] = dvdx[88];
  J->data[26] = dvdx[89];
  J->data[27] = dvdx[97];
  J->data[28] = dvdx[99];
  for (is=0; is<29; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = RCONST(0.0);
  }

  return(0);
}


 int fsx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(int Ns, realtype t, N_Vector x, N_Vector xdot, 
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
  fsu_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(data, t);
  dvdx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(t, x, data);
  dvdu_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(t, x, data);
  dvdp_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(t, x, data);
  for (is=0; is<10; is++) {
    sv[is] = 0.0;
  }
  sv[1] = dvdu[1]*su[(ip*5)+0]+dvdx[1]*sx_tmp[0]+dvdx[11]*sx_tmp[1]+dvdx[21]*sx_tmp[2];
  sv[2] = dvdu[2]*su[(ip*5)+0]+dvdx[12]*sx_tmp[1]+dvdx[22]*sx_tmp[2];
  sv[3] = dvdx[23]*sx_tmp[2]+dvdx[33]*sx_tmp[3]+dvdx[43]*sx_tmp[4]+dvdx[53]*sx_tmp[5];
  sv[4] = dvdx[34]*sx_tmp[3]+dvdx[44]*sx_tmp[4]+dvdx[54]*sx_tmp[5]+dvdx[64]*sx_tmp[6];
  sv[5] = dvdx[35]*sx_tmp[3]+dvdx[45]*sx_tmp[4]+dvdx[55]*sx_tmp[5];
  sv[6] = dvdx[56]*sx_tmp[5]+dvdx[66]*sx_tmp[6]+dvdx[76]*sx_tmp[7]+dvdx[86]*sx_tmp[8];
  sv[7] = dvdx[67]*sx_tmp[6]+dvdx[77]*sx_tmp[7]+dvdx[87]*sx_tmp[8]+dvdx[97]*sx_tmp[9];
  sv[8] = dvdx[68]*sx_tmp[6]+dvdx[78]*sx_tmp[7]+dvdx[88]*sx_tmp[8];
  sv[9] = dvdx[89]*sx_tmp[8]+dvdx[99]*sx_tmp[9];
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

    } break;
    case 5: {
      sv[1] += dvdp[51];
      sv[2] += dvdp[52];
    } break;
    case 6: {
      sv[1] += dvdp[61];
      sv[2] += dvdp[62];
    } break;
    case 7: {
      sv[2] += dvdp[72];
      sv[3] += dvdp[73];
    } break;
    case 8: {
      sv[3] += dvdp[83];
    } break;
    case 9: {
      sv[3] += dvdp[93];
      sv[4] += dvdp[94];
      sv[5] += dvdp[95];
    } break;
    case 10: {
      sv[3] += dvdp[103];
      sv[4] += dvdp[104];
      sv[5] += dvdp[105];
    } break;
    case 11: {
      sv[3] += dvdp[113];
      sv[5] += dvdp[115];
      sv[6] += dvdp[116];
    } break;
    case 12: {
      sv[4] += dvdp[124];
      sv[6] += dvdp[126];
    } break;
    case 13: {
      sv[6] += dvdp[136];
      sv[7] += dvdp[137];
      sv[8] += dvdp[138];
    } break;
    case 14: {
      sv[6] += dvdp[146];
      sv[7] += dvdp[147];
      sv[8] += dvdp[148];
    } break;
    case 15: {
      sv[6] += dvdp[156];
      sv[8] += dvdp[158];
      sv[9] += dvdp[159];
    } break;
    case 16: {
      sv[7] += dvdp[167];
      sv[9] += dvdp[169];
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
  for (is=0; is<10; is++) {
    if(mxIsNaN(sxdot_tmp[is])) sxdot_tmp[is] = 0.0;
  }

  return(0);
}


 int subfsx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(int Ns, realtype t, N_Vector x, N_Vector xdot, 
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
 {
   UserData data = (UserData) user_data;
   return fsx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(Ns, t, x, xdot, data->sensIndices[ip], sx, sxdot, user_data, tmp1, tmp2);
 };

 void csv_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data)
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
  for (is=0; is<10; is++) {
    sv[is] = 0.0;
  }
  sv[1] = dvdu[1]*su[(ip*5)+0]+dvdx[1]*sx_tmp[0]+dvdx[11]*sx_tmp[1]+dvdx[21]*sx_tmp[2];
  sv[2] = dvdu[2]*su[(ip*5)+0]+dvdx[12]*sx_tmp[1]+dvdx[22]*sx_tmp[2];
  sv[3] = dvdx[23]*sx_tmp[2]+dvdx[33]*sx_tmp[3]+dvdx[43]*sx_tmp[4]+dvdx[53]*sx_tmp[5];
  sv[4] = dvdx[34]*sx_tmp[3]+dvdx[44]*sx_tmp[4]+dvdx[54]*sx_tmp[5]+dvdx[64]*sx_tmp[6];
  sv[5] = dvdx[35]*sx_tmp[3]+dvdx[45]*sx_tmp[4]+dvdx[55]*sx_tmp[5];
  sv[6] = dvdx[56]*sx_tmp[5]+dvdx[66]*sx_tmp[6]+dvdx[76]*sx_tmp[7]+dvdx[86]*sx_tmp[8];
  sv[7] = dvdx[67]*sx_tmp[6]+dvdx[77]*sx_tmp[7]+dvdx[87]*sx_tmp[8]+dvdx[97]*sx_tmp[9];
  sv[8] = dvdx[68]*sx_tmp[6]+dvdx[78]*sx_tmp[7]+dvdx[88]*sx_tmp[8];
  sv[9] = dvdx[89]*sx_tmp[8]+dvdx[99]*sx_tmp[9];
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

    } break;
    case 5: {
      sv[1] += dvdp[51];
      sv[2] += dvdp[52];
    } break;
    case 6: {
      sv[1] += dvdp[61];
      sv[2] += dvdp[62];
    } break;
    case 7: {
      sv[2] += dvdp[72];
      sv[3] += dvdp[73];
    } break;
    case 8: {
      sv[3] += dvdp[83];
    } break;
    case 9: {
      sv[3] += dvdp[93];
      sv[4] += dvdp[94];
      sv[5] += dvdp[95];
    } break;
    case 10: {
      sv[3] += dvdp[103];
      sv[4] += dvdp[104];
      sv[5] += dvdp[105];
    } break;
    case 11: {
      sv[3] += dvdp[113];
      sv[5] += dvdp[115];
      sv[6] += dvdp[116];
    } break;
    case 12: {
      sv[4] += dvdp[124];
      sv[6] += dvdp[126];
    } break;
    case 13: {
      sv[6] += dvdp[136];
      sv[7] += dvdp[137];
      sv[8] += dvdp[138];
    } break;
    case 14: {
      sv[6] += dvdp[146];
      sv[7] += dvdp[147];
      sv[8] += dvdp[148];
    } break;
    case 15: {
      sv[6] += dvdp[156];
      sv[8] += dvdp[158];
      sv[9] += dvdp[159];
    } break;
    case 16: {
      sv[7] += dvdp[167];
      sv[9] += dvdp[169];
    } break;
  }

  return;
}


 void fsx0_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  switch (ip) {
    case 1: {
      sx0_tmp[1] = 1.0;
    } break;
  }

  return;
}


 void subfsx0_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(int ip, N_Vector sx0, void *user_data)
 {
   UserData data = (UserData) user_data;
   fsx0_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(data->sensIndices[ip], sx0, user_data);
 };

 void dfxdp0_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, double *dfxdp0, void *user_data)
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
  dfxdp0[11] = dvdx[11];
  dfxdp0[12] = dvdx[12];
  dfxdp0[51] = dvdp[51];
  dfxdp0[52] = dvdp[52];
  dfxdp0[61] = dvdp[61];
  dfxdp0[62] = dvdp[62];
  dfxdp0[72] = dvdp[72];
  dfxdp0[73] = dvdp[73];
  dfxdp0[83] = dvdp[83];
  dfxdp0[93] = dvdp[93];
  dfxdp0[94] = dvdp[94];
  dfxdp0[95] = dvdp[95];
  dfxdp0[103] = dvdp[103];
  dfxdp0[104] = dvdp[104];
  dfxdp0[105] = dvdp[105];
  dfxdp0[113] = dvdp[113];
  dfxdp0[115] = dvdp[115];
  dfxdp0[116] = dvdp[116];
  dfxdp0[124] = dvdp[124];
  dfxdp0[126] = dvdp[126];
  dfxdp0[136] = dvdp[136];
  dfxdp0[137] = dvdp[137];
  dfxdp0[138] = dvdp[138];
  dfxdp0[146] = dvdp[146];
  dfxdp0[147] = dvdp[147];
  dfxdp0[148] = dvdp[148];
  dfxdp0[156] = dvdp[156];
  dfxdp0[158] = dvdp[158];
  dfxdp0[159] = dvdp[159];
  dfxdp0[167] = dvdp[167];
  dfxdp0[169] = dvdp[169];
  for (is=0; is<170; is++) {
    if(mxIsNaN(dfxdp0[is])) dfxdp0[is] = 0.0;
  }

  return;
}


 void dfxdp_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(realtype t, N_Vector x, double *dfxdp, void *user_data)
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
  dfxdp[51] = dvdp[51];
  dfxdp[52] = dvdp[52];
  dfxdp[61] = dvdp[61];
  dfxdp[62] = dvdp[62];
  dfxdp[72] = dvdp[72];
  dfxdp[73] = dvdp[73];
  dfxdp[83] = dvdp[83];
  dfxdp[93] = dvdp[93];
  dfxdp[94] = dvdp[94];
  dfxdp[95] = dvdp[95];
  dfxdp[103] = dvdp[103];
  dfxdp[104] = dvdp[104];
  dfxdp[105] = dvdp[105];
  dfxdp[113] = dvdp[113];
  dfxdp[115] = dvdp[115];
  dfxdp[116] = dvdp[116];
  dfxdp[124] = dvdp[124];
  dfxdp[126] = dvdp[126];
  dfxdp[136] = dvdp[136];
  dfxdp[137] = dvdp[137];
  dfxdp[138] = dvdp[138];
  dfxdp[146] = dvdp[146];
  dfxdp[147] = dvdp[147];
  dfxdp[148] = dvdp[148];
  dfxdp[156] = dvdp[156];
  dfxdp[158] = dvdp[158];
  dfxdp[159] = dvdp[159];
  dfxdp[167] = dvdp[167];
  dfxdp[169] = dvdp[169];
  for (is=0; is<170; is++) {
    if(mxIsNaN(dfxdp[is])) dfxdp[is] = 0.0;
  }

  return;
}


 void fz_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x){
  z[nz*nt*iruns+it+nt*0] = (x[nx*nt*iruns+it+nt*3]+x[nx*nt*iruns+it+nt*5])*p[3];
  z[nz*nt*iruns+it+nt*1] = (x[nx*nt*iruns+it+nt*6]+x[nx*nt*iruns+it+nt*8])*p[2];
  z[nz*nt*iruns+it+nt*2] = p[4]*x[nx*nt*iruns+it+nt*9];

  return;
}


 void fsz_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx){
  int jp;
  for (jp=0; jp<np; jp++) {
      sz[it + nt*3*jp + nt*0] = p[3]*sx[it + nt*10*jp + nt*3]+p[3]*sx[it + nt*10*jp + nt*5];
      sz[it + nt*3*jp + nt*1] = p[2]*sx[it + nt*10*jp + nt*6]+p[2]*sx[it + nt*10*jp + nt*8];
      sz[it + nt*3*jp + nt*2] = p[4]*sx[it + nt*10*jp + nt*9];
  };

  sz[it+nt*7] += x[it+nt*6]+x[it+nt*8];
  sz[it+nt*9] += x[it+nt*3]+x[it+nt*5];
  sz[it+nt*14] += x[it+nt*9];

  return;
}


 void dfzdx_AktPathwayFujita_51F1A11B81F1FAE2328E5D9D4206AF54(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x){
      dfzdxs[nx*nt*iruns+it+nt*9] = p[3];
      dfzdxs[nx*nt*iruns+it+nt*15] = p[3];
      dfzdxs[nx*nt*iruns+it+nt*19] = p[2];
      dfzdxs[nx*nt*iruns+it+nt*25] = p[2];
      dfzdxs[nx*nt*iruns+it+nt*29] = p[4];


  return;
}


