#include "RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892.h"
#include "ELISA_Nigericin_formatted_WT_DF6FE2ED4623F9C78A4E869344080A31.h"
#include "ELISA_Nigericin_formatted_WT_AC1655417E166DD17CB2887C9716336C.h"
#include "ELISA_Nigericin_formatted_WT_51C280AA503F24CE4AEC3074E055D372.h"
#include "ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6.h"
#include "ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2.h"
#include "ELISA_Nigericin_formatted_WT_C027752051CE350E0D7725D6F30DB58B.h"
#include "ELISA_Nigericin_formatted_WT_B8208439898CC2FA44F23532772DECF1.h"

 int AR_CVodeInit(void *cvode_mem, N_Vector x, double t, int im, int ic){
  if((im==0) & (ic==0)) return CVodeInit(cvode_mem, fx_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892, RCONST(t), x);
  return(-1);
}

 void fx(realtype t, N_Vector x, double *xdot, void *user_data, int im, int ic){
  if((im==0) & (ic==0)) fxdouble_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, x, xdot, user_data);
}

 void fx0(N_Vector x0, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fx0_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(x0, data);
}

 int AR_CVDlsSetDenseJacFn(void *cvode_mem, int im, int ic, int setSparse){
  if((im==0) & (ic==0) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892);
 
 }else if((im==0) & (ic==0) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892);

}
  return(-1);
}

 void getdfxdx(int im, int ic, realtype t, N_Vector x, realtype *J, void *user_data){
  if((im==0) & (ic==0)) { dfxdx_out_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, x, J, user_data); return; }

}

 void fsx0(int is, N_Vector sx_is, void *user_data, int im, int ic, int sensitivitySubset) {
  UserData data = (UserData) user_data;
  if ( sensitivitySubset == 0 ) {
    if((im==0) & (ic==0)) fsx0_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(is, sx_is, data);
  } else {
    if((im==0) & (ic==0)) subfsx0_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(is, sx_is, data);
  }
}

 void csv(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) csv_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, x, ip, sx, data);
}

 int AR_CVodeSensInit1(void *cvode_mem, int nps, int sensi_meth, int sensirhs, N_Vector *sx, int im, int ic, int sensitivitySubset){
  if (sensirhs == 1) {
    if (sensitivitySubset == 0) {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892, sx);
    } else {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892, sx);
  }
  } else {
    if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
  }
  return(-1);
}

 void fu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fu_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(data, t);
}

 void fsu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fsu_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(data, t);
}

 void fv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fv_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, x, data);
}

 void fsv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
	dvdp_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, x, data);
	dvdu_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, x, data);
	dvdx_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, x, data);
}
}

 void dfxdp0(void *user_data, double t, N_Vector x, double *dfxdp0, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp0_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, x, dfxdp0, data);
}

 void dfxdp(void *user_data, double t, N_Vector x, double *dfxdp, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, x, dfxdp, data);
}

void fz(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) fz_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
}

void dfzdx(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdx, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) dfzdx_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
}

void fsz(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx, int im, int ic){
  if((im==0) & (ic==0)) fsz_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, nt, it, np, sz, p, u, x, z, su, sx);
}

 void fy(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fy_ELISA_Nigericin_formatted_WT_DF6FE2ED4623F9C78A4E869344080A31(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==1)) fy_ELISA_Nigericin_formatted_WT_AC1655417E166DD17CB2887C9716336C(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==2)) fy_ELISA_Nigericin_formatted_WT_51C280AA503F24CE4AEC3074E055D372(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==3)) fy_ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==4)) fy_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==5)) fy_ELISA_Nigericin_formatted_WT_C027752051CE350E0D7725D6F30DB58B(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==6)) fy_ELISA_Nigericin_formatted_WT_B8208439898CC2FA44F23532772DECF1(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
}

 void fy_scale(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx, int im, int id){
  if((im==0) & (id==0)) fy_scale_ELISA_Nigericin_formatted_WT_DF6FE2ED4623F9C78A4E869344080A31(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==1)) fy_scale_ELISA_Nigericin_formatted_WT_AC1655417E166DD17CB2887C9716336C(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==2)) fy_scale_ELISA_Nigericin_formatted_WT_51C280AA503F24CE4AEC3074E055D372(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==3)) fy_scale_ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==4)) fy_scale_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==5)) fy_scale_ELISA_Nigericin_formatted_WT_C027752051CE350E0D7725D6F30DB58B(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==6)) fy_scale_ELISA_Nigericin_formatted_WT_B8208439898CC2FA44F23532772DECF1(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
}

 void fystd(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fystd_ELISA_Nigericin_formatted_WT_DF6FE2ED4623F9C78A4E869344080A31(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==1)) fystd_ELISA_Nigericin_formatted_WT_AC1655417E166DD17CB2887C9716336C(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==2)) fystd_ELISA_Nigericin_formatted_WT_51C280AA503F24CE4AEC3074E055D372(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==3)) fystd_ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==4)) fystd_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==5)) fystd_ELISA_Nigericin_formatted_WT_C027752051CE350E0D7725D6F30DB58B(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==6)) fystd_ELISA_Nigericin_formatted_WT_B8208439898CC2FA44F23532772DECF1(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
}

 void fsy(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsy_ELISA_Nigericin_formatted_WT_DF6FE2ED4623F9C78A4E869344080A31(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==1)) fsy_ELISA_Nigericin_formatted_WT_AC1655417E166DD17CB2887C9716336C(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==2)) fsy_ELISA_Nigericin_formatted_WT_51C280AA503F24CE4AEC3074E055D372(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==3)) fsy_ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==4)) fsy_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==5)) fsy_ELISA_Nigericin_formatted_WT_C027752051CE350E0D7725D6F30DB58B(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==6)) fsy_ELISA_Nigericin_formatted_WT_B8208439898CC2FA44F23532772DECF1(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
}

 void fsystd(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsystd_ELISA_Nigericin_formatted_WT_DF6FE2ED4623F9C78A4E869344080A31(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==1)) fsystd_ELISA_Nigericin_formatted_WT_AC1655417E166DD17CB2887C9716336C(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==2)) fsystd_ELISA_Nigericin_formatted_WT_51C280AA503F24CE4AEC3074E055D372(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==3)) fsystd_ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==4)) fsystd_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==5)) fsystd_ELISA_Nigericin_formatted_WT_C027752051CE350E0D7725D6F30DB58B(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==6)) fsystd_ELISA_Nigericin_formatted_WT_B8208439898CC2FA44F23532772DECF1(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
}

/* for arSSACalc.c */

 void fvSSA(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
    fu_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(data, t);
    fv_RTF_DoseDependent_WT_F298A94C5A6864FEB754E83604898892(t, x, data);
  }
}

