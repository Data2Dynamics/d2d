#include "ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3.h"
#include "ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400.h"

 int AR_CVodeInit(void *cvode_mem, N_Vector x, double t, int im, int ic){
  if((im==0) & (ic==0)) return CVodeInit(cvode_mem, fx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3, RCONST(t), x);
  return(-1);
}

 void fx(realtype t, N_Vector x, double *xdot, void *user_data, int im, int ic){
  if((im==0) & (ic==0)) fxdouble_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, xdot, user_data);
}

 void fx0(N_Vector x0, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fx0_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(x0, data);
}

 int AR_CVDlsSetDenseJacFn(void *cvode_mem, int im, int ic, int setSparse){
  if((im==0) & (ic==0) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3);
 
 }else if((im==0) & (ic==0) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3);

}
  return(-1);
}

 void getdfxdx(int im, int ic, realtype t, N_Vector x, realtype *J, void *user_data){
  if((im==0) & (ic==0)) { dfxdx_out_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, J, user_data); return; }

}

 void fsx0(int is, N_Vector sx_is, void *user_data, int im, int ic, int sensitivitySubset) {
  UserData data = (UserData) user_data;
  if ( sensitivitySubset == 0 ) {
    if((im==0) & (ic==0)) fsx0_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(is, sx_is, data);
  } else {
    if((im==0) & (ic==0)) subfsx0_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(is, sx_is, data);
  }
}

 void csv(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) csv_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, ip, sx, data);
}

 int AR_CVodeSensInit1(void *cvode_mem, int nps, int sensi_meth, int sensirhs, N_Vector *sx, int im, int ic, int sensitivitySubset){
  if (sensirhs == 1) {
    if (sensitivitySubset == 0) {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3, sx);
    } else {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3, sx);
  }
  } else {
    if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
  }
  return(-1);
}

 void fu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(data, t);
}

 void fsu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fsu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(data, t);
}

 void fv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fv_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
}

 void fsv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
	dvdp_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
	dvdu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
	dvdx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
}
}

 void dfxdp0(void *user_data, double t, N_Vector x, double *dfxdp0, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp0_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, dfxdp0, data);
}

 void dfxdp(void *user_data, double t, N_Vector x, double *dfxdp, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, dfxdp, data);
}

void fz(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) fz_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
}

void dfzdx(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdx, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) dfzdx_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
}

void fsz(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx, int im, int ic){
  if((im==0) & (ic==0)) fsz_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, nt, it, np, sz, p, u, x, z, su, sx);
}

 void fy(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fy_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
}

 void fy_scale(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx, int im, int id){
  if((im==0) & (id==0)) fy_scale_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
}

 void fystd(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fystd_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
}

 void fsy(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsy_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
}

 void fsystd(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsystd_ABC_data_2_44E88DC0B59393FDAE87A7F44E7C6400(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
}

/* for arSSACalc.c */

 void fvSSA(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
    fu_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(data, t);
    fv_ABC_model_C53B93EFA95D48EDC27A14D83EEA27A3(t, x, data);
  }
}

