#include "ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61.h"
#include "ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98.h"
#include "ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE.h"
#include "ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723.h"
#include "experimentaldata1_F2E72D0B72CC4E6D278D13BD39066899.h"
#include "experimentaldata2_831558B7E3ED5C226179982B73A46012.h"
#include "experimentaldata3_28DF2F5C2299B091259C19607121A9E9.h"
#include "experimentaldata4_ACB2471524E3BC799814F6785E0C4CC0.h"

 int AR_CVodeInit(void *cvode_mem, N_Vector x, double t, int im, int ic){
  if((im==0) & (ic==0)) return CVodeInit(cvode_mem, fx_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61, RCONST(t), x);
  if((im==0) & (ic==1)) return CVodeInit(cvode_mem, fx_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98, RCONST(t), x);
  if((im==0) & (ic==2)) return CVodeInit(cvode_mem, fx_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE, RCONST(t), x);
  if((im==0) & (ic==3)) return CVodeInit(cvode_mem, fx_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723, RCONST(t), x);
  return(-1);
}

 void fx(realtype t, N_Vector x, double *xdot, void *user_data, int im, int ic){
  if((im==0) & (ic==0)) fxdouble_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, x, xdot, user_data);
  if((im==0) & (ic==1)) fxdouble_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, x, xdot, user_data);
  if((im==0) & (ic==2)) fxdouble_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, x, xdot, user_data);
  if((im==0) & (ic==3)) fxdouble_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, x, xdot, user_data);
}

 void fx0(N_Vector x0, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fx0_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(x0, data);
  if((im==0) & (ic==1)) fx0_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(x0, data);
  if((im==0) & (ic==2)) fx0_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(x0, data);
  if((im==0) & (ic==3)) fx0_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(x0, data);
}

 int AR_CVDlsSetDenseJacFn(void *cvode_mem, int im, int ic, int setSparse){
  if((im==0) & (ic==0) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61);
 
 }else if((im==0) & (ic==0) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61);

}
  if((im==0) & (ic==1) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98);
 
 }else if((im==0) & (ic==1) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98);

}
  if((im==0) & (ic==2) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE);
 
 }else if((im==0) & (ic==2) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE);

}
  if((im==0) & (ic==3) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723);
 
 }else if((im==0) & (ic==3) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723);

}
  return(-1);
}

 void getdfxdx(int im, int ic, realtype t, N_Vector x, realtype *J, void *user_data){
  if((im==0) & (ic==0)) { dfxdx_out_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, x, J, user_data); return; }
  if((im==0) & (ic==1)) { dfxdx_out_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, x, J, user_data); return; }
  if((im==0) & (ic==2)) { dfxdx_out_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, x, J, user_data); return; }
  if((im==0) & (ic==3)) { dfxdx_out_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, x, J, user_data); return; }

}

 void fsx0(int is, N_Vector sx_is, void *user_data, int im, int ic, int sensitivitySubset) {
  UserData data = (UserData) user_data;
  if ( sensitivitySubset == 0 ) {
    if((im==0) & (ic==0)) fsx0_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(is, sx_is, data);
    if((im==0) & (ic==1)) fsx0_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(is, sx_is, data);
    if((im==0) & (ic==2)) fsx0_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(is, sx_is, data);
    if((im==0) & (ic==3)) fsx0_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(is, sx_is, data);
  } else {
    if((im==0) & (ic==0)) subfsx0_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(is, sx_is, data);
    if((im==0) & (ic==1)) subfsx0_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(is, sx_is, data);
    if((im==0) & (ic==2)) subfsx0_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(is, sx_is, data);
    if((im==0) & (ic==3)) subfsx0_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(is, sx_is, data);
  }
}

 void csv(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) csv_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, x, ip, sx, data);
  if((im==0) & (ic==1)) csv_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, x, ip, sx, data);
  if((im==0) & (ic==2)) csv_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, x, ip, sx, data);
  if((im==0) & (ic==3)) csv_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, x, ip, sx, data);
}

 int AR_CVodeSensInit1(void *cvode_mem, int nps, int sensi_meth, int sensirhs, N_Vector *sx, int im, int ic, int sensitivitySubset){
  if (sensirhs == 1) {
    if (sensitivitySubset == 0) {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61, sx);
      if((im==0) & (ic==1)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98, sx);
      if((im==0) & (ic==2)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE, sx);
      if((im==0) & (ic==3)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723, sx);
    } else {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61, sx);
      if((im==0) & (ic==1)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98, sx);
      if((im==0) & (ic==2)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE, sx);
      if((im==0) & (ic==3)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723, sx);
  }
  } else {
    if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
    if((im==0) & (ic==1)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
    if((im==0) & (ic==2)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
    if((im==0) & (ic==3)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
  }
  return(-1);
}

 void fu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fu_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(data, t);
  if((im==0) & (ic==1)) fu_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(data, t);
  if((im==0) & (ic==2)) fu_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(data, t);
  if((im==0) & (ic==3)) fu_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(data, t);
}

 void fsu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fsu_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(data, t);
  if((im==0) & (ic==1)) fsu_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(data, t);
  if((im==0) & (ic==2)) fsu_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(data, t);
  if((im==0) & (ic==3)) fsu_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(data, t);
}

 void fv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fv_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, x, data);
  if((im==0) & (ic==1)) fv_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, x, data);
  if((im==0) & (ic==2)) fv_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, x, data);
  if((im==0) & (ic==3)) fv_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, x, data);
}

 void fsv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
	dvdp_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, x, data);
	dvdu_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, x, data);
	dvdx_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, x, data);
}
  if((im==0) & (ic==1)) {
	dvdp_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, x, data);
	dvdu_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, x, data);
	dvdx_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, x, data);
}
  if((im==0) & (ic==2)) {
	dvdp_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, x, data);
	dvdu_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, x, data);
	dvdx_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, x, data);
}
  if((im==0) & (ic==3)) {
	dvdp_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, x, data);
	dvdu_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, x, data);
	dvdx_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, x, data);
}
}

 void dfxdp0(void *user_data, double t, N_Vector x, double *dfxdp0, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp0_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, x, dfxdp0, data);
  if((im==0) & (ic==1)) dfxdp0_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, x, dfxdp0, data);
  if((im==0) & (ic==2)) dfxdp0_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, x, dfxdp0, data);
  if((im==0) & (ic==3)) dfxdp0_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, x, dfxdp0, data);
}

 void dfxdp(void *user_data, double t, N_Vector x, double *dfxdp, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, x, dfxdp, data);
  if((im==0) & (ic==1)) dfxdp_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, x, dfxdp, data);
  if((im==0) & (ic==2)) dfxdp_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, x, dfxdp, data);
  if((im==0) & (ic==3)) dfxdp_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, x, dfxdp, data);
}

void fz(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) fz_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==1)) fz_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==2)) fz_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==3)) fz_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
}

void dfzdx(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdx, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) dfzdx_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==1)) dfzdx_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==2)) dfzdx_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==3)) dfzdx_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
}

void fsz(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx, int im, int ic){
  if((im==0) & (ic==0)) fsz_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==1)) fsz_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==2)) fsz_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==3)) fsz_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, nt, it, np, sz, p, u, x, z, su, sx);
}

 void fy(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fy_experimentaldata1_F2E72D0B72CC4E6D278D13BD39066899(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==1)) fy_experimentaldata2_831558B7E3ED5C226179982B73A46012(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==2)) fy_experimentaldata3_28DF2F5C2299B091259C19607121A9E9(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==3)) fy_experimentaldata4_ACB2471524E3BC799814F6785E0C4CC0(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
}

 void fy_scale(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx, int im, int id){
  if((im==0) & (id==0)) fy_scale_experimentaldata1_F2E72D0B72CC4E6D278D13BD39066899(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==1)) fy_scale_experimentaldata2_831558B7E3ED5C226179982B73A46012(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==2)) fy_scale_experimentaldata3_28DF2F5C2299B091259C19607121A9E9(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==3)) fy_scale_experimentaldata4_ACB2471524E3BC799814F6785E0C4CC0(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
}

 void fystd(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fystd_experimentaldata1_F2E72D0B72CC4E6D278D13BD39066899(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==1)) fystd_experimentaldata2_831558B7E3ED5C226179982B73A46012(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==2)) fystd_experimentaldata3_28DF2F5C2299B091259C19607121A9E9(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==3)) fystd_experimentaldata4_ACB2471524E3BC799814F6785E0C4CC0(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
}

 void fsy(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsy_experimentaldata1_F2E72D0B72CC4E6D278D13BD39066899(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==1)) fsy_experimentaldata2_831558B7E3ED5C226179982B73A46012(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==2)) fsy_experimentaldata3_28DF2F5C2299B091259C19607121A9E9(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==3)) fsy_experimentaldata4_ACB2471524E3BC799814F6785E0C4CC0(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
}

 void fsystd(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsystd_experimentaldata1_F2E72D0B72CC4E6D278D13BD39066899(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==1)) fsystd_experimentaldata2_831558B7E3ED5C226179982B73A46012(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==2)) fsystd_experimentaldata3_28DF2F5C2299B091259C19607121A9E9(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==3)) fsystd_experimentaldata4_ACB2471524E3BC799814F6785E0C4CC0(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
}

/* for arSSACalc.c */

 void fvSSA(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
    fu_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(data, t);
    fv_ChenMSB2009_39A0DC04B59A44371D5CD48EE360BB61(t, x, data);
  }
  if((im==0) & (ic==1)) {
    fu_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(data, t);
    fv_ChenMSB2009_C87A636D5C1DC992F2A30307B2AA8F98(t, x, data);
  }
  if((im==0) & (ic==2)) {
    fu_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(data, t);
    fv_ChenMSB2009_420810DDAC6219DC1B2ADE86FCE31FCE(t, x, data);
  }
  if((im==0) & (ic==3)) {
    fu_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(data, t);
    fv_ChenMSB2009_2DB103E32FF3D2CE4B752CC531A1B723(t, x, data);
  }
}

