#include "model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0.h"
#include "model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE.h"
#include "model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277.h"
#include "model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD.h"
#include "model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11.h"
#include "model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A.h"
#include "experimentaldata1_EED7A0B0CCCAD902928521279111B00B.h"
#include "experimentaldata2_436D55A61FF8D552D4D75D610B44C192.h"
#include "experimentaldata3_B25559DB36F7BC324746F0FDABDBE57D.h"
#include "experimentaldata4_9C5AE63E55BAF34680B30E7FE8B5D589.h"
#include "experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5.h"
#include "experimentaldata6_17E8ABD116F4EB7F6F5147246F5E3109.h"

 int AR_CVodeInit(void *cvode_mem, N_Vector x, double t, int im, int ic){
  if((im==0) & (ic==0)) return CVodeInit(cvode_mem, fx_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0, RCONST(t), x);
  if((im==0) & (ic==1)) return CVodeInit(cvode_mem, fx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE, RCONST(t), x);
  if((im==0) & (ic==2)) return CVodeInit(cvode_mem, fx_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277, RCONST(t), x);
  if((im==0) & (ic==3)) return CVodeInit(cvode_mem, fx_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD, RCONST(t), x);
  if((im==0) & (ic==4)) return CVodeInit(cvode_mem, fx_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11, RCONST(t), x);
  if((im==0) & (ic==5)) return CVodeInit(cvode_mem, fx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A, RCONST(t), x);
  return(-1);
}

 void fx(realtype t, N_Vector x, double *xdot, void *user_data, int im, int ic){
  if((im==0) & (ic==0)) fxdouble_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, x, xdot, user_data);
  if((im==0) & (ic==1)) fxdouble_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, xdot, user_data);
  if((im==0) & (ic==2)) fxdouble_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, x, xdot, user_data);
  if((im==0) & (ic==3)) fxdouble_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, x, xdot, user_data);
  if((im==0) & (ic==4)) fxdouble_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, x, xdot, user_data);
  if((im==0) & (ic==5)) fxdouble_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, x, xdot, user_data);
}

 void fx0(N_Vector x0, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fx0_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(x0, data);
  if((im==0) & (ic==1)) fx0_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(x0, data);
  if((im==0) & (ic==2)) fx0_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(x0, data);
  if((im==0) & (ic==3)) fx0_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(x0, data);
  if((im==0) & (ic==4)) fx0_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(x0, data);
  if((im==0) & (ic==5)) fx0_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(x0, data);
}

 int AR_CVDlsSetDenseJacFn(void *cvode_mem, int im, int ic, int setSparse){
  if((im==0) & (ic==0) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0);
 
 }else if((im==0) & (ic==0) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0);

}
  if((im==0) & (ic==1) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE);
 
 }else if((im==0) & (ic==1) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE);

}
  if((im==0) & (ic==2) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277);
 
 }else if((im==0) & (ic==2) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277);

}
  if((im==0) & (ic==3) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD);
 
 }else if((im==0) & (ic==3) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD);

}
  if((im==0) & (ic==4) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11);
 
 }else if((im==0) & (ic==4) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11);

}
  if((im==0) & (ic==5) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A);
 
 }else if((im==0) & (ic==5) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A);

}
  return(-1);
}

 void getdfxdx(int im, int ic, realtype t, N_Vector x, realtype *J, void *user_data){
  if((im==0) & (ic==0)) { dfxdx_out_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, x, J, user_data); return; }
  if((im==0) & (ic==1)) { dfxdx_out_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, J, user_data); return; }
  if((im==0) & (ic==2)) { dfxdx_out_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, x, J, user_data); return; }
  if((im==0) & (ic==3)) { dfxdx_out_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, x, J, user_data); return; }
  if((im==0) & (ic==4)) { dfxdx_out_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, x, J, user_data); return; }
  if((im==0) & (ic==5)) { dfxdx_out_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, x, J, user_data); return; }

}

 void fsx0(int is, N_Vector sx_is, void *user_data, int im, int ic, int sensitivitySubset) {
  UserData data = (UserData) user_data;
  if ( sensitivitySubset == 0 ) {
    if((im==0) & (ic==0)) fsx0_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(is, sx_is, data);
    if((im==0) & (ic==1)) fsx0_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(is, sx_is, data);
    if((im==0) & (ic==2)) fsx0_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(is, sx_is, data);
    if((im==0) & (ic==3)) fsx0_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(is, sx_is, data);
    if((im==0) & (ic==4)) fsx0_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(is, sx_is, data);
    if((im==0) & (ic==5)) fsx0_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(is, sx_is, data);
  } else {
    if((im==0) & (ic==0)) subfsx0_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(is, sx_is, data);
    if((im==0) & (ic==1)) subfsx0_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(is, sx_is, data);
    if((im==0) & (ic==2)) subfsx0_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(is, sx_is, data);
    if((im==0) & (ic==3)) subfsx0_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(is, sx_is, data);
    if((im==0) & (ic==4)) subfsx0_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(is, sx_is, data);
    if((im==0) & (ic==5)) subfsx0_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(is, sx_is, data);
  }
}

 void csv(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) csv_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, x, ip, sx, data);
  if((im==0) & (ic==1)) csv_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, ip, sx, data);
  if((im==0) & (ic==2)) csv_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, x, ip, sx, data);
  if((im==0) & (ic==3)) csv_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, x, ip, sx, data);
  if((im==0) & (ic==4)) csv_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, x, ip, sx, data);
  if((im==0) & (ic==5)) csv_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, x, ip, sx, data);
}

 int AR_CVodeSensInit1(void *cvode_mem, int nps, int sensi_meth, int sensirhs, N_Vector *sx, int im, int ic, int sensitivitySubset){
  if (sensirhs == 1) {
    if (sensitivitySubset == 0) {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0, sx);
      if((im==0) & (ic==1)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE, sx);
      if((im==0) & (ic==2)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277, sx);
      if((im==0) & (ic==3)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD, sx);
      if((im==0) & (ic==4)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11, sx);
      if((im==0) & (ic==5)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A, sx);
    } else {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0, sx);
      if((im==0) & (ic==1)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE, sx);
      if((im==0) & (ic==2)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277, sx);
      if((im==0) & (ic==3)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD, sx);
      if((im==0) & (ic==4)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11, sx);
      if((im==0) & (ic==5)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A, sx);
  }
  } else {
    if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
    if((im==0) & (ic==1)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
    if((im==0) & (ic==2)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
    if((im==0) & (ic==3)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
    if((im==0) & (ic==4)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
    if((im==0) & (ic==5)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
  }
  return(-1);
}

 void fu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fu_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(data, t);
  if((im==0) & (ic==1)) fu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(data, t);
  if((im==0) & (ic==2)) fu_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(data, t);
  if((im==0) & (ic==3)) fu_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(data, t);
  if((im==0) & (ic==4)) fu_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(data, t);
  if((im==0) & (ic==5)) fu_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(data, t);
}

 void fsu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fsu_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(data, t);
  if((im==0) & (ic==1)) fsu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(data, t);
  if((im==0) & (ic==2)) fsu_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(data, t);
  if((im==0) & (ic==3)) fsu_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(data, t);
  if((im==0) & (ic==4)) fsu_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(data, t);
  if((im==0) & (ic==5)) fsu_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(data, t);
}

 void fv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fv_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, x, data);
  if((im==0) & (ic==1)) fv_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
  if((im==0) & (ic==2)) fv_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, x, data);
  if((im==0) & (ic==3)) fv_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, x, data);
  if((im==0) & (ic==4)) fv_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, x, data);
  if((im==0) & (ic==5)) fv_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, x, data);
}

 void fsv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
	dvdp_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, x, data);
	dvdu_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, x, data);
	dvdx_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, x, data);
}
  if((im==0) & (ic==1)) {
	dvdp_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
	dvdu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
	dvdx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
}
  if((im==0) & (ic==2)) {
	dvdp_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, x, data);
	dvdu_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, x, data);
	dvdx_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, x, data);
}
  if((im==0) & (ic==3)) {
	dvdp_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, x, data);
	dvdu_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, x, data);
	dvdx_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, x, data);
}
  if((im==0) & (ic==4)) {
	dvdp_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, x, data);
	dvdu_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, x, data);
	dvdx_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, x, data);
}
  if((im==0) & (ic==5)) {
	dvdp_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, x, data);
	dvdu_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, x, data);
	dvdx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, x, data);
}
}

 void dfxdp0(void *user_data, double t, N_Vector x, double *dfxdp0, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp0_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, x, dfxdp0, data);
  if((im==0) & (ic==1)) dfxdp0_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, dfxdp0, data);
  if((im==0) & (ic==2)) dfxdp0_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, x, dfxdp0, data);
  if((im==0) & (ic==3)) dfxdp0_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, x, dfxdp0, data);
  if((im==0) & (ic==4)) dfxdp0_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, x, dfxdp0, data);
  if((im==0) & (ic==5)) dfxdp0_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, x, dfxdp0, data);
}

 void dfxdp(void *user_data, double t, N_Vector x, double *dfxdp, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, x, dfxdp, data);
  if((im==0) & (ic==1)) dfxdp_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, dfxdp, data);
  if((im==0) & (ic==2)) dfxdp_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, x, dfxdp, data);
  if((im==0) & (ic==3)) dfxdp_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, x, dfxdp, data);
  if((im==0) & (ic==4)) dfxdp_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, x, dfxdp, data);
  if((im==0) & (ic==5)) dfxdp_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, x, dfxdp, data);
}

void fz(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) fz_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==1)) fz_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==2)) fz_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==3)) fz_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==4)) fz_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==5)) fz_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
}

void dfzdx(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdx, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) dfzdx_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==1)) dfzdx_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==2)) dfzdx_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==3)) dfzdx_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==4)) dfzdx_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==5)) dfzdx_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
}

void fsz(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx, int im, int ic){
  if((im==0) & (ic==0)) fsz_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==1)) fsz_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==2)) fsz_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==3)) fsz_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==4)) fsz_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==5)) fsz_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, nt, it, np, sz, p, u, x, z, su, sx);
}

 void fy(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fy_experimentaldata1_EED7A0B0CCCAD902928521279111B00B(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==1)) fy_experimentaldata2_436D55A61FF8D552D4D75D610B44C192(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==2)) fy_experimentaldata3_B25559DB36F7BC324746F0FDABDBE57D(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==3)) fy_experimentaldata4_9C5AE63E55BAF34680B30E7FE8B5D589(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==4)) fy_experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==5)) fy_experimentaldata6_17E8ABD116F4EB7F6F5147246F5E3109(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
}

 void fy_scale(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx, int im, int id){
  if((im==0) & (id==0)) fy_scale_experimentaldata1_EED7A0B0CCCAD902928521279111B00B(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==1)) fy_scale_experimentaldata2_436D55A61FF8D552D4D75D610B44C192(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==2)) fy_scale_experimentaldata3_B25559DB36F7BC324746F0FDABDBE57D(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==3)) fy_scale_experimentaldata4_9C5AE63E55BAF34680B30E7FE8B5D589(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==4)) fy_scale_experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==5)) fy_scale_experimentaldata6_17E8ABD116F4EB7F6F5147246F5E3109(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
}

 void fystd(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fystd_experimentaldata1_EED7A0B0CCCAD902928521279111B00B(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==1)) fystd_experimentaldata2_436D55A61FF8D552D4D75D610B44C192(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==2)) fystd_experimentaldata3_B25559DB36F7BC324746F0FDABDBE57D(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==3)) fystd_experimentaldata4_9C5AE63E55BAF34680B30E7FE8B5D589(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==4)) fystd_experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==5)) fystd_experimentaldata6_17E8ABD116F4EB7F6F5147246F5E3109(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
}

 void fsy(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsy_experimentaldata1_EED7A0B0CCCAD902928521279111B00B(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==1)) fsy_experimentaldata2_436D55A61FF8D552D4D75D610B44C192(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==2)) fsy_experimentaldata3_B25559DB36F7BC324746F0FDABDBE57D(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==3)) fsy_experimentaldata4_9C5AE63E55BAF34680B30E7FE8B5D589(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==4)) fsy_experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==5)) fsy_experimentaldata6_17E8ABD116F4EB7F6F5147246F5E3109(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
}

 void fsystd(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsystd_experimentaldata1_EED7A0B0CCCAD902928521279111B00B(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==1)) fsystd_experimentaldata2_436D55A61FF8D552D4D75D610B44C192(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==2)) fsystd_experimentaldata3_B25559DB36F7BC324746F0FDABDBE57D(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==3)) fsystd_experimentaldata4_9C5AE63E55BAF34680B30E7FE8B5D589(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==4)) fsystd_experimentaldata5_9047403079DCEEF0057FC05B5A3AC8C5(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==5)) fsystd_experimentaldata6_17E8ABD116F4EB7F6F5147246F5E3109(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
}

/* for arSSACalc.c */

 void fvSSA(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
    fu_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(data, t);
    fv_model_AktPathwayFujita_CACD23CE1D9670C0C9EFAAA0D88341E0(t, x, data);
  }
  if((im==0) & (ic==1)) {
    fu_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(data, t);
    fv_model_AktPathwayFujita_727DD824F350367B8ABC8B349AADFABE(t, x, data);
  }
  if((im==0) & (ic==2)) {
    fu_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(data, t);
    fv_model_AktPathwayFujita_19691A74565C8A7A565C743C011C7277(t, x, data);
  }
  if((im==0) & (ic==3)) {
    fu_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(data, t);
    fv_model_AktPathwayFujita_7DC996A9736DCFC3E6C787B3DF25DFFD(t, x, data);
  }
  if((im==0) & (ic==4)) {
    fu_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(data, t);
    fv_model_AktPathwayFujita_F5FDCDF25EEC60D7DB4C9798044F3D11(t, x, data);
  }
  if((im==0) & (ic==5)) {
    fu_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(data, t);
    fv_model_AktPathwayFujita_8518F4FA0CE00579B4676CC94FA4356A(t, x, data);
  }
}

