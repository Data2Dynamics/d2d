#include "AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7.h"
#include "AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB.h"
#include "AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6.h"
#include "AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B.h"
#include "AktPathwayFujita_682D2488568EE2F473172B1526A835D5.h"
#include "AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB.h"
#include "experimentaldata1_258544EB6B3AA8CBD3145AB82F974965.h"
#include "experimentaldata2_D90EAD8994BDEFB35ABE1905F22F8E93.h"
#include "experimentaldata3_C230EC340F7FFF25812B99F99870D51A.h"
#include "experimentaldata4_5F66C2E0235D12505358B893E4F722CF.h"
#include "experimentaldata5_5B4FCB1D8192A23C5B469A3F75AC4E18.h"
#include "experimentaldata6_29AF108A47A7B9DE3D5BF91CBED9712A.h"

 int AR_CVodeInit(void *cvode_mem, N_Vector x, double t, int im, int ic){
  if((im==0) & (ic==0)) return CVodeInit(cvode_mem, fx_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7, RCONST(t), x);
  if((im==0) & (ic==1)) return CVodeInit(cvode_mem, fx_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB, RCONST(t), x);
  if((im==0) & (ic==2)) return CVodeInit(cvode_mem, fx_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6, RCONST(t), x);
  if((im==0) & (ic==3)) return CVodeInit(cvode_mem, fx_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B, RCONST(t), x);
  if((im==0) & (ic==4)) return CVodeInit(cvode_mem, fx_AktPathwayFujita_682D2488568EE2F473172B1526A835D5, RCONST(t), x);
  if((im==0) & (ic==5)) return CVodeInit(cvode_mem, fx_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB, RCONST(t), x);
  return(-1);
}

 void fx(realtype t, N_Vector x, double *xdot, void *user_data, int im, int ic){
  if((im==0) & (ic==0)) fxdouble_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, x, xdot, user_data);
  if((im==0) & (ic==1)) fxdouble_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, x, xdot, user_data);
  if((im==0) & (ic==2)) fxdouble_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, x, xdot, user_data);
  if((im==0) & (ic==3)) fxdouble_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, x, xdot, user_data);
  if((im==0) & (ic==4)) fxdouble_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, x, xdot, user_data);
  if((im==0) & (ic==5)) fxdouble_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, x, xdot, user_data);
}

 void fx0(N_Vector x0, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fx0_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(x0, data);
  if((im==0) & (ic==1)) fx0_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(x0, data);
  if((im==0) & (ic==2)) fx0_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(x0, data);
  if((im==0) & (ic==3)) fx0_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(x0, data);
  if((im==0) & (ic==4)) fx0_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(x0, data);
  if((im==0) & (ic==5)) fx0_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(x0, data);
}

 int AR_CVDlsSetDenseJacFn(void *cvode_mem, int im, int ic, int setSparse){
  if((im==0) & (ic==0) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7);
 
 }else if((im==0) & (ic==0) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7);

}
  if((im==0) & (ic==1) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB);
 
 }else if((im==0) & (ic==1) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB);

}
  if((im==0) & (ic==2) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6);
 
 }else if((im==0) & (ic==2) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6);

}
  if((im==0) & (ic==3) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B);
 
 }else if((im==0) & (ic==3) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B);

}
  if((im==0) & (ic==4) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_AktPathwayFujita_682D2488568EE2F473172B1526A835D5);
 
 }else if((im==0) & (ic==4) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_AktPathwayFujita_682D2488568EE2F473172B1526A835D5);

}
  if((im==0) & (ic==5) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB);
 
 }else if((im==0) & (ic==5) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB);

}
  return(-1);
}

 void getdfxdx(int im, int ic, realtype t, N_Vector x, realtype *J, void *user_data){
  if((im==0) & (ic==0)) { dfxdx_out_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, x, J, user_data); return; }
  if((im==0) & (ic==1)) { dfxdx_out_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, x, J, user_data); return; }
  if((im==0) & (ic==2)) { dfxdx_out_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, x, J, user_data); return; }
  if((im==0) & (ic==3)) { dfxdx_out_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, x, J, user_data); return; }
  if((im==0) & (ic==4)) { dfxdx_out_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, x, J, user_data); return; }
  if((im==0) & (ic==5)) { dfxdx_out_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, x, J, user_data); return; }

}

 void fsx0(int is, N_Vector sx_is, void *user_data, int im, int ic, int sensitivitySubset) {
  UserData data = (UserData) user_data;
  if ( sensitivitySubset == 0 ) {
    if((im==0) & (ic==0)) fsx0_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(is, sx_is, data);
    if((im==0) & (ic==1)) fsx0_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(is, sx_is, data);
    if((im==0) & (ic==2)) fsx0_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(is, sx_is, data);
    if((im==0) & (ic==3)) fsx0_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(is, sx_is, data);
    if((im==0) & (ic==4)) fsx0_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(is, sx_is, data);
    if((im==0) & (ic==5)) fsx0_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(is, sx_is, data);
  } else {
    if((im==0) & (ic==0)) subfsx0_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(is, sx_is, data);
    if((im==0) & (ic==1)) subfsx0_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(is, sx_is, data);
    if((im==0) & (ic==2)) subfsx0_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(is, sx_is, data);
    if((im==0) & (ic==3)) subfsx0_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(is, sx_is, data);
    if((im==0) & (ic==4)) subfsx0_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(is, sx_is, data);
    if((im==0) & (ic==5)) subfsx0_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(is, sx_is, data);
  }
}

 void csv(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) csv_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, x, ip, sx, data);
  if((im==0) & (ic==1)) csv_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, x, ip, sx, data);
  if((im==0) & (ic==2)) csv_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, x, ip, sx, data);
  if((im==0) & (ic==3)) csv_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, x, ip, sx, data);
  if((im==0) & (ic==4)) csv_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, x, ip, sx, data);
  if((im==0) & (ic==5)) csv_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, x, ip, sx, data);
}

 int AR_CVodeSensInit1(void *cvode_mem, int nps, int sensi_meth, int sensirhs, N_Vector *sx, int im, int ic, int sensitivitySubset){
  if (sensirhs == 1) {
    if (sensitivitySubset == 0) {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7, sx);
      if((im==0) & (ic==1)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB, sx);
      if((im==0) & (ic==2)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6, sx);
      if((im==0) & (ic==3)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B, sx);
      if((im==0) & (ic==4)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_AktPathwayFujita_682D2488568EE2F473172B1526A835D5, sx);
      if((im==0) & (ic==5)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB, sx);
    } else {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7, sx);
      if((im==0) & (ic==1)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB, sx);
      if((im==0) & (ic==2)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6, sx);
      if((im==0) & (ic==3)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B, sx);
      if((im==0) & (ic==4)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_AktPathwayFujita_682D2488568EE2F473172B1526A835D5, sx);
      if((im==0) & (ic==5)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB, sx);
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
  if((im==0) & (ic==0)) fu_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(data, t);
  if((im==0) & (ic==1)) fu_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(data, t);
  if((im==0) & (ic==2)) fu_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(data, t);
  if((im==0) & (ic==3)) fu_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(data, t);
  if((im==0) & (ic==4)) fu_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(data, t);
  if((im==0) & (ic==5)) fu_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(data, t);
}

 void fsu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fsu_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(data, t);
  if((im==0) & (ic==1)) fsu_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(data, t);
  if((im==0) & (ic==2)) fsu_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(data, t);
  if((im==0) & (ic==3)) fsu_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(data, t);
  if((im==0) & (ic==4)) fsu_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(data, t);
  if((im==0) & (ic==5)) fsu_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(data, t);
}

 void fv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fv_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, x, data);
  if((im==0) & (ic==1)) fv_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, x, data);
  if((im==0) & (ic==2)) fv_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, x, data);
  if((im==0) & (ic==3)) fv_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, x, data);
  if((im==0) & (ic==4)) fv_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, x, data);
  if((im==0) & (ic==5)) fv_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, x, data);
}

 void fsv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
	dvdp_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, x, data);
	dvdu_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, x, data);
	dvdx_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, x, data);
}
  if((im==0) & (ic==1)) {
	dvdp_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, x, data);
	dvdu_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, x, data);
	dvdx_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, x, data);
}
  if((im==0) & (ic==2)) {
	dvdp_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, x, data);
	dvdu_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, x, data);
	dvdx_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, x, data);
}
  if((im==0) & (ic==3)) {
	dvdp_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, x, data);
	dvdu_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, x, data);
	dvdx_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, x, data);
}
  if((im==0) & (ic==4)) {
	dvdp_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, x, data);
	dvdu_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, x, data);
	dvdx_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, x, data);
}
  if((im==0) & (ic==5)) {
	dvdp_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, x, data);
	dvdu_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, x, data);
	dvdx_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, x, data);
}
}

 void dfxdp0(void *user_data, double t, N_Vector x, double *dfxdp0, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp0_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, x, dfxdp0, data);
  if((im==0) & (ic==1)) dfxdp0_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, x, dfxdp0, data);
  if((im==0) & (ic==2)) dfxdp0_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, x, dfxdp0, data);
  if((im==0) & (ic==3)) dfxdp0_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, x, dfxdp0, data);
  if((im==0) & (ic==4)) dfxdp0_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, x, dfxdp0, data);
  if((im==0) & (ic==5)) dfxdp0_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, x, dfxdp0, data);
}

 void dfxdp(void *user_data, double t, N_Vector x, double *dfxdp, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, x, dfxdp, data);
  if((im==0) & (ic==1)) dfxdp_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, x, dfxdp, data);
  if((im==0) & (ic==2)) dfxdp_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, x, dfxdp, data);
  if((im==0) & (ic==3)) dfxdp_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, x, dfxdp, data);
  if((im==0) & (ic==4)) dfxdp_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, x, dfxdp, data);
  if((im==0) & (ic==5)) dfxdp_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, x, dfxdp, data);
}

void fz(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) fz_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==1)) fz_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==2)) fz_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==3)) fz_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==4)) fz_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==5)) fz_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
}

void dfzdx(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdx, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) dfzdx_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==1)) dfzdx_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==2)) dfzdx_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==3)) dfzdx_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==4)) dfzdx_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==5)) dfzdx_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
}

void fsz(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx, int im, int ic){
  if((im==0) & (ic==0)) fsz_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==1)) fsz_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==2)) fsz_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==3)) fsz_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==4)) fsz_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==5)) fsz_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, nt, it, np, sz, p, u, x, z, su, sx);
}

 void fy(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fy_experimentaldata1_258544EB6B3AA8CBD3145AB82F974965(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==1)) fy_experimentaldata2_D90EAD8994BDEFB35ABE1905F22F8E93(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==2)) fy_experimentaldata3_C230EC340F7FFF25812B99F99870D51A(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==3)) fy_experimentaldata4_5F66C2E0235D12505358B893E4F722CF(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==4)) fy_experimentaldata5_5B4FCB1D8192A23C5B469A3F75AC4E18(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==5)) fy_experimentaldata6_29AF108A47A7B9DE3D5BF91CBED9712A(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
}

 void fy_scale(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx, int im, int id){
  if((im==0) & (id==0)) fy_scale_experimentaldata1_258544EB6B3AA8CBD3145AB82F974965(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==1)) fy_scale_experimentaldata2_D90EAD8994BDEFB35ABE1905F22F8E93(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==2)) fy_scale_experimentaldata3_C230EC340F7FFF25812B99F99870D51A(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==3)) fy_scale_experimentaldata4_5F66C2E0235D12505358B893E4F722CF(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==4)) fy_scale_experimentaldata5_5B4FCB1D8192A23C5B469A3F75AC4E18(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==5)) fy_scale_experimentaldata6_29AF108A47A7B9DE3D5BF91CBED9712A(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
}

 void fystd(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fystd_experimentaldata1_258544EB6B3AA8CBD3145AB82F974965(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==1)) fystd_experimentaldata2_D90EAD8994BDEFB35ABE1905F22F8E93(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==2)) fystd_experimentaldata3_C230EC340F7FFF25812B99F99870D51A(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==3)) fystd_experimentaldata4_5F66C2E0235D12505358B893E4F722CF(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==4)) fystd_experimentaldata5_5B4FCB1D8192A23C5B469A3F75AC4E18(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==5)) fystd_experimentaldata6_29AF108A47A7B9DE3D5BF91CBED9712A(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
}

 void fsy(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsy_experimentaldata1_258544EB6B3AA8CBD3145AB82F974965(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==1)) fsy_experimentaldata2_D90EAD8994BDEFB35ABE1905F22F8E93(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==2)) fsy_experimentaldata3_C230EC340F7FFF25812B99F99870D51A(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==3)) fsy_experimentaldata4_5F66C2E0235D12505358B893E4F722CF(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==4)) fsy_experimentaldata5_5B4FCB1D8192A23C5B469A3F75AC4E18(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==5)) fsy_experimentaldata6_29AF108A47A7B9DE3D5BF91CBED9712A(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
}

 void fsystd(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsystd_experimentaldata1_258544EB6B3AA8CBD3145AB82F974965(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==1)) fsystd_experimentaldata2_D90EAD8994BDEFB35ABE1905F22F8E93(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==2)) fsystd_experimentaldata3_C230EC340F7FFF25812B99F99870D51A(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==3)) fsystd_experimentaldata4_5F66C2E0235D12505358B893E4F722CF(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==4)) fsystd_experimentaldata5_5B4FCB1D8192A23C5B469A3F75AC4E18(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==5)) fsystd_experimentaldata6_29AF108A47A7B9DE3D5BF91CBED9712A(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
}

/* for arSSACalc.c */

 void fvSSA(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
    fu_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(data, t);
    fv_AktPathwayFujita_411D2F712EB6E486B5E602E77F7FA1F7(t, x, data);
  }
  if((im==0) & (ic==1)) {
    fu_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(data, t);
    fv_AktPathwayFujita_1029F17CBF497C455F622528CD8B13FB(t, x, data);
  }
  if((im==0) & (ic==2)) {
    fu_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(data, t);
    fv_AktPathwayFujita_D49E5244DFFC02E15C3D15EE45C242B6(t, x, data);
  }
  if((im==0) & (ic==3)) {
    fu_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(data, t);
    fv_AktPathwayFujita_3A5AAFEA6D669B362C447AD33A15191B(t, x, data);
  }
  if((im==0) & (ic==4)) {
    fu_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(data, t);
    fv_AktPathwayFujita_682D2488568EE2F473172B1526A835D5(t, x, data);
  }
  if((im==0) & (ic==5)) {
    fu_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(data, t);
    fv_AktPathwayFujita_5E034222D21D0AC5A4F50E75188F52AB(t, x, data);
  }
}

