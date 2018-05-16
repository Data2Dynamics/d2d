#include "erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8.h"
#include "erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1.h"
#include "erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E.h"
#include "erbb_signaling_C63B01C16263D42E082581726F3DC867.h"
#include "experimentaldata1_DF16BB96AF0E4B734CE722B20F7F1F94.h"
#include "experimentaldata2_7759CA296E423A7852EF498129D9C9C4.h"
#include "experimentaldata3_6CF63E57231D682CB322A83168E431BF.h"
#include "experimentaldata4_CF2545A583D43C1FCB12219952CBA49D.h"

 int AR_CVodeInit(void *cvode_mem, N_Vector x, double t, int im, int ic){
  if((im==0) & (ic==0)) return CVodeInit(cvode_mem, fx_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8, RCONST(t), x);
  if((im==0) & (ic==1)) return CVodeInit(cvode_mem, fx_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1, RCONST(t), x);
  if((im==0) & (ic==2)) return CVodeInit(cvode_mem, fx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E, RCONST(t), x);
  if((im==0) & (ic==3)) return CVodeInit(cvode_mem, fx_erbb_signaling_C63B01C16263D42E082581726F3DC867, RCONST(t), x);
  return(-1);
}

 void fx(realtype t, N_Vector x, double *xdot, void *user_data, int im, int ic){
  if((im==0) & (ic==0)) fxdouble_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, x, xdot, user_data);
  if((im==0) & (ic==1)) fxdouble_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, x, xdot, user_data);
  if((im==0) & (ic==2)) fxdouble_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, x, xdot, user_data);
  if((im==0) & (ic==3)) fxdouble_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, x, xdot, user_data);
}

 void fx0(N_Vector x0, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fx0_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(x0, data);
  if((im==0) & (ic==1)) fx0_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(x0, data);
  if((im==0) & (ic==2)) fx0_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(x0, data);
  if((im==0) & (ic==3)) fx0_erbb_signaling_C63B01C16263D42E082581726F3DC867(x0, data);
}

 int AR_CVDlsSetDenseJacFn(void *cvode_mem, int im, int ic, int setSparse){
  if((im==0) & (ic==0) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8);
 
 }else if((im==0) & (ic==0) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8);

}
  if((im==0) & (ic==1) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1);
 
 }else if((im==0) & (ic==1) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1);

}
  if((im==0) & (ic==2) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E);
 
 }else if((im==0) & (ic==2) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E);

}
  if((im==0) & (ic==3) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_erbb_signaling_C63B01C16263D42E082581726F3DC867);
 
 }else if((im==0) & (ic==3) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_erbb_signaling_C63B01C16263D42E082581726F3DC867);

}
  return(-1);
}

 void getdfxdx(int im, int ic, realtype t, N_Vector x, realtype *J, void *user_data){
  if((im==0) & (ic==0)) { dfxdx_out_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, x, J, user_data); return; }
  if((im==0) & (ic==1)) { dfxdx_out_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, x, J, user_data); return; }
  if((im==0) & (ic==2)) { dfxdx_out_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, x, J, user_data); return; }
  if((im==0) & (ic==3)) { dfxdx_out_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, x, J, user_data); return; }

}

 void fsx0(int is, N_Vector sx_is, void *user_data, int im, int ic, int sensitivitySubset) {
  UserData data = (UserData) user_data;
  if ( sensitivitySubset == 0 ) {
    if((im==0) & (ic==0)) fsx0_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(is, sx_is, data);
    if((im==0) & (ic==1)) fsx0_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(is, sx_is, data);
    if((im==0) & (ic==2)) fsx0_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(is, sx_is, data);
    if((im==0) & (ic==3)) fsx0_erbb_signaling_C63B01C16263D42E082581726F3DC867(is, sx_is, data);
  } else {
    if((im==0) & (ic==0)) subfsx0_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(is, sx_is, data);
    if((im==0) & (ic==1)) subfsx0_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(is, sx_is, data);
    if((im==0) & (ic==2)) subfsx0_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(is, sx_is, data);
    if((im==0) & (ic==3)) subfsx0_erbb_signaling_C63B01C16263D42E082581726F3DC867(is, sx_is, data);
  }
}

 void csv(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) csv_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, x, ip, sx, data);
  if((im==0) & (ic==1)) csv_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, x, ip, sx, data);
  if((im==0) & (ic==2)) csv_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, x, ip, sx, data);
  if((im==0) & (ic==3)) csv_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, x, ip, sx, data);
}

 int AR_CVodeSensInit1(void *cvode_mem, int nps, int sensi_meth, int sensirhs, N_Vector *sx, int im, int ic, int sensitivitySubset){
  if (sensirhs == 1) {
    if (sensitivitySubset == 0) {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8, sx);
      if((im==0) & (ic==1)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1, sx);
      if((im==0) & (ic==2)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E, sx);
      if((im==0) & (ic==3)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_erbb_signaling_C63B01C16263D42E082581726F3DC867, sx);
    } else {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8, sx);
      if((im==0) & (ic==1)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1, sx);
      if((im==0) & (ic==2)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E, sx);
      if((im==0) & (ic==3)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_erbb_signaling_C63B01C16263D42E082581726F3DC867, sx);
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
  if((im==0) & (ic==0)) fu_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(data, t);
  if((im==0) & (ic==1)) fu_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(data, t);
  if((im==0) & (ic==2)) fu_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(data, t);
  if((im==0) & (ic==3)) fu_erbb_signaling_C63B01C16263D42E082581726F3DC867(data, t);
}

 void fsu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fsu_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(data, t);
  if((im==0) & (ic==1)) fsu_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(data, t);
  if((im==0) & (ic==2)) fsu_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(data, t);
  if((im==0) & (ic==3)) fsu_erbb_signaling_C63B01C16263D42E082581726F3DC867(data, t);
}

 void fv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fv_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, x, data);
  if((im==0) & (ic==1)) fv_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, x, data);
  if((im==0) & (ic==2)) fv_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, x, data);
  if((im==0) & (ic==3)) fv_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, x, data);
}

 void fsv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
	dvdp_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, x, data);
	dvdu_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, x, data);
	dvdx_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, x, data);
}
  if((im==0) & (ic==1)) {
	dvdp_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, x, data);
	dvdu_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, x, data);
	dvdx_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, x, data);
}
  if((im==0) & (ic==2)) {
	dvdp_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, x, data);
	dvdu_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, x, data);
	dvdx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, x, data);
}
  if((im==0) & (ic==3)) {
	dvdp_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, x, data);
	dvdu_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, x, data);
	dvdx_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, x, data);
}
}

 void dfxdp0(void *user_data, double t, N_Vector x, double *dfxdp0, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp0_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, x, dfxdp0, data);
  if((im==0) & (ic==1)) dfxdp0_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, x, dfxdp0, data);
  if((im==0) & (ic==2)) dfxdp0_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, x, dfxdp0, data);
  if((im==0) & (ic==3)) dfxdp0_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, x, dfxdp0, data);
}

 void dfxdp(void *user_data, double t, N_Vector x, double *dfxdp, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, x, dfxdp, data);
  if((im==0) & (ic==1)) dfxdp_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, x, dfxdp, data);
  if((im==0) & (ic==2)) dfxdp_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, x, dfxdp, data);
  if((im==0) & (ic==3)) dfxdp_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, x, dfxdp, data);
}

void fz(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) fz_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==1)) fz_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==2)) fz_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
  if((im==0) & (ic==3)) fz_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
}

void dfzdx(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdx, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) dfzdx_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==1)) dfzdx_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==2)) dfzdx_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
  if((im==0) & (ic==3)) dfzdx_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
}

void fsz(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx, int im, int ic){
  if((im==0) & (ic==0)) fsz_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==1)) fsz_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==2)) fsz_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, nt, it, np, sz, p, u, x, z, su, sx);
  if((im==0) & (ic==3)) fsz_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, nt, it, np, sz, p, u, x, z, su, sx);
}

 void fy(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fy_experimentaldata1_DF16BB96AF0E4B734CE722B20F7F1F94(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==1)) fy_experimentaldata2_7759CA296E423A7852EF498129D9C9C4(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==2)) fy_experimentaldata3_6CF63E57231D682CB322A83168E431BF(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==3)) fy_experimentaldata4_CF2545A583D43C1FCB12219952CBA49D(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
}

 void fy_scale(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx, int im, int id){
  if((im==0) & (id==0)) fy_scale_experimentaldata1_DF16BB96AF0E4B734CE722B20F7F1F94(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==1)) fy_scale_experimentaldata2_7759CA296E423A7852EF498129D9C9C4(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==2)) fy_scale_experimentaldata3_6CF63E57231D682CB322A83168E431BF(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==3)) fy_scale_experimentaldata4_CF2545A583D43C1FCB12219952CBA49D(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
}

 void fystd(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fystd_experimentaldata1_DF16BB96AF0E4B734CE722B20F7F1F94(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==1)) fystd_experimentaldata2_7759CA296E423A7852EF498129D9C9C4(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==2)) fystd_experimentaldata3_6CF63E57231D682CB322A83168E431BF(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==3)) fystd_experimentaldata4_CF2545A583D43C1FCB12219952CBA49D(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
}

 void fsy(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsy_experimentaldata1_DF16BB96AF0E4B734CE722B20F7F1F94(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==1)) fsy_experimentaldata2_7759CA296E423A7852EF498129D9C9C4(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==2)) fsy_experimentaldata3_6CF63E57231D682CB322A83168E431BF(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==3)) fsy_experimentaldata4_CF2545A583D43C1FCB12219952CBA49D(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
}

 void fsystd(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsystd_experimentaldata1_DF16BB96AF0E4B734CE722B20F7F1F94(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==1)) fsystd_experimentaldata2_7759CA296E423A7852EF498129D9C9C4(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==2)) fsystd_experimentaldata3_6CF63E57231D682CB322A83168E431BF(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==3)) fsystd_experimentaldata4_CF2545A583D43C1FCB12219952CBA49D(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
}

/* for arSSACalc.c */

 void fvSSA(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
    fu_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(data, t);
    fv_erbb_signaling_D1D470B8BE2BAA4CBBCB033DC9BD3EE8(t, x, data);
  }
  if((im==0) & (ic==1)) {
    fu_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(data, t);
    fv_erbb_signaling_F4DBEF6898F3E929B8764DD34A2756C1(t, x, data);
  }
  if((im==0) & (ic==2)) {
    fu_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(data, t);
    fv_erbb_signaling_BD90C6D916A5B71EA8DA7808C055872E(t, x, data);
  }
  if((im==0) & (ic==3)) {
    fu_erbb_signaling_C63B01C16263D42E082581726F3DC867(data, t);
    fv_erbb_signaling_C63B01C16263D42E082581726F3DC867(t, x, data);
  }
}

