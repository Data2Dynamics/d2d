#include "ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6.h"
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





 void fy_ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z){
  y[ny*nt*iruns+it+nt*0] = p[6]-(pow(3.0,p[7])*(exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))-1.0)*p[3])/(pow(3.0,p[7])+pow(p[0],p[7]));

  return;
}


 void fystd_ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z){
  ystd[it+nt*0] = p[10];

  return;
}


 void fsy_ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz){
  sy[it+nt*0] = pow(3.0,p[7])*(exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))-1.0)*1.0/pow(pow(3.0,p[7])+pow(p[0],p[7]),2.0)*pow(p[0],p[7]-1.0)*p[3]*p[7];
  sy[it+nt*1] = (pow(3.0,p[7])*pow(3.0,p[8])*exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*1.0/pow(pow(3.0,p[8])+pow(p[1],p[8]),2.0)*pow(p[1],p[8]-1.0)*p[3]*p[4]*p[8])/(pow(3.0,p[7])+pow(p[0],p[7]));
  sy[it+nt*2] = -(pow(3.0,p[7])*pow(3.0,p[8])*exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))*((pow(3.0,p[9])*pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))*1.0/pow(pow(3.0,p[9])+pow(p[2],p[9]),2.0)*pow(p[2],p[9]-1.0)*p[5]*p[9])/(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)-(pow(3.0,p[9])*pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))*1.0/pow(pow(3.0,p[9])+pow(p[2],p[9]),2.0)*pow(p[2],p[9]-1.0)*p[5]*p[9])/(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))))*p[3]*p[4])/((pow(3.0,p[7])+pow(p[0],p[7]))*(pow(3.0,p[8])+pow(p[1],p[8])));
  sy[it+nt*3] = -(pow(3.0,p[7])*(exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))-1.0))/(pow(3.0,p[7])+pow(p[0],p[7]));
  sy[it+nt*4] = -(pow(3.0,p[7])*pow(3.0,p[8])*exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[3])/((pow(3.0,p[7])+pow(p[0],p[7]))*(pow(3.0,p[8])+pow(p[1],p[8])));
  sy[it+nt*5] = (pow(3.0,p[7])*pow(3.0,p[8])*exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))*((pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))/(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)-(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))/(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))))*p[3]*p[4])/((pow(3.0,p[7])+pow(p[0],p[7]))*(pow(3.0,p[8])+pow(p[1],p[8])));
  sy[it+nt*6] = 1.0;
  sy[it+nt*7] = -(pow(3.0,p[7])*log(3.0)*(exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))-1.0)*p[3])/(pow(3.0,p[7])+pow(p[0],p[7]))+pow(3.0,p[7])*(exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))-1.0)*1.0/pow(pow(3.0,p[7])+pow(p[0],p[7]),2.0)*p[3]*(log(p[0])*pow(p[0],p[7])+pow(3.0,p[7])*log(3.0));
  sy[it+nt*8] = -(pow(3.0,p[7])*exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))*p[3]*((pow(3.0,p[8])*log(3.0)*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8]))-pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*1.0/pow(pow(3.0,p[8])+pow(p[1],p[8]),2.0)*p[4]*(log(p[1])*pow(p[1],p[8])+pow(3.0,p[8])*log(3.0))))/(pow(3.0,p[7])+pow(p[0],p[7]));
  sy[it+nt*9] = -(pow(3.0,p[7])*pow(3.0,p[8])*exp((pow(3.0,p[8])*(log(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)/log(1.0E+1)-log(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0)))/log(1.0E+1))*p[4])/(pow(3.0,p[8])+pow(p[1],p[8])))*((pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))*(pow(3.0,p[9])*1.0/pow(pow(3.0,p[9])+pow(p[2],p[9]),2.0)*(log(p[2])*pow(p[2],p[9])+pow(3.0,p[9])*log(3.0))-(pow(3.0,p[9])*log(3.0))/(pow(3.0,p[9])+pow(p[2],p[9])))*p[5])/(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))+1.0)-(pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))*(pow(3.0,p[9])*1.0/pow(pow(3.0,p[9])+pow(p[2],p[9]),2.0)*(log(p[2])*pow(p[2],p[9])+pow(3.0,p[9])*log(3.0))-(pow(3.0,p[9])*log(3.0))/(pow(3.0,p[9])+pow(p[2],p[9])))*p[5])/(pow(1.0E+1,t*(5.0/3.0))+pow(1.0E+1,-p[5]*(pow(3.0,p[9])/(pow(3.0,p[9])+pow(p[2],p[9]))-1.0))))*p[3]*p[4])/((pow(3.0,p[7])+pow(p[0],p[7]))*(pow(3.0,p[8])+pow(p[1],p[8])));
  sy[it+nt*10] = 0.0;

  return;
}


 void fsystd_ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz){
  systd[it+nt*0] = 0.0;
  systd[it+nt*1] = 0.0;
  systd[it+nt*2] = 0.0;
  systd[it+nt*3] = 0.0;
  systd[it+nt*4] = 0.0;
  systd[it+nt*5] = 0.0;
  systd[it+nt*6] = 0.0;
  systd[it+nt*7] = 0.0;
  systd[it+nt*8] = 0.0;
  systd[it+nt*9] = 0.0;
  systd[it+nt*10] = 1.0;

  return;
}


 void fy_scale_ELISA_Nigericin_formatted_WT_4827847680731DB1788022E3CE15D3E6(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx){


  return;
}


