#ifndef _MY_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2
#define _MY_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2

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



 void fy_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z);
 void fystd_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z);
 void fsy_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz);
 void fsystd_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz);

 void fy_scale_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx);
#endif /* _MY_ELISA_Nigericin_formatted_WT_6BF70393A8A997A3F72D6F1812FCB6A2 */



