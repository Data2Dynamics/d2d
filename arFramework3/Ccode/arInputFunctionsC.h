#ifndef _ARINPUTFUNCTIONS_C_
#define _ARINPUTFUNCTIONS_C_

#include <math.h>

/* general input functions */
double heaviside(double t);
double dirac(double t);

double LUT_bilinear( double x, double y, int NX, int NY, const double data[] );
double DLUT_bilinear( double x, double y, int NX, int NY, const double data[], int deriv );
double getData2D( const int NX, const int NY, const double data[], int x, int y );

double step1(double t, double u1, double t1, double u2);
double dstep1(double t, double u1, double t1, double u2, int p_index);

double step2(double t, double u1, double t1, double u2, double t2, double u3);
double dstep2(double t, double u1, double t1, double u2, double t2, double u3, int p_index);

/* Spline which optionally allows for caching */
double inputspline( double t, const int n, const double ts[], const double us[]);
double inputfastspline( double t, int ID, double **splineCache, int *idCache, const int n, const double ts[], const double us[]);

int cspline(int n, int end1, int end2, double slope1, double slope2, const double x[], const double y[], double b[], double c[], double d[], int cacheID, double **splineCache, int *IDcache);
int cmonotoneSpline( int n, const double x[], const double y[], double b[], double c[], double d[], int cacheID, double **splineCache, int *IDcache );
int clongmonotoneSpline( int n, const double x[], const double y[], double b[], double c[], double d[], int cacheID, double **splineCache, int *IDcache );

/* splines */
double spline3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
double spline_pos3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
double monospline3(double t, double t1, double p1, double t2, double p2, double t3, double p3);

double spline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
double spline_pos4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
double monospline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4);

double spline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
double spline_pos5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
double monospline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5);

double spline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
double spline_pos10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
double monospline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10);

double spline20(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10,double t11, double p11, double t12, double p12, double t13, double p13, double t14, double p14, double t15, double p15, double t16, double p16, double t17, double p17, double t18, double p18, double t19, double p19, double t20, double p20, int ss, double dudt);

/* Implementations which store the coefficients */
double fastspline3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
double fastspline_pos3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
double monofastspline3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3);

double fastspline4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
double fastspline_pos4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
double monofastspline4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4);

double fastspline5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
double fastspline_pos5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
double monofastspline5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5);

double fastspline10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
double fastspline_pos10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
double monofastspline10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10);

/* spline derivatives */
double Dspline3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt, int id);
double Dspline_pos3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt, int id);
double Dmonospline3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int id);

double Dspline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt, int id);
double Dspline_pos4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt, int id);
double Dmonospline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int id);

double Dspline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt, int id);
double Dspline_pos5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt, int id);
double Dmonospline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int id);

double Dspline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt, int id);
double Dspline_pos10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt, int id);
double Dmonospline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int id);


double Dspline20(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, double t11, double p11, double t12, double p12, double t13, double p13, double t14, double p14, double t15, double p15, double t16, double p16, double t17, double p17, double t18, double p18, double t19, double p19, double t20, double p20, int ss, double dudt, int id);

/* Implementations which store the coefficients */
double Dfastspline3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt, int id);
double Dfastspline_pos3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt, int id);
double Dmonofastspline3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, int id);

double Dfastspline4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt, int id);
double Dfastspline_pos4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt, int id);
double Dmonofastspline4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int id);

double Dfastspline5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt, int id);
double Dfastspline_pos5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt, int id);
double Dmonofastspline5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int id);

double Dfastspline10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt, int id);
double Dfastspline_pos10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt, int id);
double Dmonofastspline10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int id);

/* A special spline implementation where coefficients are given directly */
double splineFixCoeffs( double t, int n, const double time[], const double data[] );
double interpolateLinear( double t, int n, const double time[], const double data[] );

/* custom rate laws */
double mmenten(double x, double vmax, double km);
double mmenten_alt(double x, double klin, double ksat);

double hill_kd(double x, double h, double kd);
double hill_ka(double x, double h, double ka);

#endif
