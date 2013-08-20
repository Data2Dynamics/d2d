#include <math.h>

/* general input functions */

double step1(double t, double u1, double t1, double u2);
double dstep1(double t, double u1, double t1, double u2, int p_index);

double step2(double t, double u1, double t1, double u2, double t2, double u3);
double dstep2(double t, double u1, double t1, double u2, double t2, double u3, int p_index);

/* custom rate laws */

double mmenten(double x, double vmax, double km);
double mmenten_alt(double x, double klin, double ksat);

double hill_kd(double x, double h, double kd);
double hill_ka(double x, double h, double ka);
