#include <math.h>

/* general input functions */

double step1(double t, double u1, double t1, double u2) {
    if(t <= t1) {
        return(u1);
    } else {
        return(u2);
    }
}
double dstep1(double t, double u1, double t1, double u2, int p_index) {
    if(t <= t1) {
        if(p_index == 1)
            return(1);
        else
            return(0);
    } else {
        if(p_index == 3)
            return(1);
        else
            return(0);
    }
}

double step2(double t, double u1, double t1, double u2, double t2, double u3) {
    if(t <= t1) {
        return(u1);
    } else if(t > t1 & t <= t2) {
        return(u2);
    } else {
        return(u3);
    }
}
double dstep2(double t, double u1, double t1, double u2, double t2, double u3, int p_index) {
    if(t <= t1) {
        if(p_index == 1)
            return(1);
        else
            return(0);
    } else if(t > t1 & t <= t2) {
        if(p_index == 2)
            return(1);
        else
            return(0);
    } else {
        if(p_index == 3)
            return(1);
        else
            return(0);
    }
}

/* custom rate laws */

double mmenten(double x, double vmax, double km){
    return(vmax * x / (km + x));
}

double mmenten_alt(double x, double klin, double ksat){
    return(0);
}

double hill_kd(double x, double h, double kd){
    return(pow(x,h) / (kd + pow(x,h)));
}

double hill_ka(double x, double h, double ka){
    /* return(pow(x,h) / (pow(ka,h) + pow(x,h))); */
    return(1 / (pow((ka/x),h) + 1));
}


