#include <math.h>
#include "spline.c"
#include "monotone.c"
#include "arInputFunctionsC.h"

/* general input functions */
double heaviside(double t) {
    if (t < 0)
        return(0.0);
    else
        return(1.0);
}

/* delta functions are silently dropped */
double dirac(double t) {
    return 0.0;
}

/* Helper function for clamped 2D array access */
double getData2D( const int NX, const int NY, const double data[], int x, int y )
{
    /* Prevent exceeding the data range */
    x = x < NX ? x : NX - 1;
    y = y < NY ? y : NY - 1;
    x = x > 0 ? x : 0;
    y = y > 0 ? y : 0;

    return data[ x + NX * y ];
}

/* Bilinear LUT */
double LUT_bilinear( double x, double y, int NX, int NY, const double data[] )
{
    double xmi = floor( x * (NX-1) );
    double xma = ceil( x * (NX-1) );
    double ymi = floor( y * (NY-1) );
    double yma = ceil( y * (NY-1) );
    double dx = x * (NX-1) - xmi;
    double dy = y * (NY-1) - ymi;
    
    double y1 = ( 1 - dx ) * getData2D( NX, NY, data, xmi, ymi ) + dx * getData2D( NX, NY, data, xma, ymi );
    double y2 = ( 1 - dx ) * getData2D( NX, NY, data, xmi, yma ) + dx * getData2D( NX, NY, data, xma, yma );
    
    double out = ( 1 - dy ) * y1 + dy * y2;
    
    return out;
}

/* Derivatives of the bilinear LUT */
double DLUT_bilinear( double x, double y, int NX, int NY, const double data[], int deriv )
{
    double xmi = floor( x * (NX-1) );
    double xma = ceil( x * (NX-1) );
    double ymi = floor( y * (NY-1) );
    double yma = ceil( y * (NY-1) );
    
    double dx = x * (NX-1) - xmi;
    double dy = y * (NY-1) - ymi;
    
    double fxmiymi = getData2D( NX, NY, data, xmi, ymi );
    double fxmiyma = getData2D( NX, NY, data, xmi, yma );
    double fxmaymi = getData2D( NX, NY, data, xma, ymi );
    double fxmayma = getData2D( NX, NY, data, xma, yma );
        
    /* Derivative w.r.t. x (1) or y (2)? */
    if ( deriv == 1 )
        return (fxmaymi*(NX - 1) - fxmiymi*(NX - 1))*(ymi - y*(NY - 1) + 1) - (fxmayma*(NX - 1) - fxmiyma*(NX - 1))*(ymi - y*(NY - 1));
    else
        return (fxmaymi*(xmi - x*(NX - 1)) - fxmiymi*(xmi - x*(NX - 1) + 1))*(NY - 1) - (fxmayma*(xmi - x*(NX - 1)) - fxmiyma*(xmi - x*(NX - 1) + 1))*(NY - 1);
}

/* Spline with fixed time points and coefficients */
double inputfastspline( double t, int ID, double **splineCache, int *idCache, const int n, const double ts[], const double us[])
{
    double uout;
    double b[MAX_LONG_SPLINE];
    double c[MAX_LONG_SPLINE];
    double d[MAX_LONG_SPLINE];
    
    clongmonotoneSpline( n, ts, us, b, c, d, ID, splineCache, idCache );
    uout = seval_fixed( n, t, ts, us, b, c, d, &(idCache[ID]));
    
    return uout;
}

double interpolateLinear( double t, int n, const double time[], const double data[] )
{
    double *deriv2;
    double val;
    
    int max             = n-1;
    int min             = 0;
    int cur             = max/2;

    double a, b, h, tmp;

    /* Out of bounds? */
    if ( t >= time[ max ] )
        return data[ n-1 ];
    else if ( t <= time[ min ] )
        return data[ 0 ];

    /* Where are we? */
    while ( min != max - 1 )
    {
        if ( time[ cur ] <= t )
            min = cur;
        else
            max = cur;

        cur = min + ( max - min ) / 2;
    }

    val = data[min] + ( t - time[min] ) * ( data[max] - data[min] ) / ( time[max] - time[min] );
    return val;
}

/* Fixed coefficient cubic spline */
double splineFixCoeffs( double t, int n, const double time[], const double data[] )
{
    double *deriv2;
    double val;
    
    int max             = n-1;
    int min             = 0;
    int cur             = max/2;

    double a, b, h, tmp;

    /* Out of bounds? */
    if ( t >= time[ max ] )
        return data[ n-1 ];
    else if ( t <= time[ min ] )
        return data[ 0 ];

    /* Where are we? */
    while ( min != max - 1 )
    {
        if ( time[ cur ] <= t )
            min = cur;
        else
            max = cur;

        cur = min + ( max - min ) / 2;
    }

    deriv2 = &data[n];
    h = time[max] - time[min];
    a = ( time[ max ] - t ) / h;
    b = ( t - time[ min ] ) / h;
    val = a * data[ min ] + b * data[ max ] + ( ( a * a * a - a ) * deriv2[ min ] + ( b * b * b - b ) * deriv2[ max ] ) * ( h * h ) / 6.0;
    
    return val;
}

/* Spline with fixed time points and coefficients */
double inputspline( double t, const int n, const double ts[], const double us[])
{
    double uout;
    
    double b[MAX_LONG_SPLINE];
    double c[MAX_LONG_SPLINE];
    double d[MAX_LONG_SPLINE];
 
    longMonotoneSpline( n, ts, us, b, c, d );
    uout = seval( n, t, ts, us, b, c, d );
    
    return uout;
}

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
            return(1.0);
        else
            return(0.0);
    } else {
        if(p_index == 3)
            return(1.0);
        else
            return(0.0);
    }
}

double step2(double t, double u1, double t1, double u2, double t2, double u3) {
    if(t <= t1) {
        return(u1);
    } else if( (t > t1) & (t <= t2) ) {
        return(u2);
    } else {
        return(u3);
    }
}
double dstep2(double t, double u1, double t1, double u2, double t2, double u3, int p_index) {
    if(t <= t1) {
        if(p_index == 1)
            return(1.0);
        else
            return(0);
    } else if( (t > t1) & (t <= t2) ) {
        if(p_index == 2)
            return(1.0);
        else
            return(0.0);
    } else {
        if(p_index == 3)
            return(1.0);
        else
            return(0.0);
    }
}

/* splines */
double spline3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt) {   
    double uout;
    
    double ts[3];
    double us[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    
    spline(3, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(3, t, ts, us, b, c, d);
    
    return(uout);
}

double monospline3(double t, double t1, double p1, double t2, double p2, double t3, double p3) {   
    double uout;
    
    double ts[3];
    double us[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;

    monotoneSpline( 3, ts, us, b, c, d );
    uout = seval(3, t, ts, us, b, c, d);
    
    return(uout);
}

double spline_pos3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt) {   
    int is;
    double uout;
    
    double ts[3];
    double us[3];
    double uslog[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    
    for (is = 0; is<3; is++){
        uslog[is] = log(us[is]);
    }
    
    spline(3, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uout = seval(3, t, ts, uslog, b, c, d);
    
    return(exp(uout));
}

double spline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt) {   
    double uout;
    
    double ts[4];
    double us[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    
    spline(4, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(4, t, ts, us, b, c, d);
    
    return(uout);
}

double monospline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4) {   
    double uout;
    
    double ts[4];
    double us[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;

    monotoneSpline( 4, ts, us, b, c, d );
    uout = seval(4, t, ts, us, b, c, d);
    
    return(uout);
}

double spline_pos4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt) {   
    int is;
    double uout;
    
    double ts[4];
    double us[4];
    double uslog[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    
    for (is = 0; is<4; is++){
        uslog[is] = log(us[is]);
    }
    
    spline(4, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uout = seval(4, t, ts, uslog, b, c, d);
    
    return(exp(uout));
}

double spline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt) {   
    double uout;
    
    double ts[5];
    double us[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    
    spline(5, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(5, t, ts, us, b, c, d);
    
    return(uout);
}

double monospline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5) {   
    double uout;
    
    double ts[5];
    double us[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;

    monotoneSpline( 5, ts, us, b, c, d );
    uout = seval(5, t, ts, us, b, c, d);
    
    return(uout);
}

double spline_pos5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt) {   
    int is;
    double uout;
    
    double ts[5];
    double us[5];
    double uslog[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    
    for (is = 0; is<5; is++){
        uslog[is] = log(us[is]);
    }
    
    spline(5, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uout = seval(5, t, ts, uslog, b, c, d);
    
    return(exp(uout));
}

double spline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt) {   
    double uout;
    
    double ts[10];
    double us[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10;
    
    spline(10, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(10, t, ts, us, b, c, d);
    
    return(uout);
}

double spline15(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, double t11, double p11, double t12, double p12, double t13, double p13, double t14, double p14, double t15, double p15, int ss, double dudt) {   
    double uout;

    double ts[15];
    double us[15];

    double b[15];
    double c[15];
    double d[15];

    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    ts[10] = t11;
    ts[11] = t12;
    ts[12] = t13;
    ts[13] = t14;
    ts[14] = t15;

    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10;
    us[10] = p11;
    us[11] = p12;
    us[12] = p13;
    us[13] = p14;
    us[14] = p15;

    spline(15, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(15, t, ts, us, b, c, d);

    return(uout);
}

double spline20(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, double t11, double p11, double t12, double p12, double t13, double p13, double t14, double p14, double t15, double p15, double t16, double p16, double t17, double p17, double t18, double p18, double t19, double p19, double t20, double p20, int ss, double dudt) {   
    double uout;
    
    double ts[20];
    double us[20];
    
    double b[20];
    double c[20];
    double d[20];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    ts[10] = t11;
    ts[11] = t12;
    ts[12] = t13;
    ts[13] = t14;
    ts[14] = t15;
    ts[15] = t16;
    ts[16] = t17;
    ts[17] = t18;
    ts[18] = t19;
    ts[19] = t20;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10;
    us[10] = p11;
    us[11] = p12;
    us[12] = p13;
    us[13] = p14;
    us[14] = p15;
    us[15] = p16;
    us[16] = p17;
    us[17] = p18;
    us[18] = p19;
    us[19] = p20;
    
    spline(20, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(20, t, ts, us, b, c, d);
    
    return(uout);
}

double spline25(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, double t11, double p11, double t12, double p12, double t13, double p13, double t14, double p14, double t15, double p15, double t16, double p16, double t17, double p17, double t18, double p18, double t19, double p19, double t20, double p20, double t21, double p21, double t22, double p22, double t23, double p23, double t24, double p24, double t25, double p25, int ss, double dudt) {   
    double uout;
    
    double ts[25];
    double us[25];
    
    double b[25];
    double c[25];
    double d[25];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    ts[10] = t11;
    ts[11] = t12;
    ts[12] = t13;
    ts[13] = t14;
    ts[14] = t15;
    ts[15] = t16;
    ts[16] = t17;
    ts[17] = t18;
    ts[18] = t19;
    ts[19] = t20;
    ts[20] = t21;
    ts[21] = t22;
    ts[22] = t23;
    ts[23] = t24;
    ts[24] = t25;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10;
    us[10] = p11;
    us[11] = p12;
    us[12] = p13;
    us[13] = p14;
    us[14] = p15;
    us[15] = p16;
    us[16] = p17;
    us[17] = p18;
    us[18] = p19;
    us[19] = p20;
    us[20] = p21;
    us[21] = p22;
    us[22] = p23;
    us[23] = p24;
    us[24] = p25;
    
    spline(25, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(25, t, ts, us, b, c, d);
    
    return(uout);
}

double spline30(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, double t11, double p11, double t12, double p12, double t13, double p13, double t14, double p14, double t15, double p15, double t16, double p16, double t17, double p17, double t18, double p18, double t19, double p19, double t20, double p20, double t21, double p21, double t22, double p22, double t23, double p23, double t24, double p24, double t25, double p25, double t26, double p26, double t27, double p27, double t28, double p28, double t29, double p29, double t30, double p30, int ss, double dudt) {   
    double uout;
    
    double ts[30];
    double us[30];
    
    double b[30];
    double c[30];
    double d[30];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    ts[10] = t11;
    ts[11] = t12;
    ts[12] = t13;
    ts[13] = t14;
    ts[14] = t15;
    ts[15] = t16;
    ts[16] = t17;
    ts[17] = t18;
    ts[18] = t19;
    ts[19] = t20;
    ts[20] = t21;
    ts[21] = t22;
    ts[22] = t23;
    ts[23] = t24;
    ts[24] = t25;
    ts[25] = t26;
    ts[26] = t27;
    ts[27] = t28;
    ts[28] = t29;
    ts[29] = t30;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10;
    us[10] = p11;
    us[11] = p12;
    us[12] = p13;
    us[13] = p14;
    us[14] = p15;
    us[15] = p16;
    us[16] = p17;
    us[17] = p18;
    us[18] = p19;
    us[19] = p20;
    us[20] = p21;
    us[21] = p22;
    us[22] = p23;
    us[23] = p24;
    us[24] = p25;
    us[25] = p26;
    us[26] = p27;
    us[27] = p28;
    us[28] = p29;
    us[29] = p30;
    
    spline(30, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(30, t, ts, us, b, c, d);
    
    return(uout);
}

double monospline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10) {   
    double uout;
    
    double ts[10];
    double us[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10;
    
    monotoneSpline( 10, ts, us, b, c, d );
    uout = seval(10, t, ts, us, b, c, d);
    
    return(uout);
}

double spline_pos10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt) {   
    int is;
    double uout;
    
    double ts[10];
    double us[10];
    double uslog[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10; 
    
    for (is = 0; is<10; is++){
        uslog[is] = log(us[is]);
    }
    
    spline(10, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uout = seval(10, t, ts, uslog, b, c, d);
    
    return(exp(uout));
}

/* spline derivatives */
double Dspline3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt, int id) {   
    double uout;
    
    double ts[3];
    double us[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    
    us[id-1] = 1.0;
    
    spline(3, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(3, t, ts, us, b, c, d);
    
    return(uout);
}

double Dmonospline3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int id) {   
    double uout;
    
    double ts[3];
    double us[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    
    us[id-1] = 1.0;
    
    monotoneSpline( 3, ts, us, b, c, d );
    uout = seval(3, t, ts, us, b, c, d);
    
    return(uout);
}

double Dspline_pos3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt, int id) {   
    
    double uout;
    double uspline_pos3;
    double suspline3;
            
    double ps[3];
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    
    uspline_pos3 = spline_pos3(t, t1, p1, t2, p2, t3, p3, ss, dudt);
    suspline3 = Dspline3(t, t1, p1, t2, p2, t3, p3, ss, dudt, id);
    uout = suspline3 * uspline_pos3 * log(10.0);
    uout = uout / ps[id-1] / log(10.0);
    
    return(uout);
}

double Dspline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt, int id) {   
    double uout;
    
    double ts[4];
    double us[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    
    us[id-1] = 1.0;
    
    spline(4, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(4, t, ts, us, b, c, d);
    
    return(uout);
}

double Dmonospline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int id) {   
    double uout;
    
    double ts[4];
    double us[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    
    us[id-1] = 1.0;
    
    monotoneSpline( 4, ts, us, b, c, d );
    uout = seval(4, t, ts, us, b, c, d);
    
    return(uout);
}

double Dspline_pos4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt, int id) {   
    
    double uout;
    double uspline_pos4;
    double suspline4;
            
    double ps[4];
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    ps[3] = p4;
    
    uspline_pos4 = spline_pos4(t, t1, p1, t2, p2, t3, p3, t4, p4, ss, dudt);
    suspline4 = Dspline4(t, t1, p1, t2, p2, t3, p3, t4, p4, ss, dudt, id);
    uout = suspline4 * uspline_pos4 * log(10.0);
    uout = uout / ps[id-1] / log(10.0);
    
    return(uout);
}

double Dspline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt, int id) {   
    double uout;
    
    double ts[5];
    double us[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    
    us[id-1] = 1.0;
    
    spline(5, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(5, t, ts, us, b, c, d);
    
    return(uout);
}

double Dmonospline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int id) {   
    double uout;
    
    double ts[5];
    double us[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    
    us[id-1] = 1.0;
    
    monotoneSpline( 5, ts, us, b, c, d );
    uout = seval(5, t, ts, us, b, c, d);
    
    return(uout);
}

double Dspline_pos5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt, int id) {   
    
    double uout;
    double uspline_pos5;
    double suspline5;
            
    double ps[5];
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    ps[3] = p4;
    ps[4] = p5;
    
    uspline_pos5 = spline_pos5(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, ss, dudt);
    suspline5 = Dspline5(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, ss, dudt, id);
    uout = suspline5 * uspline_pos5 * log(10.0);
    uout = uout / ps[id-1] / log(10.0);
    
    return(uout);
}

double Dspline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt, int id) {   
    double uout;
    
    double ts[10];
    double us[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    us[5] = 0.0;
    us[6] = 0.0;
    us[7] = 0.0;
    us[8] = 0.0;
    us[9] = 0.0;
    
    us[id-1] = 1.0;
    
    spline(10, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(10, t, ts, us, b, c, d);
    
    return(uout);
}

double Dspline15(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, double t11, double p11, double t12, double p12, double t13, double p13, double t14, double p14, double t15, double p15, int ss, double dudt, int id) {   
    double uout;

    double ts[15];
    double us[15];

    double b[15];
    double c[15];
    double d[15];

    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    ts[10] = t11;
    ts[11] = t12;
    ts[12] = t13;
    ts[13] = t14;
    ts[14] = t15;

    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    us[5] = 0.0;
    us[6] = 0.0;
    us[7] = 0.0;
    us[8] = 0.0;
    us[9] = 0.0;
    us[10] = 0.0;
    us[11] = 0.0;
    us[12] = 0.0;
    us[13] = 0.0;
    us[14] = 0.0;

    us[id-1] = 1.0;

    spline(15, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(15, t, ts, us, b, c, d);

    return(uout);
}

double Dspline20(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, double t11, double p11, double t12, double p12, double t13, double p13, double t14, double p14, double t15, double p15, double t16, double p16, double t17, double p17, double t18, double p18, double t19, double p19, double t20, double p20, int ss, double dudt, int id) {   
    double uout;
    
    double ts[20];
    double us[20];
    
    double b[20];
    double c[20];
    double d[20];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    ts[10] = t11;
    ts[11] = t12;
    ts[12] = t13;
    ts[13] = t14;
    ts[14] = t15;
    ts[15] = t16;
    ts[16] = t17;
    ts[17] = t18;
    ts[18] = t19;
    ts[19] = t20;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    us[5] = 0.0;
    us[6] = 0.0;
    us[7] = 0.0;
    us[8] = 0.0;
    us[9] = 0.0;
    us[10] = 0.0;
    us[11] = 0.0;
    us[12] = 0.0;
    us[13] = 0.0;
    us[14] = 0.0;
    us[15] = 0.0;
    us[16] = 0.0;
    us[17] = 0.0;
    us[18] = 0.0;
    us[19] = 0.0;
    
    us[id-1] = 1.0;
    
    spline(20, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(20, t, ts, us, b, c, d);
    
    return(uout);
}

double Dspline25(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, double t11, double p11, double t12, double p12, double t13, double p13, double t14, double p14, double t15, double p15, double t16, double p16, double t17, double p17, double t18, double p18, double t19, double p19, double t20, double p20, double t21, double p21, double t22, double p22, double t23, double p23, double t24, double p24, double t25, double p25, int ss, double dudt, int id) {   
    double uout;
    
    double ts[25];
    double us[25];
    
    double b[25];
    double c[25];
    double d[25];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    ts[10] = t11;
    ts[11] = t12;
    ts[12] = t13;
    ts[13] = t14;
    ts[14] = t15;
    ts[15] = t16;
    ts[16] = t17;
    ts[17] = t18;
    ts[18] = t19;
    ts[19] = t20;
    ts[20] = t21;
    ts[21] = t22;
    ts[22] = t23;
    ts[23] = t24;
    ts[24] = t25;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    us[5] = 0.0;
    us[6] = 0.0;
    us[7] = 0.0;
    us[8] = 0.0;
    us[9] = 0.0;
    us[10] = 0.0;
    us[11] = 0.0;
    us[12] = 0.0;
    us[13] = 0.0;
    us[14] = 0.0;
    us[15] = 0.0;
    us[16] = 0.0;
    us[17] = 0.0;
    us[18] = 0.0;
    us[19] = 0.0;
    us[20] = 0.0;
    us[21] = 0.0;
    us[22] = 0.0;
    us[23] = 0.0;
    us[24] = 0.0;

    us[id-1] = 1.0;
    
    spline(25, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(25, t, ts, us, b, c, d);
    
    return(uout);
}

double Dspline30(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, double t11, double p11, double t12, double p12, double t13, double p13, double t14, double p14, double t15, double p15, double t16, double p16, double t17, double p17, double t18, double p18, double t19, double p19, double t20, double p20, double t21, double p21, double t22, double p22, double t23, double p23, double t24, double p24, double t25, double p25, double t26, double p26, double t27, double p27, double t28, double p28, double t29, double p29, double t30, double p30, int ss, double dudt, int id) {   
    double uout;
    
    double ts[30];
    double us[30];
    
    double b[30];
    double c[30];
    double d[30];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    ts[10] = t11;
    ts[11] = t12;
    ts[12] = t13;
    ts[13] = t14;
    ts[14] = t15;
    ts[15] = t16;
    ts[16] = t17;
    ts[17] = t18;
    ts[18] = t19;
    ts[19] = t20;
    ts[20] = t21;
    ts[21] = t22;
    ts[22] = t23;
    ts[23] = t24;
    ts[24] = t25;
    ts[25] = t26;
    ts[26] = t27;
    ts[27] = t28;
    ts[28] = t29;
    ts[29] = t30;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    us[5] = 0.0;
    us[6] = 0.0;
    us[7] = 0.0;
    us[8] = 0.0;
    us[9] = 0.0;
    us[10] = 0.0;
    us[11] = 0.0;
    us[12] = 0.0;
    us[13] = 0.0;
    us[14] = 0.0;
    us[15] = 0.0;
    us[16] = 0.0;
    us[17] = 0.0;
    us[18] = 0.0;
    us[19] = 0.0;
    us[20] = 0.0;
    us[21] = 0.0;
    us[22] = 0.0;
    us[23] = 0.0;
    us[24] = 0.0;
    us[25] = 0.0;
    us[26] = 0.0;
    us[27] = 0.0;
    us[28] = 0.0;
    us[29] = 0.0;

    us[id-1] = 1.0;
    
    spline(30, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(30, t, ts, us, b, c, d);
    
    return(uout);
}

double Dmonospline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int id) {   
    double uout;
    
    double ts[10];
    double us[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    us[5] = 0.0;
    us[6] = 0.0;
    us[7] = 0.0;
    us[8] = 0.0;
    us[9] = 0.0;
    
    us[id-1] = 1.0;
    
    monotoneSpline( 10, ts, us, b, c, d );
    uout = seval(10, t, ts, us, b, c, d);
    
    return(uout);
}

double Dspline_pos10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt, int id) {   
    
    double uout;
    double uspline_pos10;
    double suspline10;
            
    double ps[10];
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    ps[3] = p4;
    ps[4] = p5;
    ps[5] = p6;
    ps[6] = p7;
    ps[7] = p8;
    ps[8] = p9;
    ps[9] = p10;
    
    uspline_pos10 = spline_pos10(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9, t10, p10, ss, dudt);
    suspline10 = Dspline10(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9, t10, p10, ss, dudt, id);
    uout = suspline10 * uspline_pos10 * log(10.0);
    uout = uout / ps[id-1] / log(10.0);
    
    return(uout);
}


/* Faster implementation of the splines */
/* This version caches the spline in a userdata struct so that the coefficients don't have to be determined every RHS evaluation */

/* Function which can cache the computed spline */
int cspline(int n, int end1, int end2, double slope1, double slope2, const double x[], const double y[], double b[], double c[], double d[], int cacheID, double **splineCache, int *IDcache)
{
    int j;
    
    /* Nothing in the cache? Allocate memory and fill (free-ing is handled in arSimuCalc.c) */
    if ( splineCache[cacheID] == NULL )
    {
        splineCache[cacheID] = (double*) malloc(3 * n * sizeof(double));
        spline(n, end1, end2, slope1, slope2, x, y, b, c, d);
        for ( j = 0; j < n; j++ )
        {
            splineCache[cacheID][j]     = b[j];
            splineCache[cacheID][j+n]   = c[j];
            splineCache[cacheID][j+2*n] = d[j];
        }
        IDcache[cacheID] = 0;
    }
    else
    {
        for ( j = 0; j < n; j++ )
        {
            b[j] = splineCache[cacheID][j];
            c[j] = splineCache[cacheID][j+n];
            d[j] = splineCache[cacheID][j+2*n];
        }        
    }
}

int cmonotoneSpline( int n, const double x[], const double y[], double b[], double c[], double d[], int cacheID, double **splineCache, int *IDcache )
{
    int j;
    
    /* Nothing in the cache? Allocate memory and fill (free-ing is handled in arSimuCalc.c) */
    if ( splineCache[cacheID] == NULL )
    {
        splineCache[cacheID] = (double*) malloc(3 * n * sizeof(double));
        monotoneSpline(n, x, y, b, c, d);
        for ( j = 0; j < n; j++ )
        {
            splineCache[cacheID][j]     = b[j];
            splineCache[cacheID][j+n]   = c[j];
            splineCache[cacheID][j+2*n] = d[j];
        }
        IDcache[cacheID] = 0;
    }
    else
    {
        for ( j = 0; j < n; j++ )
        {
            b[j] = splineCache[cacheID][j];
            c[j] = splineCache[cacheID][j+n];
            d[j] = splineCache[cacheID][j+2*n];
        }        
    } 
}

int clongmonotoneSpline( int n, const double x[], const double y[], double b[], double c[], double d[], int cacheID, double **splineCache, int *IDcache )
{
    int j;
    
    /* Nothing in the cache? Allocate memory and fill (free-ing is handled in arSimuCalc.c) */
    if ( splineCache[cacheID] == NULL )
    {
        splineCache[cacheID] = (double*) malloc(3 * n * sizeof(double));
        longMonotoneSpline(n, x, y, b, c, d);
        for ( j = 0; j < n; j++ )
        {
            splineCache[cacheID][j]     = b[j];
            splineCache[cacheID][j+n]   = c[j];
            splineCache[cacheID][j+2*n] = d[j];
        }
        IDcache[cacheID] = 0;
    }
    else
    {
        for ( j = 0; j < n; j++ )
        {
            b[j] = splineCache[cacheID][j];
            c[j] = splineCache[cacheID][j+n];
            d[j] = splineCache[cacheID][j+2*n];
        }        
    } 
}

double fastspline3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt) {   
    double uout;
    
    double ts[3];
    double us[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    
    cspline(3, ss, 0, dudt, 0.0, ts, us, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(3, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double monofastspline3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3) {   
    double uout;
    
    double ts[3];
    double us[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;

    cmonotoneSpline( 3, ts, us, b, c, d, ID, splineCache, idCache );
    uout = seval_fixed(3, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double fastspline_pos3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt) {   
    int is;
    double uout;
    
    double ts[3];
    double us[3];
    double uslog[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    
    for (is = 0; is<3; is++){
        uslog[is] = log(us[is]);
    }
    
    cspline(3, ss, 0, dudt, 0.0, ts, uslog, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(3, t, ts, uslog, b, c, d, &(idCache[ID]));
    
    return(exp(uout));
}

double fastspline4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt) {   
    double uout;
    
    double ts[4];
    double us[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    
    cspline(4, ss, 0, dudt, 0.0, ts, us, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(4, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double monofastspline4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4) {   
    double uout;
    
    double ts[4];
    double us[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;

    cmonotoneSpline( 4, ts, us, b, c, d, ID, splineCache, idCache );
    uout = seval_fixed(4, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double fastspline_pos4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt) {   
    int is;
    double uout;
    
    double ts[4];
    double us[4];
    double uslog[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    
    for (is = 0; is<4; is++){
        uslog[is] = log(us[is]);
    }
    
    cspline(4, ss, 0, dudt, 0.0, ts, uslog, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(4, t, ts, uslog, b, c, d, &(idCache[ID]));
    
    return(exp(uout));
}

double fastspline5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt) {   
    double uout;
    
    double ts[5];
    double us[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    
    cspline(5, ss, 0, dudt, 0.0, ts, us, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(5, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double monofastspline5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5) {   
    double uout;
    
    double ts[5];
    double us[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;

    cmonotoneSpline( 5, ts, us, b, c, d, ID, splineCache, idCache );
    uout = seval_fixed(5, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double fastspline_pos5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt) {   
    int is;
    double uout;
    
    double ts[5];
    double us[5];
    double uslog[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    
    for (is = 0; is<5; is++){
        uslog[is] = log(us[is]);
    }
    
    cspline(5, ss, 0, dudt, 0.0, ts, uslog, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(5, t, ts, uslog, b, c, d, &(idCache[ID]));
    
    return(exp(uout));
}

double fastspline10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt) {   
    double uout;
    
    double ts[10];
    double us[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10;
    
    cspline(10, ss, 0, dudt, 0.0, ts, us, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(10, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double monofastspline10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10) {   
    double uout;
    
    double ts[10];
    double us[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10;
    
    cmonotoneSpline( 10, ts, us, b, c, d, ID, splineCache, idCache );
    uout = seval_fixed(10, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double fastspline_pos10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt) {   
    int is;
    double uout;
    
    double ts[10];
    double us[10];
    double uslog[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10; 
    
    for (is = 0; is<10; is++){
        uslog[is] = log(us[is]);
    }
    
    cspline(10, ss, 0, dudt, 0.0, ts, uslog, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(10, t, ts, uslog, b, c, d, &(idCache[ID]));
    
    return(exp(uout));
}

/* fast spline derivatives */
double Dfastspline3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt, int id) {   
    double uout;
    
    double ts[3];
    double us[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    
    us[id-1] = 1.0;
    
    cspline(3, ss, 0, dudt, 0.0, ts, us, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(3, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double Dmonofastspline3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, int id) {   
    double uout;
    
    double ts[3];
    double us[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    
    us[id-1] = 1.0;
    
    cmonotoneSpline( 3, ts, us, b, c, d, ID, splineCache, idCache );
    uout = seval_fixed(3, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double Dfastspline_pos3(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt, int id) {   
    
    double uout;
    double uspline_pos3;
    double suspline3;
            
    double ps[3];
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    
    uspline_pos3 = spline_pos3(t, t1, p1, t2, p2, t3, p3, ss, dudt);
    suspline3 = Dspline3(t, t1, p1, t2, p2, t3, p3, ss, dudt, id);
    uout = suspline3 * uspline_pos3 * log(10.0);
    uout = uout / ps[id-1] / log(10.0);
    
    return(uout);
}

double Dfastspline4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt, int id) {   
    double uout;
    
    double ts[4];
    double us[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    
    us[id-1] = 1.0;
    
    cspline(4, ss, 0, dudt, 0.0, ts, us, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(4, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double Dmonofastspline4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int id) {   
    double uout;
    
    double ts[4];
    double us[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    
    us[id-1] = 1.0;
    
    cmonotoneSpline( 4, ts, us, b, c, d, ID, splineCache, idCache );
    uout = seval_fixed(4, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double Dfastspline_pos4(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt, int id) {   
    
    double uout;
    double uspline_pos4;
    double suspline4;
            
    double ps[4];
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    ps[3] = p4;
    
    uspline_pos4 = spline_pos4(t, t1, p1, t2, p2, t3, p3, t4, p4, ss, dudt);
    suspline4 = Dspline4(t, t1, p1, t2, p2, t3, p3, t4, p4, ss, dudt, id);
    uout = suspline4 * uspline_pos4 * log(10.0);
    uout = uout / ps[id-1] / log(10.0);
    
    return(uout);
}

double Dfastspline5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt, int id) {   
    double uout;
    
    double ts[5];
    double us[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    
    us[id-1] = 1.0;
    
    cspline(5, ss, 0, dudt, 0.0, ts, us, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(5, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double Dmonofastspline5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int id) {   
    double uout;
    
    double ts[5];
    double us[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    
    us[id-1] = 1.0;
    
    cmonotoneSpline( 5, ts, us, b, c, d, ID, splineCache, idCache );
    uout = seval_fixed(5, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double Dfastspline_pos5(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt, int id) {   
    
    double uout;
    double uspline_pos5;
    double suspline5;
            
    double ps[5];
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    ps[3] = p4;
    ps[4] = p5;
    
    uspline_pos5 = spline_pos5(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, ss, dudt);
    suspline5 = Dspline5(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, ss, dudt, id);
    uout = suspline5 * uspline_pos5 * log(10.0);
    uout = uout / ps[id-1] / log(10.0);
    
    return(uout);
}

double Dfastspline10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt, int id) {   
    double uout;
    
    double ts[10];
    double us[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    us[5] = 0.0;
    us[6] = 0.0;
    us[7] = 0.0;
    us[8] = 0.0;
    us[9] = 0.0;
    
    us[id-1] = 1.0;
    
    cspline(10, ss, 0, dudt, 0.0, ts, us, b, c, d, ID, splineCache, idCache);
    uout = seval_fixed(10, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double Dmonofastspline10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int id) {   
    double uout;
    
    double ts[10];
    double us[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    us[5] = 0.0;
    us[6] = 0.0;
    us[7] = 0.0;
    us[8] = 0.0;
    us[9] = 0.0;
    
    us[id-1] = 1.0;
    
    cmonotoneSpline( 10, ts, us, b, c, d, ID, splineCache, idCache );
    uout = seval_fixed(10, t, ts, us, b, c, d, &(idCache[ID]));
    
    return(uout);
}

double Dfastspline_pos10(double t, int ID, double **splineCache, int *idCache, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt, int id) {   
    
    double uout;
    double uspline_pos10;
    double suspline10;
            
    double ps[10];
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    ps[3] = p4;
    ps[4] = p5;
    ps[5] = p6;
    ps[6] = p7;
    ps[7] = p8;
    ps[8] = p9;
    ps[9] = p10;
    
    uspline_pos10 = spline_pos10(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9, t10, p10, ss, dudt);
    suspline10 = Dspline10(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9, t10, p10, ss, dudt, id);
    uout = suspline10 * uspline_pos10 * log(10.0);
    uout = uout / ps[id-1] / log(10.0);
    
    return(uout);
}
    
/* custom rate laws */

double mmenten(double x, double vmax, double km){
    return(vmax * x / (km + x));
}

double mmenten_alt(double x, double klin, double ksat){
    return(0.0);
}

double hill_kd(double x, double h, double kd){
    return(pow(x,h) / (kd + pow(x,h)));
}

double hill_ka(double x, double h, double ka){
    /* return(pow(x,h) / (pow(ka,h) + pow(x,h))); */
    return(1 / (pow((ka/x),h) + 1));
}


