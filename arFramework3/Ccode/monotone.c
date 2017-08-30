/* Implementation of a monotonic spline                                               */
/* Contact: joep.vanlier@fdm.uni-freiburg.de                                          */
/* Function to compute spline coefficients of a spline that will always be monotonous */
/* Usage example: monotoneSpline( 10, ts, uslog, b, c, d )                            */

#include <stdlib.h>

/* TO DOs: */
/* There is a reason the long splines use a fixed define for the length. 
 * Allocating this memory on the heap would be slower. The define is undesirable 
 * and should ideally be replaced with a compilation-time generated number for the 
 * maximal spline length. Another potential large performance gain could be achieved 
 * by keeping dx, dy and ms in static memory; but this would require passing a thread 
 * ID, as every thread should get its own static variable to avoid dangerous race 
 * conditions */

#define MAX_LONG_SPLINE 10000

int longMonotoneSpline( const int n, const double x[], const double y[], double b[], double c[], double d[] )
{
    int i;
    double coeff;
    double m;
    double invDx;
    double common;
    
    /* For long splines, MAX_LONG_SPLINE is the maximum length */
    double dy[MAX_LONG_SPLINE];
    double dx[MAX_LONG_SPLINE];
    double ms[MAX_LONG_SPLINE];
    
    /* Slopes and differences */
    for ( i = 0; i < (n-1); i++ )
    {
        dx[i] = x[i+1] - x[i];
        dy[i] = y[i+1] - y[i];
        ms[i] = dy[i]  / dx[i];
    }
    
    /* Degree one coefficients */
    b[0] = ms[0];
    for ( i = 0; i < (n-2); i++ )
    {
        if ( ms[i]*ms[i+1] <= 0 )
            b[i+1] = 0;
        else
        {
            common = dx[i] + dx[i+1];
            b[i+1] = (3*common/((common + dx[i+1])/ms[i] + (common+dx[i])/(ms[i+1])));
        }
    }
    b[n-1] = ms[n-2];
    
	/* Degree two and three coefficients */
	for ( i = 0; i < (n-1); i++ )
    {
        coeff    = b[i];
        m        = ms[i];
        invDx    = 1/dx[i];
        common   = coeff + b[i + 1] - m - m;
        
        c[i] = (m-coeff-common)*invDx;
        d[i] = (common*invDx*invDx);
	}
    c[n-1] = 0;
    d[n-1] = 0;
        
    return 1;
}

int monotoneSpline( const int n, const double x[], const double y[], double b[], double c[], double d[] )
{
    int i;
    double coeff;
    double m;
    double invDx;
    double common;
    
    /* For normal splines, 10 is the maximum length. */
    double dy[10];
    double dx[10];
    double ms[10];
    
    /* Slopes and differences */
    for ( i = 0; i < (n-1); i++ )
    {
        dx[i] = x[i+1] - x[i];
        dy[i] = y[i+1] - y[i];
        ms[i] = dy[i]  / dx[i];
    }
    
    /* Degree one coefficients */
    b[0] = ms[0];
    for ( i = 0; i < (n-2); i++ )
    {
        if ( ms[i]*ms[i+1] <= 0 )
            b[i+1] = 0;
        else
        {
            common = dx[i] + dx[i+1];
            b[i+1] = (3*common/((common + dx[i+1])/ms[i] + (common+dx[i])/(ms[i+1])));
        }
    }
    b[n-1] = ms[n-2];
    
	/* Degree two and three coefficients */
	for ( i = 0; i < (n-1); i++ )
    {
        coeff    = b[i];
        m        = ms[i];
        invDx    = 1/dx[i];
        common   = coeff + b[i + 1] - m - m;
        
        c[i] = (m-coeff-common)*invDx;
        d[i] = (common*invDx*invDx);
	}
    c[n-1] = 0;
    d[n-1] = 0;
        
    return 1;
}
