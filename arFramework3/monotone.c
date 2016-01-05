/* Implementation of a monotonic spline                                               */
/* Contact: joep.vanlier@fdm.uni-freiburg.de                                          */
/* Function to compute spline coefficients of a spline that will always be monotonous */
/* Usage example: monotoneSpline( 10, ts, uslog, b, c, d )                            */

#include <stdlib.h>

int monotoneSpline( int n, double x[], double y[], double b[], double c[], double d[] )
{
    int i;
    double coeff;
    double m;
    double invDx;
    double common;
    
    /* Since we have a maximum number of 10 knots, might as well put it on the stack */
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
