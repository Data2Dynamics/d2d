#include "inverseC.h"
#include "blas.h"
#include "lapack.h"

#define DEBUGPRINT0(DBGMODE, LVL, STR) { if ( DBGMODE > LVL ) { mexPrintf( " [D] " ); mexPrintf( STR ); mexEvalString("drawnow;"); } }
#define DEBUGPRINT1(DBGMODE, LVL, STR, ARG1) { if ( DBGMODE > LVL ) { mexPrintf( " [D] " ); mexPrintf( STR, ARG1 ); mexEvalString("drawnow;"); } }
#define DEBUGPRINT2(DBGMODE, LVL, STR, ARG1, ARG2) { if ( DBGMODE > LVL ) { mexPrintf( " [D] " ); mexPrintf( STR, ARG1, ARG2 ); mexEvalString("drawnow;"); } }
#define DEBUGPRINT3(DBGMODE, LVL, STR, ARG1, ARG2, ARG3) { if ( DBGMODE > LVL ) { mexPrintf( " [D] " ); mexPrintf( STR, ARG1, ARG2, ARG3 ); mexEvalString("drawnow;"); } }
#define DEBUGPRINT4(DBGMODE, LVL, STR, ARG1, ARG2, ARG3, ARG4) { if ( DBGMODE > LVL ) { mexPrintf( " [D] " ); mexPrintf( STR, ARG1, ARG2, ARG3, ARG4 ); mexEvalString("drawnow;"); } }

/* Invert a matrix using LAPACK */
void invert(double* mat, mwSignedIndex N)
{
    mwSignedIndex *ipiv;
    double *workmem;
    mwSignedIndex workSize = N*N;
    mwSignedIndex errorCode;
    
    ipiv = (mwSignedIndex *) malloc(sizeof(mwSignedIndex) * (N+1));
    workmem = (double *) malloc( sizeof(double) * workSize );

    dgetrf(&N,&N,mat,&N,ipiv,&errorCode);
    dgetri(&N,mat,&N,ipiv,workmem,&workSize,&errorCode);

    free(ipiv);
    free(workmem);
};

/* Matrix multiply for the rootfinding using BLAS */
void mmultiply(double* invdfdx, double* f, double *result, mwSignedIndex N)
{
    char *chn1 = "N";
    char *chn2 = "N";
    mwSignedIndex intone = 1;
    double one = 1.0;
    double zero = 0.0;
    
    dgemm(chn1, chn2, &N, &intone, &N, &one, invdfdx, &N, f, &N, &zero, result, &N);
};

/* Single Newton-Raphson step */
double NRstep( double *xptr, double* dfdx, double* f, double *tmp, int N )
{
    double maxAbs;
    int i;
           
    invert(dfdx, N);
    mmultiply( dfdx, f, tmp, N );

    maxAbs = 0.0;
    for ( i = 0; i < N; i++ )
    {
        xptr[i] = xptr[i] - tmp[i];
        maxAbs = ( maxAbs > fabs(tmp[i]) ) ? maxAbs : fabs(tmp[i]);
    }
    
    return maxAbs;
};

/* Simple rootfinding function */
double solveSS( int debugMode, mxArray *arcondition, int im, int ic, int isim, double t, N_Vector x, void *user_data, double tol )
{
    int i;
    int N;
    int iterationLimit;
    double maxAbs;
    double *xptr;
	double *dfdx;
    double *f;
    double *tmp;
    
    xptr            = N_VGetArrayPointer(x);
    dfdx            = mxGetData(mxGetField(arcondition, ic, "dfdxNum"));
    f               = mxGetData(mxGetField(arcondition, ic, "dxdt"));
    N               = mxGetNumberOfElements(mxGetField(arcondition, ic, "dxdt"));
    tmp             = (double *)malloc( sizeof( double ) * N );
    iterationLimit  = 100;
    maxAbs          = 100000000.0;
    
    /* Perform Newton Raphson until convergence or abort */
    i = 0;
    while( ( maxAbs > tol ) && ( i < iterationLimit ) )
    {
        DEBUGPRINT3( debugMode, 9, "Rootfinding iteration %d: Maximum absolute dx: %g (tolerance = %g)\n", i, maxAbs, tol );
        fx( t, x, f, user_data, im, isim );                           /* Compute RHS and store in f */
        if ( debugMode > 10 )
        {
            int j;
            for ( j = 0; j < N; j++ )
                DEBUGPRINT2( debugMode, 10, "fx[%d] = %g", j, f[j] );
            
            for ( j = 0; j < N; j++ )
                DEBUGPRINT2( debugMode, 10, "x[%d] = %g", j, x[j] );            
        }
        DEBUGPRINT0( debugMode, 10, "postfx\n" );
        getdfxdx( im, isim, t, x, dfdx, user_data );                  /* Updates dvdx and stores dfxdx */
        DEBUGPRINT0( debugMode, 10, "postgetdfdx\n" );
        maxAbs = NRstep( xptr, dfdx, f, tmp, N );                     /* Perform Newton-Raphson step */
        DEBUGPRINT0( debugMode, 10, "postNRstep\n" );
        i++;
    }
    DEBUGPRINT3( debugMode, 9, "Final absolute at interation %d = dx: %g (tolerance = %g)\n", i, maxAbs, tol );
    
    free( tmp );
    
    return maxAbs;
}

