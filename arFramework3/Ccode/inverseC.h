#include <cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDENSE fcts. and constants */
#include <cvodes/cvodes_sparse.h>     /* prototype for CVSPARSE fcts. and constants */
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */
#include <mex.h>

/* Simple rootfinding function */
double solveSS( int debugMode, mxArray *arcondition, int im, int ic, int isim, double t, N_Vector x, void *user_data, double tol );

/* Single Newton-Raphson step */
double NRstep( double *xptr, double* dfdx, double* f, double *tmp, mwSignedIndex workSize, double *workmem, mwSignedIndex *ipiv, int N );

/* Invert a matrix using LAPACK */
void invert(double* mat, mwSignedIndex workSize, double *workmem, mwSignedIndex *ipiv, mwSignedIndex N);

/* Matrix multiply for the rootfinding using BLAS */
void mmultiply(double* invdfdx, double* f, double *result, mwSignedIndex N);
