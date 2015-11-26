/* MEX interface for NL2SOL             */
/* (C) Joep Vanlier 2015                */
/* Contact: joep.vanlier at gmail.com   */

#include "mex.h"
#include <string.h> /* For memcpy */

#define parFUN              prhs[0]
#define parX0               prhs[1]
#define parLB               prhs[2]
#define parUB               prhs[3]
#define parOPTS             prhs[4]
#define parPRINTLEVEL       prhs[5]

double getStatus(int stat, int printLevel);
void validateInput( const mxArray *prhs[], int nrhs, int *npars, int *bounded );
void loadOptions( const mxArray *prhs[], int* maxFun, int* maxIter, double* atol, double* rtol, double* xtol, int* printLevel, int* useJacobian );
static void CALCR(int *n, int *p, double *x, int *nf, double *r, int *uiparm, double *ydata, void *ufparm);
static void CALCJ(int *n, int *p, double *x, int *nf, double *j, int *uiparm, double *ydata, void *ufparm);

typedef struct {
     char f[128], g[128];
     mxArray *plhs[2];
     mxArray *callData[2];
     int xrhs, nrhs;
     int nJac, nFunc;
} UserData;

/* Linux FORTRAN calling convention */
#ifdef __linux__
    #define dn2g    dn2g_
    #define dn2gb   dn2gb_
    #define dn2f    dn2f_
    #define dn2fb   dn2fb_
    #define divset  divset_
    #define d1mach  d1mach_
#endif

/* Windows FORTRAN calling convention */
#ifdef _WIN32
    #define dn2g    DN2G
    #define dn2gb   DN2GB
    #define dn2f    DN2F
    #define dn2fb   DN2FB
    #define divset  DIVSET
    #define d1mach  D1MACH
#endif

/* Bounded NL2SOL routine linkages */
extern void dn2g(int *n, int *p, double *x, 
                 void(*fun)(int*,int*,double*,int*,double*,int*,double*,void*),
                 void(*jac)(int*,int*,double*,int*,double*,int*,double*,void*),
                 int *iv, int *liv, int *lv, double *v, int *uiparm, double *urparm, void *ufparm);

extern void dn2gb(int *n, int *p, double *x, double *b,
                 void(*fun)(int*,int*,double*,int*,double*,int*,double*,void*),
                 void(*jac)(int*,int*,double*,int*,double*,int*,double*,void*),
                 int *iv, int *liv, int *lv, double *v, int *uiparm, double *urparm, void *ufparm);

extern void dn2f(int *n, int *p, double *x, 
                 void(*fun)(int*,int*,double*,int*,double*,int*,double*,void*),
                 int *iv, int *liv, int *lv, double *v, int *uiparm, double *urparm, void *ufparm);

extern void dn2fb(int *n, int *p, double *x, double *b,
                 void(*fun)(int*,int*,double*,int*,double*,int*,double*,void*),
                 int *iv, int *liv, int *lv, double *v, int *uiparm, double *urparm, void *ufparm);

/* Default parameters for NL2SOL */
extern void divset(int *alg, int *iv, int *liv, int *lv, double *v);

/* Function to get certain number representations for specific machines */
extern double d1mach(int *val);

/* X = lsqnonlin(FUN,X0,LB,UB,OPTIONS) */
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    int     bounded;
    int     useJacobian;
    
    UserData userData;
    
    /* Temp variables to test residual function */
    mxArray *lhs[2];
    int     result;
    
    /* Iterators */
    int     i;
    int     j;
    
    /* NL2SOL storage / info */
    int     n;                   /* Length residuals             */   
    int     p;                   /* Length parameters            */
    int     *iv;                 /* Storage                      */
    double  *v;                  /* More storage                 */
    int     liv;                 /* Length of IV                 */
    int     lv;                  /* Length of V                  */
    int     uiparm[2];           /* User parameters              */    
    double  *bounds;             /* Bounds                       */
    int     one;
    int     two;
    
    /* Solver input */
    double  *x0;                 /* Initial parameter vector     */
    double  *yData;              /* Vector with initial yData    */
    double  *lb;                 /* Lower bound                  */
    double  *ub;                 /* Upper bound                  */
    int     printLevel;          /* User specified printlevel    */
    int     maxFun;              /* Maximum number of evals      */
    int     maxIter;             /* Maximum number of iterations */
    double  atol;                /* Absolute tolerance           */
    double  rtol;                /* Relative tolerance           */
    double  xtol;                /* X-tolerance                  */
    
    /* Solver output */
    double  *x;                  /* Output                       */
    double  *fVal;               /* Final function value         */
    double  *fEval;              /* Number of function evals     */
    double  *iter;               /* Iteration count              */
    double  *exitFlag;           /* Exit flag                    */
    double  *lambda;             /* Lambda                       */
    
    /* Defaults */
    maxFun          = 1000;
    maxIter         = 1000;
    atol            = 1e-7;
    rtol            = 1e-6;
    xtol            = -1.0;
    printLevel      = 0;
    useJacobian     = 0;
    
    /* Integers so we can pass them by address to FORTRAN */
    one             = 1;
    two             = 2;
       
    /* Validate input, check whether we have bounds and determine number of parameters */
    validateInput( prhs, nrhs, &p, &bounded );

    /* Options struct */
    if ( nrhs > 3 )
        loadOptions( prhs, &maxFun, &maxIter, &atol, &rtol, &xtol, &printLevel, &useJacobian );

    /* Fetch print level */
    if ( nrhs > 5 )
    {
        /* Temporary use x0 as temp storage */
        x0 = mxGetPr(parPRINTLEVEL);
        printLevel = (int) x0[0];
    }
    
    uiparm[0]       = 1;
    uiparm[1]       = printLevel;
    
    /* Fetch data */
    x0 = mxGetPr(parX0);
    
    bounds = NULL;
    if ( bounded == 1 )
    {
        lb      = mxGetPr(parLB);
        ub      = mxGetPr(parUB);
        bounds  = (double*) mxCalloc( (mwSize) 2 * p, sizeof(double) );
        
        /* Copy the bounds in the format NL2SOL expects */
        for ( i=0; i < p; i++ )
        {
            if ( mxIsInf(lb[j]) )
                bounds[2*i] = -d1mach(&two);
            else
                bounds[2*i] = lb[i];
            
            if ( mxIsInf(ub[j]) )
                bounds[2*i+1] = d1mach(&two);
            else
                bounds[2*i+1] = ub[i];            
        }
    }
       
    /* Allocate function memory */
    userData.callData[1] = mxCreateDoubleMatrix(p, 1, mxREAL);
    userData.nJac = 0;
    userData.nFunc = 0;
    
    /* Grab function and prepare function memory */
    userData.callData[0] = (mxArray*) parFUN;
    
    /* Insert initial parameters */
    memcpy( mxGetPr(userData.callData[1]), x0, p*sizeof(double) );
    userData.nrhs = 2;
    
    /* Evaluate function to get length of residuals and see if it is sensitivity capable */
    if ( useJacobian == 0 )
        result = mexCallMATLAB( 1, lhs, userData.nrhs, userData.callData, "feval" );
    else
        result = mexCallMATLAB( 2, lhs, userData.nrhs, userData.callData, "feval" );
    
    if ( result )
        mexErrMsgTxt("Error calling objective function");
    
    /* Fetch residual (to check the size) */
    if(!mxIsDouble(lhs[0]) || mxIsComplex(lhs[0]) || mxIsEmpty(lhs[0]))
        mexErrMsgTxt("Error: Objective function did not return residual vector");
    else
        n = (int) mxGetNumberOfElements( lhs[0] );
    
    /* Check size of the returned Jacobian */
    if ( useJacobian == 1 )
    {
        if(!mxIsDouble(lhs[1]) || mxIsComplex(lhs[1]) || mxIsEmpty(lhs[1]))
        {
            mexPrintf("Error: Objective function did not return Jacobian. Switching to finite differencing!");
            useJacobian = 0;
        } else
        {
            /* Check if Jacobian is the correct size before we segfault all over the place */
            if ( !( (mxGetN( lhs[1] ) == p) && ( mxGetM( lhs[1] ) == n ) ) )
            {
                mexPrintf("Error: Jacobian wrong size. Switching to finite differencing!");
                useJacobian = 0;
            } else
                yData = (double*) mxGetPr(lhs[1]);
        }
    }
       
    /* Grab data memory */
    yData = (double*) mxGetPr(lhs[0]);
    
    /* Parameter output handle */
    plhs[0]     = mxCreateDoubleMatrix(p,1, mxREAL);
    x           = mxGetPr( plhs[0] );
       
    /* Copy initial guess in case everything fails */
    memcpy(x, x0, p*sizeof(double));

    /* Start assigning the solver options and allocate its storage */
    if ( bounded ) {
        liv = 82 + 4*p;
        lv  = 105 + p*(n + 2*p + 21) + 2*n;
    } else {
        liv = 82 + p;
        lv  = 105 + p*(n + 2*p + 17) + 2*n;
    }
    
    iv  = NULL;
    v   = NULL;
    iv  = (int*) mxCalloc( liv + 10, sizeof(int) );
    v   = (double*) mxCalloc( lv + 10, sizeof(double) );
    
    /* Fetch NLSOL defaults */
    divset(&one, iv, &liv, &lv, v);
    
    /* Fill our own custom functions */
    iv[13] = 0;         /* no covariance */
    iv[14] = 0;         /* no covariance */
    iv[16] = maxFun;    /* Max function and gradient evaluations */
    iv[17] = maxIter;   /* Maximum number of iterations */
    iv[18] = 0;         /* No iteration printing */
    iv[19] = 0;         /* Default printing */
    iv[20] = 0;         /* Output unit printing */
    iv[21] = 0;         /* X printing */
    iv[22] = 0;         /* Summary printing */
    iv[23] = 0;         /* Initial printing */
    
    /* Tolerances */
    v[30] = atol;       /* Absolute function tolerance */
    v[31] = rtol;       /* Relative function tolerance */
    if ( xtol > -1.0f )
        v[32] = xtol;   /* X-tolerance */
    
    /* Ready to go, print some options */
    if ( printLevel > 1 )
    {
        mexPrintf("NL2SOL v2.13 - Settings\n");
        mexPrintf("  Maximum function evaluations: %d\n", maxFun );
        mexPrintf("  Maximum iterations:           %d\n", maxIter );
        mexPrintf("  Absolute function tolerance:  %e\n", v[30] );
        mexPrintf("  Relative function tolerance:  %e\n", v[31] );
        mexPrintf("  X tolerance:                  %e\n", v[32] );
        if ( useJacobian == 1 )
            mexPrintf("  Jacobian:                     User supplied\n" );
        else
            mexPrintf("  Jacobian:                     Finite differencing\n" );
        
        if ( bounded == 1 )
            mexPrintf("  Bounds:                       Enabled\n\n" );
        else
            mexPrintf("  Bounds:                       Not specified\n\n" );
    }        
    
    /* Run theg algorithm */
    if ( useJacobian )
    {
        if(bounded)
            dn2gb(&n,&p,x,bounds,CALCR,CALCJ,iv,&liv,&lv,v,uiparm,yData,&userData);
        else
            dn2g(&n,&p,x,CALCR,CALCJ,iv,&liv,&lv,v,uiparm,yData,&userData);
    } else {
        if(bounded)
            dn2fb(&n,&p,x,bounds,CALCR,iv,&liv,&lv,v,uiparm,yData,&userData);
        else
            dn2f(&n,&p,x,CALCR,iv,&liv,&lv,v,uiparm,yData,&userData);
    }
    
    /* Copy final values if succesful */
    getStatus(iv[0], printLevel);

    /* Handle optimization output */
    if ( nlhs > 1 ) {
        plhs[1]     = mxCreateDoubleMatrix(1,1, mxREAL);
        fVal        = mxGetPr( plhs[1] );
        *fVal       = 2*v[9];             /* output = 0.5 sum(r^2) */
    }
    
    if ( nlhs > 2 ) {
        plhs[2]     = mxCreateDoubleMatrix(1,1, mxREAL);
        exitFlag    = mxGetPr( plhs[2] );
        *exitFlag   = iv[0];
    }
    
    if ( nlhs > 3 ) {
        plhs[3]     = mxCreateDoubleMatrix(1,1, mxREAL);
        iter        = mxGetPr( plhs[3] );
        *iter       = (double)iv[30];
    }
    
    if ( nlhs > 4 ) {
        plhs[4]     = mxCreateDoubleMatrix(1,1, mxREAL);
        fEval       = mxGetPr( plhs[4] );
        *fEval      = (double)uiparm[0];
    }
    
    if ( nlhs > 5 ) {
        plhs[5]     = mxCreateDoubleMatrix(1,1, mxREAL);
        lambda      = mxGetPr( plhs[5] );
        *lambda     = 0;
    }
    
    /* Jacobian requested. Produce it again. */
    if ( nlhs > 6 ) {
        memcpy( mxGetPr(userData.callData[1]), x, (p)*sizeof(double) );
        userData.nrhs = 2;
        
        /* Call the MATLAB function one last time */
        result = mexCallMATLAB( 2, userData.plhs, userData.nrhs, userData.callData, "feval" );
        
        /* Point memory to the Jacobian */
        plhs[6]     = userData.plhs[1];
        
        /* Clean up memory for residual vector */
        mxDestroyArray(userData.plhs[0]);
    }
    
    if ( printLevel > 1 )
        mexPrintf( "Final SSE: %f. Function evaluations: %d. Jacobian evaluations: %d\n", 2*v[9], userData.nFunc, userData.nJac );

    /* Clean up our mess */
    mxDestroyArray( userData.callData[1] );

    if ( iv )
        mxFree( iv );
    
    if ( v )
        mxFree( v );
    
    if ( bounds )
        mxFree( bounds );
}

double getValueFromStruct( const mxArray *prhs[], const char* fieldName, double oldValue )
{
    double* value;
    int fieldNumber;
    mxArray *field;
    
    fieldNumber = mxGetFieldNumber(parOPTS, fieldName);
    
    /* Does the field exist? */
    if ( fieldNumber != -1 )
    {
        field = mxGetFieldByNumber(parOPTS, 0, fieldNumber);
        /* Not a valid value? Return -1.0f */
        if( ( field == NULL ) || !mxIsDouble(field) || mxIsComplex(field) || mxIsEmpty(field) )
        {
            if (!mxIsEmpty(field))
                mexPrintf( "Invalid value for field %s\n", fieldName );
            return oldValue;
        } else
        {
            /* Fetch the value */
            value = mxGetPr(field);
            return value[0];
        }
    } else
    {
        /* Not a valid value? Return -1.0f */
        return oldValue;
    }
}

/* Load options from options struct */
void loadOptions( const mxArray *prhs[], int* maxFun, int* maxIter, double* atol, double* rtol, double* xtol, int* printLevel, int* useJacobian )
{
    char* stringBuffer;
    int fieldNumber;
    mxArray *field;
    
    *maxIter    = (int) getValueFromStruct( prhs, "MaxIter", (double) *maxIter );
    *maxFun     = (int) getValueFromStruct( prhs, "MaxFunEvals", (double) *maxFun );
    *atol       = getValueFromStruct( prhs, "TolFun", (double) *atol );
    *rtol       = getValueFromStruct( prhs, "RelTolFun", (double) *rtol );
    *xtol       = getValueFromStruct( prhs, "TolX", (double) *xtol );
    
    /* String fields are done manually */
    fieldNumber = mxGetFieldNumber(parOPTS, "Jacobian");
    if ( fieldNumber != -1 )
    {
        field = mxGetFieldByNumber(parOPTS, 0, fieldNumber);
        
        /* Not a valid value? Return -1.0f */
        if( ( field == NULL ) || !mxIsChar(field) || mxIsEmpty(field) )
            *useJacobian = 0; /* Field is of invalid format */
        else
        {            
            if ( strcmp( mxArrayToString(field), "on" ) == 0 )
                *useJacobian = 1; /* Jacobian is on */
            else
                *useJacobian = 0; /* Jacobian is off */
        }
    } else
        *useJacobian = 0;  /* Field doesn't exist */
    
}

/* Parse input */
/* Adhere to same format as lsqnonlin:    X = fit(FUN,X0,LB,UB,OPTIONS) */
void validateInput( const mxArray *prhs[], int nrhs, int *npars, int *bounded )
{
    int pars;
    
    *bounded = 0;
    *npars = 0;
    
    if(nrhs < 2)
    {
        mexPrintf("NL2SOL v2.3 (by John Dennis, David Gay and Roy Welsch)\n");
        mexPrintf("MATLAB wrapper by Joep Vanlier (contact: joep.vanlier at g mail)\n");
        mexErrMsgTxt("You must supply at least 2 arguments to nl2sol!\n\nmexnl2sol(fun,x0,lb,ub,opts,printlevel)\n");
    }
        
    /* Check Types */
    if(!mxIsFunctionHandle(parFUN) && !mxIsChar(parFUN))
        mexErrMsgTxt("First argument must be a function handle that computes the RHS!");
    if(!mxIsDouble(parX0) || mxIsComplex(parX0) || mxIsEmpty(parX0))
        mexErrMsgTxt("Second argument must be a real vector of initial parameters");

    /* Get size of the parameter vector */
    pars = (int) mxGetNumberOfElements(parX0);
    *npars = pars;
    
    /* Check whether we have bounds */
    if ( nrhs == 3 )
        mexErrMsgTxt("Specify both lower and upper bounds or none at all.");
    
    if ( nrhs > 3 )
    {
        /* If both lb and ub are empty ==> no bounds */
        if( mxIsEmpty(parLB) && mxIsEmpty(parUB) )
            *bounded = 0;
        else {
            if (mxIsEmpty(parLB) || mxIsEmpty(parUB))
                mexErrMsgTxt("Specify both lower and upper bounds or none at all.");
 
            if( !mxIsDouble(parLB) || mxIsComplex(parLB) || ( pars != (int) mxGetNumberOfElements( parLB ) ) )
                mexErrMsgTxt("Lower bound must be specified as a numeric vector of the same length as p0");
    
            if( !mxIsDouble(parUB) || mxIsComplex(parUB) || ( pars != (int) mxGetNumberOfElements( parUB ) ) )
                mexErrMsgTxt("Upper bound must be specified as a numeric vector of the same length as p0");
            
            *bounded = 1;
        }
                
        if ( nrhs > 4 )
            if ( !mxIsStruct( parOPTS ) )
                mexErrMsgTxt("Fifth argument (when specified) must be an OPTIONS struct");
        
        if ( nrhs > 5 )
            if ( mxIsEmpty( parPRINTLEVEL ) || !mxIsDouble( parPRINTLEVEL ) || mxIsComplex( parPRINTLEVEL ) )
                mexErrMsgTxt("Last argument (when specified) must be a numeric value indicating a print level");        
    }
}

/* Pass residual to FORTRAN */
static void CALCR(int *n, int *p, double *x, int *nf, double *r, int *uiparm, double *ydata, void *ufparm)
{
    int result;                                 /* mexcall success */
    int i;                                      /* iterator */
    UserData *userData = (UserData *) ufparm;   /* function handle and call memory */
    double *returnVal;                          /* pointer to returned residual */
    double SSE;                                 /* SSE */
    userData->plhs[0] = NULL;
        
    /* Place parameters into function input */
    memcpy( mxGetPr(userData->callData[1]), x, (*p)*sizeof(double) );
    userData->nrhs = 2;             /* First is the function, second the parameters */
    
    /* Call MATLAB handle */
	result = mexCallMATLAB( 1, userData->plhs, userData->nrhs, userData->callData, "feval" );
    if ( result )
        mexErrMsgTxt("Error calling objective function");
    
    /* Fetch result */
    returnVal = mxGetPr(userData->plhs[0]);

    /* Assign residuals to FORTRAN memory */
    for ( i = 0; i < (*n); i++ )
    {
        r[i] = returnVal[i];
        
        if(mxIsInf(r[i]) || mxIsNaN(r[i]))
            *nf = 0;
    }
    
    /* Clean up result */
    mxDestroyArray(userData->plhs[0]);    
    
    /* Print iter */
    SSE = 0;
    if(uiparm[1] > 2) {   
        if(uiparm[0] == 1)
            mexPrintf("fEval     fJac          SSE\n");
        
        for( i = 0; i < (*n); i++ )
            SSE += r[i]*r[i];
        
        mexPrintf("%5d   %3d    %12.5g\n", userData->nFunc, userData->nJac, SSE);
        mexEvalString("drawnow;");
    }
    
    /* Advance iteration count */
    uiparm[0]++;
    userData->nFunc++;
}

/* Pass jacobian to FORTRAN */
static void CALCJ(int *n, int *p, double *x, int *nf, double *j, int *uiparm, double *ydata, void *ufparm)
{
    int result;                                 /* mexcall success */
    int i;                                      /* iterator */
    UserData *userData = (UserData *) ufparm;   /* function handle and call memory */
    double *returnVal;                          /* pointer to returned Jacobian */
    
    /* Place parameters into function input */
    memcpy( mxGetPr(userData->callData[1]), x, (*p)*sizeof(double) );
    userData->nrhs = 2;                         /* First is the function, second the parameters */
    
    /* Call MATLAB handle */
	result = mexCallMATLAB( 2, userData->plhs, userData->nrhs, userData->callData, "feval" );
    if ( result )
        mexErrMsgTxt("Error calling objective function");
    
    /* Fetch result */
    returnVal = mxGetPr(userData->plhs[1]);

    /* Assign Jacobian to FORTRAN memory */
    for ( i = 0; i < (*n)*(*p); i++ )
        j[i] = returnVal[i];
    
    /* Clean up result */
    mxDestroyArray(userData->plhs[0]);
    mxDestroyArray(userData->plhs[1]);
    
    /* Advance Jacobian count */
    userData->nJac++;
}

void printIt( const char* message, int printLevel )
{
    if ( printLevel > 0 )
        mexPrintf( "Exit criterion: %s", message );
}

double getStatus(int stat, int printLevel)
{
    int yes;
    
    switch((int)stat)
    {     
        case 3:
            printIt( "x tolerance reached\n", printLevel );
        case 4:
            printIt( "function tolerance reached\n", printLevel );
        case 5:
            printIt( "both tolerances reached\n", printLevel );
        case 6:
            printIt( "absolute function tolerance reached\n", printLevel );
            return 1;
            break;
        case 9:
            printIt( "function evaluation limit exceeded\n", printLevel );
        case 10:
            printIt( "iteration limit exceeded\n", printLevel );
            return 0;
            break;
        case 8:
            printIt( "tolerances too small\n", printLevel );
            return -1;
        case 13:
            printIt( "bad initial function evaluation\n", printLevel );
        case 14:
            printIt( "bad parameters\n", printLevel );
        case 15:
            printIt( "bad gradient\n", printLevel );
        case 16:
            printIt( "n/p out of range\n", printLevel );
        case 17:
        case 18:
            printIt( "iv out of range\n", printLevel );
        case 19:
            printIt( "v out of range\n", printLevel );
        case 50:
        case 87:
            printIt( "v problem\n", printLevel );
            return -2;
        case 11:
            return -5;
        default:
            return -3;        
    }
}