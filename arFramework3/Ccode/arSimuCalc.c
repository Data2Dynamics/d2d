/*
 *  MATLAB usage: arSimuCalc(struct ar, int fine, int sensi)
 *
 *  Copyright Andreas Raue 2011 (andreas.raue@fdm.uni-freiburg.de)
 *
 */

#define MACRO_DEBUGPRINT
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include "inverseC.h"
#ifndef MACRO_DEBUGPRINT
#include <stdarg.h>
#endif

#ifdef HAS_PTHREAD
#include <pthread.h>
#endif

#ifdef _WIN32
#include <winsock.h>
#else
#include <sys/time.h>
#endif

/* additional sys/time functions missing under Windows */
#ifdef _WIN32
#include <sys/timeb.h>
#include <sys/types.h>
#include <winsock.h>
void gettimeofday(struct timeval* t, void* timezone)
{
    struct _timeb timebuffer;
    _ftime( &timebuffer );
    t->tv_sec=timebuffer.time;
    t->tv_usec=1000*timebuffer.millitm;
}
void timersub(struct timeval* tvp, struct timeval* uvp, struct timeval* vvp)
{
    vvp->tv_sec = tvp->tv_sec - uvp->tv_sec;
    vvp->tv_usec = tvp->tv_usec - uvp->tv_usec;
    if (vvp->tv_usec < 0)
    {
        --vvp->tv_sec;
        vvp->tv_usec += 1000000;
    }
}
#endif

/* Prototypes for the debug printer */
#ifndef MACRO_DEBUGPRINT
    void debugPrint( int debugMode, int level, const char* format, ... );
#else
    /* Macro for the debug printer */
    #define DEBUGPRINT0(DBGMODE, LVL, STR) { if ( DBGMODE > LVL ) { mexPrintf( " [D] " ); mexPrintf( STR ); mexEvalString("drawnow;"); } }
    #define DEBUGPRINT1(DBGMODE, LVL, STR, ARG1) { if ( DBGMODE > LVL ) { mexPrintf( " [D] " ); mexPrintf( STR, ARG1 ); mexEvalString("drawnow;"); } }
    #define DEBUGPRINT2(DBGMODE, LVL, STR, ARG1, ARG2) { if ( DBGMODE > LVL ) { mexPrintf( " [D] " ); mexPrintf( STR, ARG1, ARG2 ); mexEvalString("drawnow;"); } }
    #define DEBUGPRINT3(DBGMODE, LVL, STR, ARG1, ARG2, ARG3) { if ( DBGMODE > LVL ) { mexPrintf( " [D] " ); mexPrintf( STR, ARG1, ARG2, ARG3 ); mexEvalString("drawnow;"); } }
    #define DEBUGPRINT4(DBGMODE, LVL, STR, ARG1, ARG2, ARG3, ARG4) { if ( DBGMODE > LVL ) { mexPrintf( " [D] " ); mexPrintf( STR, ARG1, ARG2, ARG3, ARG4 ); mexEvalString("drawnow;"); } }
#endif
            
#include <cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDENSE fcts. and constants */
#include <cvodes/cvodes_sparse.h>     /* prototype for CVSPARSE fcts. and constants */
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */
/*#include <cvodes/cvodes_superlumt.h> */  /* prototype for CVSUPERLUMT */
#include <sundials/sundials_sparse.h> /* definitions SlsMat */
#include <cvodes/cvodes_klu.h> /* definition of CVKLU sparse solver */

/* Accessor macros */
#define Ith(v, i)     NV_Ith_S(v, i-1)        /* i-th vector component i=1..neq */
#define IJth(A, i, j) DENSE_ELEM(A, i-1, j-1) /* (i,j)-th matrix component i,j=1..neq */

#define MXSTRING		32
#define MXNCF        20
#define MXNEF        20

#ifdef HAS_PTHREAD
struct thread_data_x {
    int	id;
};
#endif

int threadStatus[NMAXTHREADS];
int threadAbortSignal[NMAXTHREADS];

mxArray *armodel;
mxArray *arthread;

int    rootFinding;
int    done;
int    fine;
int    globalsensi;
int    dynamics;
int    ssa;
int    jacobian;
int    setSparse;
int    ms;
int    events;
int    parallel;
int    sensirhs;
int    debugMode;
int    sensitivitySubset;
/*int    fiterrors;*/
int    cvodes_maxsteps;
double cvodes_maxstepsize;
int    cvodes_atolV;
double cvodes_rtol;
double cvodes_atol;
/*double fiterrors_correction;*/
/* int useFitErrorMatrix;
 double *fiterrors_matrix;
 mwSize nrows_fiterrors_matrix; */

/* Name of the substructure with conditions we're currently evaluating */
char condition_name[MXSTRING];

/* Name of the threads substructure we are currently using */
char threads_name[MXSTRING];

/* Equilibration variables (Used when time point Inf is encountered) */
int    max_eq_steps;          /* Maximal equilibration steps */
double init_eq_step;          /* Initial equilibration stepsize attempt */
double eq_step_factor;        /* Factor with which to increase the time at each equilibration step */
double eq_tol;                /* Absolute tolerance for equilibration */
double eq_rtol;               /* Relative tolerance for equilibration */

struct timeval t1;

double  mintau;
int     nruns;

/* Prototype of undocumented MATLAB function */
#ifdef ALLOW_INTERRUPTS
extern bool utIsInterruptPending(void);
#endif

/* Prototypes of private functions */
#ifdef HAS_PTHREAD
void *thread_calc(void *threadarg);
#else
void thread_calc(int id);
#endif
void x_calc(int im, int ic, int sensi, int setSparse, int *threadStatus, int *abortSignal, int rootFinding, int debugMode, int sensitivitySubset);
void z_calc(int im, int ic, int isim, mxArray *arcondition, int sensi);
void y_calc(int im, int id, mxArray *ardata, mxArray *arcondition, int sensi);

void y_checkNaN(int nt, int ny, int it, double *y, double *yexp, double *ystd);
/* void fres(int nt, int ny, int it, double *res, double *y, double *yexp, double *ystd, double *chi2, double fiterrors_correction_factor);
 void fsres(int nt, int ny, int np, int it, double *sres, double *sy, double *yexp, double *ystd, double fiterrors_correction_factor);
 void fres_error(int nt, int ny, int it, double *reserr, double *res, double *y, double *yexp, double *ystd, double *chi2);
 void fsres_error(int nt, int ny, int np, int it, double *sres, double *sreserr, double *sy, double *systd, double *y, double *yexp, double *ystd, double *res, double *reserr); */

int ewt(N_Vector y, N_Vector w, void *user_data);
void thr_error( const char* msg );
int fetch_vector( mxArray* arcondition, int ic, double **vector, const char* fieldname, int desiredLength );
int init_list( mxArray* arcondition, int ic, double tstart, int* nPoints, double** timePoints, int* currentIndex, const char* flagFieldName, const char* timePointFieldName );
void copyStates( N_Vector x, double *returnx, double *qpositivex, int neq, int nout, int offset );
void copyResult( double* data, double *returnvec, int nu, int nout, int offset );
void copyNVMatrixToDouble( N_Vector* sx, double *returnsx, int nps, int neq, int nout, int offset );
void subCopyNVMatrixToDouble( N_Vector* sx, double *returnsx, int nps, int neq, int nout, int offset, int32_T* targetIdx );

/* user functions */
#include "arSimuCalcFunctions.c"

void storeSimulation( UserData data, int im, int isim, int is, int nu, int nv, int neq, int nout, N_Vector x, double *returnx, double *returnu, double *returnv, double *qpositivex );
void storeSensitivities( UserData data, int im, int isim, int is, int np, int nu, int nv, int neq, int nout, N_Vector x, N_Vector *sx, double *returnsx, double *returnsu, double *returnsv, int sensitivitySubset, int32_T *sensitivityMapping );
void findRoots( SimMemory sim_mem, mxArray *arcondition, int im, int ic, int isim, double tstart, double eq_tol, int neq, int nu, int nv, int nout, double* returnx, double* returnu, double* returnv, double* qpositivex, double* returnsx, double* returnsu, double* returnsv, int sensi, int ysensi, int npSensi, int has_tExp );
void storeIntegrationInfo( SimMemory sim_mem, mxArray *arcondition, int ic );
void terminate_x_calc( SimMemory sim_mem, double status );
void initializeDataCVODES( SimMemory sim_mem, double tstart, int *abortSignal, mxArray *arcondition, double *qpositivex, int ic, int nsplines, int sensitivitySubset );
int allocateSimMemoryCVODES( SimMemory sim_mem, int neq, int np, int sensi, int npSensi );
int allocateSimMemorySSA( SimMemory sim_mem, int nx );
int applyInitialConditionsODE( SimMemory sim_mem, double tstart, int im, int isim, double *returndxdt, double *returndfdp0, mxArray *x0_override, int sensitivitySubset );
int initializeEvents( SimMemory sim_mem, mxArray *arcondition, int ic, double tstart );
void evaluateObservations( mxArray *arcondition, int im, int ic, int sensi, int has_tExp );

int handle_event( SimMemory sim_mem, int sensi_meth, int reinitSolver );
int equilibrate(void *cvode_mem, UserData user_data, N_Vector x, realtype t, double *equilibrated, double *returndxdt, double *teq, int neq, int im, int ic, int *abortSignal );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int nthreads, ithreads, tid;

#ifdef HAS_PTHREAD
    int rc;
    pthread_t threads_x[NMAXTHREADS];
    struct thread_data_x thread_data_x_array[NMAXTHREADS];
#endif

    mxArray    *arconfig;
    struct timeval t2, tdiff;
    double *ticks_stop;
    ticks_stop = mxGetData(mxGetField(prhs[0], 0, "stop"));
    gettimeofday(&t1, NULL);
    
    srand( (unsigned)time(NULL) );
    
    /* get ar.model */
    armodel = mxGetField(prhs[0], 0, "model");
    if(armodel==NULL){
        mexErrMsgTxt("field ar.model doesn't exist");
    }
    
    fine = (int) mxGetScalar(prhs[1]);
    globalsensi = (int) mxGetScalar(prhs[2]);
    dynamics = (int) mxGetScalar(prhs[3]);
    ssa = (int) mxGetScalar(prhs[4]);

    if ( mxGetString(prhs[5], condition_name, MXSTRING ) != 0 )
        mexErrMsgTxt("Failed to provide condition identifier to arSimuCalc. Aborting ...");
    
    if ( mxGetString(prhs[6], threads_name, MXSTRING ) != 0 )
        mexErrMsgTxt("Failed to provide name of condition list to arSimuCalc. Aborting ...");    

    /* Is the equilibrium found by rootfinding and must we immediately terminate? */
    if ( nrhs > 7 ) {
        rootFinding = (int) mxGetScalar(prhs[7]);
    } else {
        rootFinding = 0;
    }
    
    /* get ar.config */
    arconfig = mxGetField(prhs[0], 0, "config");
    parallel = (int) mxGetScalar(mxGetField(arconfig, 0, "useParallel"));
    jacobian = (int) mxGetScalar(mxGetField(arconfig, 0, "useJacobian"));
    setSparse = (int) mxGetScalar(mxGetField(arconfig, 0, "useSparseJac"));
    sensirhs = (int) mxGetScalar(mxGetField(arconfig, 0, "useSensiRHS"));
    cvodes_atolV = (int) mxGetScalar(mxGetField(arconfig, 0, "atolV"));
    cvodes_rtol = mxGetScalar(mxGetField(arconfig, 0, "rtol"));
    cvodes_atol = mxGetScalar(mxGetField(arconfig, 0, "atol"));
    cvodes_maxsteps = (int) mxGetScalar(mxGetField(arconfig, 0, "maxsteps"));
    cvodes_maxstepsize = mxGetScalar(mxGetField(arconfig, 0, "maxstepsize"));
    
    /* Do we want debug mode? */
    debugMode = 0;
    if ( mxGetField(arconfig, 0, "debug" ) )
        debugMode = (int) mxGetScalar(mxGetField(arconfig, 0, "debug"));
    
    sensitivitySubset = 0;
    if ( mxGetField(arconfig, 0, "sensitivitySubset" ) )
        sensitivitySubset = (int) mxGetScalar(mxGetField(arconfig, 0, "sensitivitySubset"));
    
    /* In debug mode we have to disable threading, since otherwise the mexPrintf can lead to a race condition which may crash MATLAB */
    if ( debugMode > 0 )
    {
        mexPrintf("[DEBUG MODE: Parallelization disabled]\n");
        parallel = 0;
    }
    
    mintau = mxGetScalar(mxGetField(arconfig, 0, "ssa_min_tau"));
    nruns = (int) mxGetScalar(mxGetField(arconfig, 0, "ssa_runs"));
    ms = (int) mxGetScalar(mxGetField(arconfig, 0, "useMS"));
    events = (int) mxGetScalar(mxGetField(arconfig, 0, "useEvents"));
    max_eq_steps = (int) mxGetScalar(mxGetField(arconfig, 0, "max_eq_steps"));
    init_eq_step = (double) mxGetScalar(mxGetField(arconfig, 0, "init_eq_step"));
    eq_step_factor = (double) mxGetScalar(mxGetField(arconfig, 0, "eq_step_factor"));
    eq_tol = (double) mxGetScalar(mxGetField(arconfig, 0, "eq_tol"));
    if ( mxGetField(arconfig, 0, "eq_rtol") )
    {
        eq_rtol = (double) mxGetScalar(mxGetField(arconfig, 0, "eq_rtol"));
    } else {
        eq_rtol = 0.0;
    }
        

    DEBUGPRINT0( debugMode, 2, "Loaded configuration\n" );
    
    if ( ms == 1 ) events = 1;
        
    /* threads */
    arthread = mxGetField(arconfig, 0, threads_name);
    nthreads = (int) mxGetNumberOfElements(arthread);
    
#ifdef HAS_PTHREAD
    if(NMAXTHREADS<nthreads) mexErrMsgTxt("ERROR at NMAXTHREADS < nthreads");
#endif
        
#ifdef HAS_PTHREAD
    /* loop over threads parallel */
    for(ithreads=0; ithreads<nthreads; ++ithreads){
        threadStatus[ithreads] = 0;
        threadAbortSignal[ithreads] = 0;
        tid = (int) mxGetScalar(mxGetField(arthread, ithreads, "id"));
        
        thread_data_x_array[tid].id = tid;

        if(parallel==1){
            /* printf("creating thread %i\n", tid); */
            rc = pthread_create(&threads_x[tid], NULL, thread_calc, (void *) &thread_data_x_array[tid]);
            if (rc){
                mexErrMsgTxt("ERROR at pthread_create");
            }
        } else {
            thread_calc(&thread_data_x_array[tid]);
        }
    }
    
    /* wait for termination of condition threads, but make sure program is interruptible */
    #ifdef ALLOW_INTERRUPTS
    if (parallel==1) {
        done = 0;
        while(done == 0) {
            done = 1;
            for(ithreads=0; ithreads<nthreads; ++ithreads){
                if ( threadStatus[ithreads] == 0 ) {
                    done = 0;
                }
            }
            if (utIsInterruptPending()) {
                for(ithreads=0; ithreads<nthreads; ++ithreads){
                    threadAbortSignal[ithreads] = 1;
                }
                mexPrintf( "Interrupt detected => Aborting simulation\n" );
                done = 1;
            }
        }
    }
    #endif
    
    /* join condition threads */
    if(parallel==1){
        for(ithreads=0; ithreads<nthreads; ++ithreads){
            tid = (int) mxGetScalar(mxGetField(arthread, ithreads, "id"));
            
            rc = pthread_join(threads_x[tid], NULL);
            if (rc){
                mexErrMsgTxt("ERROR at pthread_join");
            }
        }
    }
#else
    /* loop over threads sequential */
    for(ithreads=0; ithreads<nthreads; ++ithreads){
        tid = (int) mxGetScalar(mxGetField(arthread, ithreads, "id"));
        thread_calc(tid);
    }
#endif    

    gettimeofday(&t2, NULL);
    timersub(&t2, &t1, &tdiff);
    ticks_stop[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
}

/* work of threads */
#ifdef HAS_PTHREAD
void *thread_calc(void *threadarg) {
    struct thread_data_x *my_data = (struct thread_data_x *) threadarg;
    int id = my_data->id;
#else
void thread_calc(int id) {
#endif
    
    int n = (int) mxGetScalar(mxGetField(arthread, id, "n"));
    int *ms = (int *) mxGetData(mxGetField(arthread, id, "ms"));
    int *cs = (int *) mxGetData(mxGetField(arthread, id, "cs"));
    int in;
    
    /* printf("computing thread #%i\n", id); */
    DEBUGPRINT0( debugMode, 2, "Calling conditions\n" );
    for(in=0; in<n; ++in){
        /* printf("computing thread #%i, task %i/%i (m=%i, c=%i)\n", id, in, n, ms[in], cs[in]); */
        x_calc(ms[in], cs[in], globalsensi, setSparse, &threadStatus[id], &threadAbortSignal[id], rootFinding, debugMode, sensitivitySubset);
    }
    /* printf("computing thread #%i(done)\n", id); */
    
#ifdef HAS_PTHREAD
    if(parallel==1) {pthread_exit(NULL);}
    return NULL;
#endif
}

/* Handle CVODES errors */
/* CAUTION: this function is NOT thread safe! */
void errorHandler(int error_code, const char *module, const char *func, char *msg, void *eh_data)
{
	mexPrintf( "Error code %d in module %s and function %s:\n%s\n", error_code, module, func, msg );
};

/* Function which can be used for debugging purposes */
/* CAUTION: this function is NOT thread safe! */
void waitForKey()
{
    mxArray   *junk, *str;
    str = mxCreateString("Hit [ENTER] to continue . . .");
    mexCallMATLAB(1,&junk,1,&str,"input");
};

/* calculate dynamics */
void x_calc(int im, int ic, int sensi, int setSparse, int *threadStatus, int *abortSignal, int rootFinding, int debugMode, int sensitivitySubset) {
    mxArray    *x0_override;
    mxArray    *arcondition;
    
    int nsplines;
    int ysensi;
    int has_tExp;
    int nm, nc;
    int flag;
    int is, js, ks, jss;
    int nout; /*, nyout;*/
    int nu, nv, nnz;
    
    /* Which condition to simulate */
    int isim;
       
    /* Used to override which condition to simulate */
    mxArray *src;    
    double  *isrc;
    int     only_sim;

    /* Multiple shooting and events */
    int qMS, qEvents;
    EventData event_data;

    void *cvode_mem;
    UserData data;
    
    realtype t;
    double tstart;
    double inf;
    realtype *atolV_tmp;
    realtype *sxtmp;
    
    /* SSA variables */
    double tfin, tau, meantau;
    double r1, r2;
    double alpha0, sumalpha;
    int iruns, it, itexp, ix, iv;
    double lasttau[] = {1,1,1,1,1,1,1,1,1,1};
    int ilasttau = 0;
    double *texp;
    int ntexp;
    double *xssaexp;
    
    double *qpositivex;
    double *status;
    double *ts;
    double *teq;
    double *returnu;
    double *returnsu;
    double *returnv;
    double *returnsv;
    double *returnx;
    double *returnsx;
    double *returndxdt;
    double *returndfdp0;
    double *equilibrated;
    double *y_max_scale;
    
    struct timeval t2;
    struct timeval t3;
    struct timeval t4;
    struct timeval tdiff;
    double *ticks_start;
    double *ticks_stop_data;
    double *ticks_stop;
       
    /* List of indices which map the sensitivities back to the output ones */
    int32_T *sensitivityMapping;
       
    /* Pointer to centralized container for the heap memory */
    SimMemory sim_mem = NULL;
    int np, neq, npSensi;
    
    /* Pointers to heap memory (which needs to be cleaned up!) */
    /* CVODES */
    N_Vector x = NULL;
    N_Vector atolV = NULL;
    N_Vector atols_ss = NULL;
    N_Vector *atolV_ss = NULL;
    N_Vector *sx = NULL;
    
    /* SSA */
    N_Vector x_lb = NULL;
    N_Vector x_ub = NULL;
    
    int sensi_meth = CV_SIMULTANEOUS; /* CV_SIMULTANEOUS or CV_STAGGERED */
    bool error_corr = TRUE;
    only_sim = 0;
    
    DEBUGPRINT0( debugMode, 4, "Entry point x_calc\n" );
    
    /* Grab value of infinity (used to mark steady state simulations) */
    inf = mxGetInf();

    /* check if im in range */
    nm = (int) mxGetNumberOfElements(armodel);
    if(nm<=im) {
        thr_error("im > length(ar.model)\n");
        *threadStatus = 1;
        return;
    }
    
    /* get ar.model(im).condition */
    arcondition = mxGetField(armodel, im, condition_name);
    if(arcondition==NULL){ *threadStatus = 1; return; }
    
    /* check if ic in range */
    nc = (int) mxGetNumberOfElements(arcondition);
    if(nc<=ic) { thr_error("ic > length(ar.model.condition)\n"); *threadStatus = 1; return; }
    
    /* Initialize memory to facilitate easier cleanup */
    status = mxGetData(mxGetField(arcondition, ic, "status"));
    sim_mem = simCreate( threadStatus, status );
    
    /* Get double handle to store equilibrium value */
    teq = mxGetData(mxGetField(arcondition, ic, "tEq"));
    
    has_tExp = (int) mxGetScalar(mxGetField(arcondition, ic, "has_tExp"));
    if(has_tExp == 0 && fine == 0) { terminate_x_calc( sim_mem, 0 ); DEBUGPRINT0( debugMode, 4, "Terminated simulation since tExp simulation was requested (fine = 0) and there are no experimental points.\n" ); return; }

    ticks_start = mxGetData(mxGetField(arcondition, ic, "start"));
    ticks_stop = mxGetData(mxGetField(arcondition, ic, "stop"));
    ticks_stop_data = mxGetData(mxGetField(arcondition, ic, "stop_data"));

    gettimeofday(&t2, NULL);
    
    DEBUGPRINT0( debugMode, 4, "Starting dynamic section\n" );
    
    if(dynamics == 1) {
        if(ssa == 0) {
            /**** begin of CVODES ****/
            /* NOTE: ic refers to the target condition in the target ar.model(#).condition(#) structure         */
            /* isim refers to the condition that is actually simulated. In most cases these will be the same    */
            /* but for steady state simulation, isim is redirected using the ar.model(#).condition(#).src field */
            
            /* get MATLAB values */
            qpositivex = mxGetData(mxGetField(armodel, im, "qPositiveX"));
            tstart = mxGetScalar(mxGetField(arcondition, ic, "tstart"));
            neq = (int) mxGetNumberOfElements(mxGetField(armodel, im, "xs"));
            nnz = (int) mxGetScalar(mxGetField(armodel, im, "nnz"));
     
            if(fine == 1){
                DEBUGPRINT0( debugMode, 4, "Performing fine simulation\n" );
                ts = mxGetData(mxGetField(arcondition, ic, "tFine"));
                nout = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "tFine"));
                
                returnu = mxGetData(mxGetField(arcondition, ic, "uFineSimu"));
                returnv = mxGetData(mxGetField(arcondition, ic, "vFineSimu"));
                returnx = mxGetData(mxGetField(arcondition, ic, "xFineSimu"));
                y_max_scale = mxGetData(mxGetField(arcondition, ic, "y_atol"));
                if (sensi == 1) {
                    returnsu = mxGetData(mxGetField(arcondition, ic, "suFineSimu"));
                    returnsv = mxGetData(mxGetField(arcondition, ic, "svFineSimu"));
                    returnsx = mxGetData(mxGetField(arcondition, ic, "sxFineSimu"));
                }
            }
            else{
                DEBUGPRINT0( debugMode, 4, "Performing experiment simulation\n" );
                ts = mxGetData(mxGetField(arcondition, ic, "tExp"));
                nout = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "tExp"));
                
                returnu = mxGetData(mxGetField(arcondition, ic, "uExpSimu"));
                returnv = mxGetData(mxGetField(arcondition, ic, "vExpSimu"));
                returnx = mxGetData(mxGetField(arcondition, ic, "xExpSimu"));
                
                y_max_scale = mxGetData(mxGetField(arcondition, ic, "y_atol"));
                
                if (sensi == 1) {
                    returnsu = mxGetData(mxGetField(arcondition, ic, "suExpSimu"));
                    returnsv = mxGetData(mxGetField(arcondition, ic, "svExpSimu"));
                    returnsx = mxGetData(mxGetField(arcondition, ic, "sxExpSimu"));
                }
            }
            
            returndxdt = mxGetData(mxGetField(arcondition, ic, "dxdt"));
            if (sensi == 1) {
                returndfdp0 = mxGetData(mxGetField(arcondition, ic, "ddxdtdp"));
                DEBUGPRINT0( debugMode, 4, "Sensitivities enabled\n" );
                
                /* Obtain indices which map the computed sensitivities back to the output ones */
                if ( sensitivitySubset == 1 )
                {
                    DEBUGPRINT0( debugMode, 4, "Using sensitivity subset\n" );
                    sensitivityMapping = (int32_T *) mxGetData(mxGetField(arcondition, ic, "backwardIndices"));
                }                    
            }

            /* Fetch number of inputs, parameters and fluxes */
            nu = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "uNum"));
            np = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "pNum"));
            nv = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "vNum"));
            if ( mxGetField(arcondition, ic, "splines") )
                nsplines = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "splines"));
            else
                nsplines = 0;      
            
            if ( sensitivitySubset == 1 )
                npSensi = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "sensIndices"));
            else
                npSensi = np;
            
            /* If there are no parameters, do not compute simulated sensitivities; otherwise failure at N_VCloneVectorArray_Serial */
            ysensi = sensi;
            if (npSensi==0) sensi = 0;
            
            /* Allocate heap memory required for simulation */
            DEBUGPRINT0( debugMode, 4, "Attempting to allocate CVODES memory\n" );
            if ( allocateSimMemoryCVODES( sim_mem, neq, np, sensi, npSensi ) )
            {
                /* Generate some local references to avoid having sim_mem-> littered everywhere */
                x = sim_mem->x;
                sx = sim_mem->sx;
                atolV = sim_mem->atolV;
                atols_ss = sim_mem->atols_ss;
                atolV_ss = sim_mem->atolV_ss;
                data = sim_mem->data;
                event_data = sim_mem->event_data;
                cvode_mem = sim_mem->cvode_mem;
            } else return;
            
            DEBUGPRINT3( debugMode, 4, "np: %d npSensi: %d neq: %d\n", np, npSensi, neq );
            
            /* User data structure */
            DEBUGPRINT0( debugMode, 4, "Initialize CVODES data\n" );
            initializeDataCVODES( sim_mem, tstart, abortSignal, arcondition, qpositivex, ic, nsplines, sensitivitySubset );
            
            /* Initialize event system */
            qEvents = 0;
            if ( events )
            {
                qEvents = initializeEvents( sim_mem, arcondition, ic, tstart );
                DEBUGPRINT0( debugMode, 4, "Events initialized\n" );
            }
            
            /* Initialize multiple shooting list */
            if (ms==1) 
                qMS = init_list(arcondition, ic, tstart, &(event_data->nMS), &(event_data->tMS), &(event_data->iMS), "qMS", "tMS");
            
            /* Override which condition to simulate. For normal conditions isim is equal to ic (both are condition indices) */
            /* For steady state simulations, isim is equal to another condition, since the index which indicates the location in ar.ss_conditions differs from the condition referenced in ar.conditions */
            src = mxGetField(arcondition, ic, "src");
            if (src == NULL) {
                only_sim    = 0;
                isim        = ic;
            } else {
                isrc        = mxGetData(src);
                isim        = (int) (*isrc) - 1;
                only_sim    = 1;
            }
            
            /* Is there a list which states to equilibrate? */
            if ( mxGetField( arcondition, ic, "ssStates" ) ) {
                equilibrated = mxGetData(mxGetField(arcondition, ic, "ssStates"));
            } else {
                equilibrated = NULL;
            }

            /* Apply ODE initial conditions */
            x0_override = mxGetField(arcondition, ic, "x0_override");
            if ( !applyInitialConditionsODE( sim_mem, tstart, im, isim, returndxdt, returndfdp0, x0_override, sensitivitySubset ) )
                return;

            /* Do we have a startup event? */
            if ( qEvents == 1 ) {
                DEBUGPRINT0( debugMode, 4, "Handling startup event\n" );
                if ( event_data->t[event_data->i] == tstart ) {
                    flag = handle_event( sim_mem, sensi_meth, 0 );
                    (event_data->i)++;

                    if (flag < 0) {thr_error("Failed to reinitialize solver at event"); terminate_x_calc( sim_mem, 16 ); return;}
                }
            }            
            
            /* Check if we are only simulating dxdt */
            if ( rootFinding > 0 )
            {
                findRoots( sim_mem, arcondition, im, ic, isim, tstart, eq_tol, neq, nu, nv, nout, returnx, returnu, returnv, qpositivex, returnsx, returnsu, returnsv, sensi, ysensi, npSensi, has_tExp );
                return;
            }
            
            if(neq>0){
                /* Allocate space for CVODES */
                DEBUGPRINT0( debugMode, 4, "Allocating memory for CVODES\n" );
                flag = AR_CVodeInit(cvode_mem, x, tstart, im, isim);
                if (flag < 0) {terminate_x_calc( sim_mem, 4 ); return;}
                
                DEBUGPRINT0( debugMode, 4, "Setting CVODES options\n" );
                /* Optionally enable more informative debug messages */
                if ( debugMode > 0 )
                {
                    flag = CVodeSetErrHandlerFn(cvode_mem, &errorHandler, NULL);
                    if (flag < 0) {terminate_x_calc( sim_mem, 4 ); return;}
                }
                
                /* Number of maximal internal steps */
                flag = CVodeSetMaxNumSteps(cvode_mem, cvodes_maxsteps);
                if(flag < 0) {terminate_x_calc( sim_mem, 15 ); return;}
                
                /* Maximal internal step size */
                flag = CVodeSetMaxStep(cvode_mem, cvodes_maxstepsize);
                if(flag < 0) {terminate_x_calc( sim_mem, 19 ); return;}
                	 
                if(cvodes_atolV==1) { 	
                    double tmp_tol = 1.;
                    DEBUGPRINT0( debugMode, 4, "Setting atolV\n" );
                    for(ks=0; ks < neq; ks++) {		    
                        if(y_max_scale[ks]==0 || cvodes_atol/y_max_scale[ks]>1){
                            Ith(atolV, ks+1) = 1;
                        }else if(cvodes_atol/y_max_scale[ks]<1e-8){
                            Ith(atolV, ks+1) = 1e-8;			
                            /*printf("atolV for neq=%i is %d \n", ks+1, Ith(atolV, ks+1));*/
                        }else if(cvodes_atol/y_max_scale[ks]>1e-8 && cvodes_atol/y_max_scale[ks]<1){
                            Ith(atolV, ks+1) = cvodes_atol/y_max_scale[ks];
                        }else{
                            Ith(atolV, ks+1) = cvodes_atol;
                        }
                        /*printf("atolV for neq=%i is %d \n", ks, Ith(atolV, ks+1));*/
                        tmp_tol *= Ith(atolV, ks+1);
                    }
                    tmp_tol = cvodes_atol / pow(tmp_tol,1/neq);  
                    for(ks=0; ks < neq; ks++) {	
                        Ith(atolV, ks+1) = Ith(atolV, ks+1) * tmp_tol;
                    }
                    flag = CVodeSVtolerances(cvode_mem, RCONST(cvodes_rtol), atolV);
                } else {                
                    flag = CVodeSStolerances(cvode_mem, RCONST(cvodes_rtol), RCONST(cvodes_atol));
                }
                if (flag < 0) {terminate_x_calc( sim_mem, 5 ); return;}
                
                /* Attach user data */
                DEBUGPRINT0( debugMode, 4, "Attaching Userdata\n" );
                flag = CVodeSetUserData(cvode_mem, data);
                if (flag < 0) {terminate_x_calc( sim_mem, 6 ); return;}
                
                /* Attach linear solver */
                if(setSparse == 0){
                    /* Dense solver */
                    flag = CVDense(cvode_mem, neq);
                }else{              
                    /* sparse linear solver KLU */
                    flag = CVKLU(cvode_mem, neq, nnz);
                }
                if (flag < 0) {terminate_x_calc( sim_mem, 7 ); return;}
                
                /* Jacobian-related settings */
                if (jacobian == 1) {
                    flag = AR_CVDlsSetDenseJacFn(cvode_mem, im, isim, setSparse);
                    if (flag < 0) {terminate_x_calc( sim_mem, 8 ); return;}
                }
                
                /* custom error weight function */
                /*
                    flag = CVodeWFtolerances(cvode_mem, ewt);
                    if (flag < 0) return;
                 */
            }
            
            /********************************/
            /* Sensitivity-related settings */
            /********************************/
            if (sensi == 1) {
                DEBUGPRINT0( debugMode, 4, "Initializing sensitivities\n" );
                if(neq>0){
                    flag = AR_CVodeSensInit1(cvode_mem, npSensi, sensi_meth, sensirhs, sx, im, isim, sensitivitySubset);
                     
                    if(flag < 0) {terminate_x_calc( sim_mem, 10 ); return;}
                    
                    /*
                        flag = CVodeSensEEtolerances(cvode_mem);
                        if(flag < 0) {terminate_x_calc( sim_mem, 11 ); return;}
                     */
                    
                    if ( sensitivitySubset == 1 )
                        flag = CVodeSetSensParams(cvode_mem, NULL, NULL, NULL);                     /* Parameters only need to be specified when sensis are computed */
                    else
                        flag = CVodeSetSensParams(cvode_mem, data->p, NULL, NULL);
                    
                    if (flag < 0) {terminate_x_calc( sim_mem, 13 ); return;}
                    
                    /* Set error weights */
                    for (is=0; is<npSensi; is++) Ith(atols_ss, is+1) = cvodes_atol;
                    
                    if(cvodes_atolV==1)
                    { 
                        for(js=0; js < npSensi; js++)
                        {
                            atolV_tmp = NV_DATA_S(atolV_ss[js]);
                            for(ks=0; ks < neq; ks++)
                            {
                                if(y_max_scale[ks]==0. || cvodes_atol/y_max_scale[ks]>1) {
                                    atolV_tmp[ks] = 1;
                                } else if (cvodes_atol/y_max_scale[ks]<1e-8) {   
                                    /* && Ith(atolV, ks+1)==1.e-8){*/
                                    /*printf("atolVS for neq=%i is %d \n", ks+1, atolV_tmp[ks]);*/ /*}*/
                                    atolV_tmp[ks] = 1e-8;			  
                                }else if(cvodes_atol/y_max_scale[ks]>1e-8 && cvodes_atol/y_max_scale[ks]<1) {
                                    atolV_tmp[ks] = cvodes_atol/y_max_scale[ks];
                                    /*if(atolV_tmp[ks] < Ith(atolV, ks+1)){
                                         atolV_tmp[ks] = Ith(atolV, ks+1);
                                    }*/
                                }else {
                                    atolV_tmp[ks] = cvodes_atol;
                                }			  
                                /*printf("atolV_ss for neq=%i is %f\n", ks, atolV_tmp[ks]);*/
                            }
                        }
                        flag = CVodeSensSVtolerances(cvode_mem, RCONST(cvodes_rtol), atolV_ss);
        		    }else
                    {
                      flag = CVodeSensSStolerances(cvode_mem, RCONST(cvodes_rtol), N_VGetArrayPointer(atols_ss));
                    }
                    
                    if(flag < 0) {terminate_x_calc( sim_mem, 11 ); return;}
                    
                    flag = CVodeSetSensErrCon(cvode_mem, error_corr);
                    if(flag < 0) {terminate_x_calc( sim_mem, 13 ); return;}
                }
            }

            /********************************/
            /* loop over output points      */
            /********************************/
            DEBUGPRINT0( debugMode, 4, "Starting main calculation loop\n" );
            for (is=0; is < nout; is++) {               
                /*mexPrintf("%f x-loop (im=%i ic=%i)\n", ts[is], im, ic);*/
                DEBUGPRINT3( debugMode, 7, "%f x-loop (im=%i ic=%i)\n", ts[is], im, ic );
                /* only integrate if no errors occured */
                if(status[0] == 0.0) {

                    /* only integrate after tstart */
                    if(ts[is] > tstart) {

                        if(neq>0) {
                            /* If this condition has events, make sure we don't go over them as this leads to loss of accuracy */
                            if (qEvents==1){
                                if (event_data->i < event_data->n)
                                    CVodeSetStopTime(cvode_mem, RCONST(event_data->t[event_data->i]));
                                else
                                    CVodeSetStopTime(cvode_mem, ts[nout-1]+1.0);
                            }
                            
                            if ( ts[is] == inf ) {
                                /* Equilibrate the system */
                                DEBUGPRINT1( debugMode, 4, "Equilibrating the system (t=%g)...\n", data->t );
                                flag = equilibrate(cvode_mem, data, x, t, equilibrated, returndxdt, teq, neq, im, isim, abortSignal);
                                DEBUGPRINT1( debugMode, 4, "[ OK ] (teq=%g)\n", teq[0] );
                                CVodeGetCurrentTime( cvode_mem, &(data->t) );
                                
                                if (flag < 0) {thr_error("Failed to equilibrate system"); terminate_x_calc( sim_mem, 20 ); return;}                                
                            } else {
                                /* Simulate up to the next time point */
                                flag = CVode(cvode_mem, RCONST(ts[is]), x, &t, CV_NORMAL);
                                data->t = ts[is];
                            }
                            
                            /* Found an event */
                            if ((qEvents==1) && (event_data->i < event_data->n) && (ts[is]==event_data->t[event_data->i])) /*flag==CV_TSTOP_RETURN*/
                            {
                              DEBUGPRINT0( debugMode, 5, "Handling event\n" );
                              qEvents = 2;    /* qEvents=2 denotes that an event just happened */
                              flag = 0;       /* Re-set the flag for legacy error-checking reasons */
                            }
                            
                            if ( flag==CV_TSTOP_RETURN )
                            {
                              thr_error( "Error in the event system. Did the model link properly?" );
                            }

                            status[0] = flag;
                        } else data->t = ts[is];
                    }
                }
                                
                /* Store time step results */
                DEBUGPRINT1( debugMode, 11, "Status: %g\n", status[0] );
                if(status[0] == 0.0) {
                    DEBUGPRINT0( debugMode, 11, "Storing results\n" );
                    storeSimulation( data, im, isim, is, nu, nv, neq, nout, x, returnx, returnu, returnv, qpositivex );
                } else {
                    for(js=0; js < nu; js++) returnu[js*nout+is] = 0.0;
                    for(js=0; js < nv; js++) returnv[js*nout+is] = 0.0;
                    for(js=0; js < neq; js++) returnx[js*nout+is] = 0.0;
                }
                
                /* Store output sensitivities */
                if(status[0] == 0.0) {
                    if (sensi == 1) {
                        DEBUGPRINT0( debugMode, 6, "Storing sensitivities\n" );
                        if(ts[is] > tstart) {
                            if(neq>0) {
                                flag = CVodeGetSens(cvode_mem, &t, sx);
                                if (flag < 0) {terminate_x_calc( sim_mem, 14 ); return;}
                            }
                        }
                        storeSensitivities( data, im, isim, is, np, nu, nv, neq, nout, x, sx, returnsx, returnsu, returnsv, sensitivitySubset, sensitivityMapping );
                    }
                } else {
                    /* Store empty output sensitivities in case of an error */
                    if (sensi == 1) {
                        DEBUGPRINT0( debugMode, 6, "Storing zeroed sensitivities\n" );
                        for(js=0; js < np; js++) {
                            if(neq>0) {
                                sxtmp = NV_DATA_S(sx[js]);
                                for(ks=0; ks < neq; ks++) {
                                    returnsx[js*neq*nout + ks*nout + is] = 0.0;
                                }
                            }
                            for(ks=0; ks < nu; ks++) {
                                returnsu[js*nu*nout + ks*nout + is] = 0.0;
                            }
                            for(ks=0; ks < nv; ks++) {
                                returnsv[(js*nv+ks)*nout + is] = 0.0;
                            }                            
                        }
                    }
                }
                
                /* Event handling */
                if (qEvents==2)
                {
                    DEBUGPRINT0( debugMode, 6, "Handling event + reinitializing the solver\n" );
                    flag = handle_event( sim_mem, sensi_meth, 1 );
                    if (flag < 0) {thr_error("Failed to reinitialize solver at event"); terminate_x_calc( sim_mem, 16 ); return;}
                    
                    storeSimulation( data, im, isim, is, nu, nv, neq, nout, x, returnx, returnu, returnv, qpositivex );
                    if ( sensi == 1 )
                        storeSensitivities( data, im, isim, is, np, nu, nv, neq, nout, x, sx, returnsx, returnsu, returnsv, sensitivitySubset, sensitivityMapping );
                    
                    qEvents = 1;
                    (event_data->i)++;
                }
                
              DEBUGPRINT0( debugMode, 6, "Going into next iteration cycle\n" );
            } /* End of simulation loop */
            
            /* Store dfxdx */
            {
                double *dfdx, *dfdp;
                if( cvode_mem != NULL ) CVodeGetCurrentTime( cvode_mem, &t );
                DEBUGPRINT1( debugMode, 6, "Storing final dfdx and dfdp at %g\n", t );
                dfdx = mxGetData(mxGetField(arcondition, ic, "dfdxNum"));
                dfdp = mxGetData(mxGetField(arcondition, ic, "dfdpNum"));
                
                fsv(data, t, x, im, isim);                             /* Updates dvdp, dvdu, dvdx */
                getdfxdx(im, isim, t, x, dfdx, data);                  /* Updates dvdx and stores dfxdx */
                dfxdp(data, t, x, dfdp, im, isim);                     /* Stores dfxdp. Needs dvdp, dvdx and dvdu to be up to date */
            }
            
            /* Store number of iteration steps */
            storeIntegrationInfo( sim_mem, arcondition, ic );
            
            /**** end of CVODES ****/
        } else {
            /**** begin of SSA ****/
            /* (a) generate two randon numbers r1 and r2 uniformly distributed in (0,1)
             * (b) Compute alpha0 = \sum_{i=1}^q alpha_i(t),
             *     where alpha_i(t) = propensity function of i-th reaction
             * (c) The next reaction take place at time t+tau, where
             *     tau = 1/alpha0 * log(1/r1)   ... natural logrithm
             * (d) Determine which reaction occurs. Find j such that
             *     r2 >= 1/alpha0 \sum_{i=1}^{j-1} alpha_i(t)
             *     r2 <  1/alpha0 \sum_{i=1}^j     alpha_i(t)
             *     Update the number of reactants and products of the j-th reaction
             * (e) Go to (a) with t = t + tau 
             */
                        
            /* MATLAB values */
            double *x0 = mxGetData(mxGetField(arcondition, ic, "x0_ssa"));
            double *scale_x = mxGetData(mxGetField(arcondition, ic, "scale_x_ssa"));
            double *scale_v = mxGetData(mxGetField(arcondition, ic, "scale_v_ssa"));
            
            double *N = mxGetData(mxGetField(armodel, im, "N"));
            double *tlim = mxGetData(mxGetField(armodel, im, "tLim"));
            double *tfine = mxGetData(mxGetField(arcondition, ic, "tFine"));
            double *xssa = mxGetData(mxGetField(arcondition, ic, "xFineSSA"));
            double *xssa_lb = mxGetData(mxGetField(arcondition, ic, "xFineSSA_lb"));
            double *xssa_ub = mxGetData(mxGetField(arcondition, ic, "xFineSSA_ub"));
            int nx = (int) mxGetNumberOfElements(mxGetField(armodel, im, "xs"));
            int nt = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "tFine"));
            int nv = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "vNum"));
            
            if(has_tExp==1) {
                texp = mxGetData(mxGetField(arcondition, ic, "tExp"));
                ntexp = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "tExp"));
                xssaexp = mxGetData(mxGetField(arcondition, ic, "xExpSSA"));
            }
            
            /* Allocate state memory and user data memory */
            if ( allocateSimMemorySSA( sim_mem, nx ) )
            {
                /* Make some local pointer copies to facilitate handling */
                data = sim_mem->data;
                x = sim_mem->x;
                x_lb = sim_mem->x_lb;
                x_ub = sim_mem->x_ub;
            } else return;
            
            data->abort = abortSignal;
            data->u = mxGetData(mxGetField(arcondition, ic, "uNum"));
            data->p = mxGetData(mxGetField(arcondition, ic, "pNum"));
            data->v = mxGetData(mxGetField(arcondition, ic, "vNum"));            
            
            /* nruns loop */
            for (iruns=0; iruns<nruns; iruns++) {
                it = 0;
                itexp = 0;
                t = tlim[0];
                tfin = tfine[nt-1];
                for (ix=0; ix<10; ix++) {
                    lasttau[ix] = 1;
                }
                
                for (ix=0; ix<nx; ix++) {
                    /* printf("%f ", x0[ix]); */
                    Ith(x, ix+1) = x0[ix];
                    Ith(x_lb, ix+1) = x0[ix];
                    Ith(x_ub, ix+1) = x0[ix];
                }
                
                /* time loop */
                while((t<=tfin) & (it<nt)){
                    
                    /* (a) */
                    r1 = ((double)rand()/(double)RAND_MAX);
                    r2 = ((double)rand()/(double)RAND_MAX);
                    
                    if(r1==1){r1 = 0.5;}
                    if(r2==1){r2 = 0.5;}
                    
                    /* (b modified) */
                    fvSSA(data, t, x, im, ic);
                    alpha0 = 0;
                    for (iv=0; iv<nv; iv++) {
                        data->v[iv] = data->v[iv] * scale_v[iv];
                        alpha0 = alpha0 + data->v[iv];
                    }
                    alpha0 = 1/alpha0;
                    
                    /* (c) */
                    tau = alpha0 * log(1/r1);
                    lasttau[ilasttau] = tau;
                    ilasttau = (ilasttau+1) % 10;
                    
                    /* printf("r1=%g, r2=%g, alpha0=%g, tau=%g\n", r1, r2, alpha0, tau); */
                    if(tau<=0) {
                        printf("\nmodel #%i, condition #%i, run #%i at t=%f: STOP (tau=%g < 0)\n", im+1, ic+1, iruns+1, t, tau);
                        break;
                    }
                    
                    /* (d) */
                    iv = 0;
                    sumalpha = alpha0*data->v[0];
                    while(r2 >= sumalpha){
                        iv = iv + 1;
                        sumalpha = sumalpha + (alpha0*data->v[iv]);
                    }
                    
                    /* update values before reaction occurs */
                    while((tfine[it]<(t+tau)) & (it<nt)){
                        for (ix=0; ix<nx; ix++) {
                            xssa[it+nt*ix+nt*nx*iruns] = Ith(x, ix+1);
                            xssa_lb[it+nt*ix+nt*nx*iruns] = Ith(x_lb, ix+1);
                            xssa_ub[it+nt*ix+nt*nx*iruns] = Ith(x_ub, ix+1);
                        }
                        it += 1;
                        
                        for (ix=0; ix<nx; ix++) {
                            Ith(x_lb, ix+1) = Ith(x, ix+1);
                            Ith(x_ub, ix+1) = Ith(x, ix+1);
                        }
                    }
                    if((has_tExp==1) & (ntexp>0)) {
                        while((itexp<ntexp) & (texp[itexp]<(t+tau))){
                            for (ix=0; ix<nx; ix++) {
                                xssaexp[itexp+ntexp*ix+ntexp*nx*iruns] = Ith(x, ix+1);
                            }
                            itexp += 1;
                        }
                    }
                    
                    /* i-th reaction occurs */
                    for (ix=0; ix<nx; ix++) {
                        Ith(x, ix+1) = Ith(x, ix+1) + (N[ix+iv*nx] / scale_x[ix]);
                        if(Ith(x, ix+1)<0) Ith(x, ix+1) = 0;
                        if(Ith(x, ix+1)<Ith(x_lb, ix+1)) Ith(x_lb, ix+1) = Ith(x, ix+1);
                        if(Ith(x, ix+1)>Ith(x_ub, ix+1)) Ith(x_ub, ix+1) = Ith(x, ix+1);
                    }
                    
                    meantau = 0;
                    for (ix=0; ix<10; ix++) meantau += lasttau[ix];
                    meantau /= 10;
                    
                    if(meantau < mintau) {
                        printf("\nmodel #%i, condition #%i, run #%i at t=%f: STOP (mean(tau)=%g < %g)\n", im+1, ic+1, iruns+1, t, meantau, mintau);
                        break;
                    }
                    
                    /* (e) */
                    t = t + tau;
                }
                
                /* Terminate simulation when aborted */
                if ( *abortSignal == 1 )
                    break;
            }            
            /**** end of SSA ****/
        }
        
        DEBUGPRINT0( debugMode, 5, "Finished simulation\n" );
        
        /* call z_calc */
        DEBUGPRINT0( debugMode, 5, "Calling z-calc\n" );

        z_calc(im, ic, isim, arcondition, sensi);
    }

    gettimeofday(&t3, NULL);
    
    /* printf("computing model #%i, condition #%i (done)\n", im, ic); */
    /* call y_calc */
    DEBUGPRINT0( debugMode, 5, "Evaluating observations\n" );
    evaluateObservations(arcondition, im, ic, ysensi, has_tExp);
    
    gettimeofday(&t4, NULL);
    timersub(&t2, &t1, &tdiff);
    ticks_start[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
    timersub(&t3, &t1, &tdiff);
    ticks_stop_data[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
    timersub(&t4, &t1, &tdiff);
    ticks_stop[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
    
    DEBUGPRINT0( debugMode, 5, "Terminating x_calc\n" );
    
    /* Clean up */
    terminate_x_calc( sim_mem, *status );
}

void storeSimulation( UserData data, int im, int isim, int is, int nu, int nv, int neq, int nout, N_Vector x, double *returnx, double *returnu, double *returnv, double *qpositivex )
{
    int js;
    fu(data, data->t, im, isim);
    fv(data, data->t, x, im, isim);
                    
    for(js=0; js < nu; js++) returnu[js*nout+is] = data->u[js];
    for(js=0; js < nv; js++) returnv[js*nout+is] = data->v[js];
    copyStates( x, returnx, qpositivex, neq, nout, is );
}

void storeSensitivities( UserData data, int im, int isim, int is, int np, int nu, int nv, int neq, int nout, N_Vector x, N_Vector *sx, double *returnsx, double *returnsu, double *returnsv, int sensitivitySubset, int32_T *sensitivityMapping )
{
    int js, jss, ks;
    
    fsu(data, data->t, im, isim);
    fsv(data, data->t, x, im, isim);                        

    if ( sensitivitySubset == 0 )
    {
        /*****************************
        ** RETURN ALL SENSITIVITIES **
        ******************************/
        if(neq>0) {
            /* Output state sensitivities */
            copyNVMatrixToDouble( sx, returnsx, np, neq, nout, is );

            for(js=0; js < np; js++) {
                /* Output flux sensitivities */
                csv(data->t, x, js, sx[js], data, im, isim);
                for(ks=0; ks < nv; ks++) {
                    returnsv[(js*nv+ks)*nout + is] = data->sv[ks];
                }
            }
        }

        /* Output input sensitivities */
        for(js=0; js < np; js++) {
            for(ks=0; ks < nu; ks++) {
                returnsu[(js*nu+ks)*nout + is] = data->su[(js*nu)+ks];
            }
        }
    } else {
        /********************************************
        ** RETURN ONLY SUBSET OF THE SENSITIVITIES **
        *********************************************/
        if(neq>0) {

            /* Output state sensitivities */
            subCopyNVMatrixToDouble( sx, returnsx, np, neq, nout, is, sensitivityMapping );
            for(jss=0; jss < np; jss++) {
                js = sensitivityMapping[jss];                                   
                if (js < 0) {
                    for(ks=0; ks < nv; ks++) returnsv[(jss*nv+ks)*nout + is] = 0.0;
                } else {
                    /* Output flux sensitivities */
                    csv(data->t, x, js, sx[js], data, im, isim);
                    for(ks=0; ks < nv; ks++) {
                        returnsv[(jss*nv+ks)*nout + is] = data->sv[ks];
                    }
                }
            }
        }
        /* Output input sensitivities */
        for(js=0; js < np; js++) {
            for(ks=0; ks < nu; ks++) {
                returnsu[(js*nu+ks)*nout + is] = data->su[(js*nu)+ks];
            }
        }
    }   
}

/* Root finding procedures */
/* Two rootfinding procedures have been implemented */
/* The first is to simply apply the initial condition, store intermediate arrays and terminate immediately. In this case, the rootfinding is handled on the MATLAB side */
/* The second case is to do rootfinding within C++ (rootFinding = 2) */
void findRoots( SimMemory sim_mem, mxArray *arcondition, int im, int ic, int isim, double tstart, double eq_tol, int neq, int nu, int nv, int nout, double* returnx, double* returnu, double* returnv, double* qpositivex, double* returnsx, double* returnsu, double* returnsv, int sensi, int ysensi, int npSensi, int has_tExp )
{                
    double tEq = tstart;
    UserData data = sim_mem->data;

    /* No equations. Terminate now. */
    if ( neq == 0 ) { terminate_x_calc( sim_mem, 0 ); return; };

    DEBUGPRINT1( debugMode, 4, "Rootfinding mode at t=%g\n", tEq );

    if ( rootFinding == 2 )
    {
        #ifdef ROOT_FINDING
            solveSS( debugMode, arcondition, im, ic, isim, tEq, sim_mem->x, data, eq_tol );
            DEBUGPRINT0( debugMode, 4, "Root finding terminated\n" );
        #else
            DEBUGPRINT0( debugMode, 4, "Mex file was compiled without rootfinding ...\n" );
            terminate_x_calc( sim_mem, 0 ); return;
        #endif
    }

    /* Copy states and state sensitivities */
    if ( neq > 0 )
    {
        DEBUGPRINT0( debugMode, 4, "Copying states\n" );
        copyStates( sim_mem->x, returnx, qpositivex, neq, nout, 0 );
        DEBUGPRINT0( debugMode, 4, "Copying sensis\n" );
        DEBUGPRINT4( debugMode, 4, "sx: %d, npSensi: %d, neq: %d, nout: %d\n", sim_mem->sx, npSensi, neq, nout );
        /*if ( sensi ) copyNVMatrixToDouble( sim_mem->sx, returnsx, npSensi, neq, nout, 0 );*/
        if ( sensi ) storeSensitivities( data, im, isim, 0, npSensi, nu, nv, neq, nout, sim_mem->x, sim_mem->sx, returnsx, returnsu, returnsv, 0, NULL );
    }
    DEBUGPRINT0( debugMode, 4, "Calculating z\n" );
    z_calc(im, ic, isim, arcondition, ysensi);
    DEBUGPRINT0( debugMode, 6, "Calculating fu\n" );
    fu(data, tEq, im, isim);
    DEBUGPRINT0( debugMode, 6, "Calculating fv\n" );
    fv(data, tEq, sim_mem->x, im, isim);
    DEBUGPRINT0( debugMode, 4, "Copying u and v to outputs\n" );
    copyResult( data->u, returnu, nu, nout, 0 );
    copyResult( data->v, returnv, nv, nout, 0 );

    {
        double *dfdx, *dfdp;
        dfdx = mxGetData(mxGetField(arcondition, ic, "dfdxNum"));
        dfdp = mxGetData(mxGetField(arcondition, ic, "dfdpNum"));
        fsv(data, tEq, sim_mem->x, im, isim);                  /* Updates dvdp, dvdu, dvdx */
        getdfxdx(im, isim, tEq, sim_mem->x, dfdx, data);       /* Updates dvdx and stores dfxdx */
        dfxdp(data, tEq, sim_mem->x, dfdp, im, isim);          /* Stores dfxdp. Needs dvdp, dvdx and dvdu to be up to date */                
    }

    DEBUGPRINT0( debugMode, 4, "Evaluating observations\n" );
    evaluateObservations(arcondition, im, ic, ysensi, has_tExp);
    DEBUGPRINT0( debugMode, 4, "Terminating ...\n" );
    terminate_x_calc( sim_mem, 0 );
}

/* Store some information regarding the integration */
void storeIntegrationInfo( SimMemory sim_mem, mxArray *arcondition, int ic )
{
    mxArray *stepField;
    int64_T* stepData;
    long int nsteps;
    
    if ( sim_mem && ( sim_mem->cvode_mem ) )
       CVodeGetNumSteps( sim_mem->cvode_mem, &nsteps );
    else
        nsteps = -1;
            
    stepField = mxGetField( arcondition, ic, "stepsTaken" );
    if ( stepField )
    {
        stepData = (int64_T*)mxGetData( stepField );
        stepData[0] = (int64_T) nsteps;
    }
}     

int safeGetToggle( mxArray *data, int idx, const char *fieldName )
{
    mxArray *field = mxGetField(data, idx, fieldName);
    if ( !field || ( mxIsEmpty(field) ) ) {
        return 0;
	} else {
        return (int) mxGetScalar(field);
	}
}

void evaluateObservations( mxArray *arcondition, int im, int ic, int sensi, int has_tExp )
{
    mxArray *ardata;
    mxArray *dLink;
    double  *dLinkints; 
    int     id, nd, ids;
   
    DEBUGPRINT1( debugMode, 4, "Evaluating observations for condition %d\n", ic );
    ardata = mxGetField(armodel, im, "data");
	if(ardata!=NULL){
        dLink = mxGetField(arcondition, ic, "dLink");
        dLinkints = mxGetData(dLink);
        nd = (int) mxGetNumberOfElements(dLink);
        /* loop over data */
        for(ids=0; ids<nd; ++ids){
            id = ((int) dLinkints[ids]) - 1;
            DEBUGPRINT2( debugMode, 4, "Evaluating data with idx %d for condition %d\n", id, ic );            
            has_tExp = safeGetToggle(ardata, id, "has_tExp");
                        
            if((has_tExp == 1) | (fine == 1)) {
                y_calc(im, id, ardata, arcondition, sensi);
            }
        }
    }
    DEBUGPRINT1( debugMode, 4, "Finished evaluating observations for condition %d\n", ic );
}

void copyStates( N_Vector x, double *returnx, double *qpositivex, int neq, int nout, int offset )
{
    int js;
    
	for(js=0; js < neq; js++) {
        returnx[js*nout+offset] = Ith(x, js+1);
        /* set negative values to zeros */
        if(qpositivex[js]>0.5 && returnx[js*nout+offset]<0.0) returnx[js*nout+offset] = 0.0;
	}
}

void copyResult( double* data, double *returnvec, int nu, int nout, int offset )
{
    int js;
    for(js=0; js < nu; js++) returnvec[js*nout+offset] = data[js];
}

/* Copies doubles stored in NVector* array matrix into double array with specified offset */
void copyNVMatrixToDouble( N_Vector* sx, double *returnsx, int nps, int neq, int nout, int offset )
{
    int js, ks;
    realtype* sxtmp;
    
	for(js=0; js < nps; js++) {
        sxtmp = NV_DATA_S(sx[js]);
        for(ks=0; ks < neq; ks++) {
            returnsx[(js*neq+ks)*nout + offset] = sxtmp[ks];
        }
	}
}

/* Copies a subset of doubles stored in NVector* array matrix into double array with specified offset */
void subCopyNVMatrixToDouble( N_Vector* sx, double *returnsx, int nps, int neq, int nout, int offset, int32_T* targetIdx )
{
    int js, ks, jss;
    realtype* sxtmp;
    
    for(jss=0; jss < nps; jss++)
    {
        js = targetIdx[jss];
        if ( js < 0 )
        {
            for(ks=0; ks < neq; ks++) {
                returnsx[(jss*neq+ks)*nout + offset] = 0;
            }            
        } else {
            sxtmp = NV_DATA_S(sx[js]);
            for(ks=0; ks < neq; ks++) {
                returnsx[(jss*neq+ks)*nout + offset] = sxtmp[ks];
            }
        }
    }
}

/* Equilibrate the system until the RHS is under a specified threshold */
int equilibrate(void *cvode_mem, UserData data, N_Vector x, realtype t, double *equilibrated, double *returndxdt, double *teq, int neq, int im, int ic, int *abortSignal ) {
    int    i;
    int    step;
    int    flag;
    double time;
    double current_stepsize;
    bool   converged;

    flag = 0;
    step = 0;
    current_stepsize = init_eq_step;
    converged = false;
   
    /* Set the time to the last succesful time step */
    time = data->t;
    DEBUGPRINT1( debugMode, 8, "Going into equilibration at time %g\n", time );
    while( !converged )
    {        
        time = time + current_stepsize;

        /* Simulate up to next checkpoint */
        CVodeSetStopTime(cvode_mem, RCONST(time));
        flag = CVode(cvode_mem, RCONST(time), x, &t, CV_NORMAL);

        DEBUGPRINT3( debugMode, 8, "Equilibrating ... (Step %d, current time: %g, target time %g)\n", step, t, time );
        
        /* Abort on integration failure */
        if ( flag < 0 )
        {
            converged = true;
        }

        /* Store dxdt */
        fx(time, x, returndxdt, data, im, ic);

        converged = true;
        if ( equilibrated )
        {
            for (i=0; i<neq; i++)
                converged = ( converged && ( (equilibrated[i] < 0.1) || ( (fabs(returndxdt[i])<eq_tol) || (fabs(returndxdt[i]) < fabs(eq_rtol * Ith(x, i+1))) ) ) );
        } else {
            for (i=0; i<neq; i++)
                converged = ( converged && ( (fabs(returndxdt[i])<eq_tol) || (fabs(returndxdt[i]) < fabs(eq_rtol * Ith(x, i+1))) ) );
        }
        
        /* Oh no, we didn't make it! Terminate anyway. */
        if ( step > max_eq_steps )
        {
            flag = 20;
            converged = true;
        }

        /* Increase step size */
        step++;
        current_stepsize = current_stepsize * eq_step_factor;
    }
    *teq = time;

    return flag;
}

/* This function can be used to display errors from the threaded environment 
   mexErrMsgTxt crashes on R2013b when called from a thread */
void thr_error( const char* msg ) {
    printf( msg );
}

/* Event handler */
/* Put functions that are supposed to be evaluated on events here */
int handle_event( SimMemory sim_mem, int sensi_meth, int reinitSolver )
{
    int nps     = sim_mem->np;
    int neq     = sim_mem->neq;
    int sensi   = sim_mem->sensi;
    int32_T* idx;
    
    void* cvode_mem         = sim_mem->cvode_mem;
    EventData event_data    = sim_mem->event_data;
    N_Vector x              = sim_mem->x;
    N_Vector* sx            = sim_mem->sx;
    
    double A, B;
	int state, pars, flag, cStep, tStep;
	realtype* sxtmp;

    flag = 0;
    
	/* Current events */
	cStep = event_data->i;

	/* Total events */
    tStep = event_data->n;

	/* Do we override state and sensitivity values? */
	if (event_data->overrides==1)
    {
        /* Override state variables */
        for(state=0; state<neq; state++) {
            A = event_data->value_A[state*tStep+cStep];
            B = event_data->value_B[state*tStep+cStep];
            Ith(x, state+1) = A * Ith(x, state+1) + B;
        }
   
        /* printf("t[%d/%d]=%f  A: %f, B: %f\n", cStep, tStep, event_data->t[event_data->i], A, B); */
        /* Override sensitivity equations */
        if (sensi==1) {
            if ( sim_mem->data->sensIndices == NULL )
            {
                /* Simulate with all sensitivities */
                for (pars=0; pars<nps; pars++) {
                    if (neq>0) {
                        sxtmp = NV_DATA_S(sx[pars]);
                        for (state=0; state<neq; state++)
                        {
                            A = event_data->sensValue_A[(pars*neq+state)*tStep+cStep];
                            B = event_data->sensValue_B[(pars*neq+state)*tStep+cStep];
                            sxtmp[state] = A * sxtmp[state] + B;
                        }
                    }
                }
            } else
            {
                /* Simulate with a subset of the sensitivities */
                idx = sim_mem->data->sensIndices;
                for (pars=0; pars<sim_mem->npSensi; pars++) {
                    if (neq>0) {
                        sxtmp = NV_DATA_S(sx[pars]);
                        for (state=0; state<neq; state++)
                        {
                            A = event_data->sensValue_A[(idx[pars]*neq+state)*tStep+cStep];
                            B = event_data->sensValue_B[(idx[pars]*neq+state)*tStep+cStep];
                            sxtmp[state] = A * sxtmp[state] + B;
                        }
                    }
                }                
            }
        }
	}

    /* Reinitialize the solver */
    if ( reinitSolver == 1 )
    {
        flag = CVodeReInit(cvode_mem, RCONST(event_data->t[event_data->i]), x);
        if (flag>=0) {
            if (sensi==1) {
                flag = CVodeSensReInit(cvode_mem, sensi_meth, sx);
            }
        }
    }
        
    return flag;
}

/* Apply initial conditions for solving using numerical ODE integration */
int applyInitialConditionsODE( SimMemory sim_mem, double tstart, int im, int isim, double *returndxdt, double *returndfdp0, mxArray *x0_override, int sensitivitySubset )
{
    int nPoints;
    UserData data = sim_mem->data;
    int sensi = sim_mem->sensi;
    int neq = sim_mem->neq;
    int nps = sim_mem->np;
    N_Vector x = sim_mem->x;
    N_Vector *sx = sim_mem->sx;
    double *override;
    
	int is, js, ks;
	realtype *sxtmp;

	fu(data, tstart, im, isim);
	
	if ( neq > 0 )
	{
		for (is=0; is<neq; is++) Ith(x, is+1) = 0.0;
		fx0(x, data, im, isim);
        
        /* Override initial condition */
        if ( x0_override ) {
            nPoints = (int) mxGetNumberOfElements( x0_override );
            if ( nPoints > 0 ) {
                if ( nPoints != neq ) { terminate_x_calc( sim_mem, 21 ); return 0; };
                override = (double *) mxGetData(x0_override);
                for (is=0; is<neq; is++) Ith(x, is+1) = override[is];
            }
        }
        
		fv(data, tstart, x, im, isim);
		fx(tstart, x, returndxdt, data, im, isim);

		if (sensi == 1)
		{
			fsu(data, tstart, im, isim);

            if ( sensitivitySubset == 1 ) {
    			for(js=0; js < sim_mem->npSensi; js++) {
    				sxtmp = NV_DATA_S(sx[js]);
    				for(ks=0; ks < neq; ks++) {
    					sxtmp[ks] = 0.0;
    				}
                }                
                for (is=0;is<sim_mem->npSensi;is++) fsx0(is, sx[is], data, im, isim, sensitivitySubset);
            } else {
    			for(js=0; js < nps; js++) {
    				sxtmp = NV_DATA_S(sx[js]);
    				for(ks=0; ks < neq; ks++) {
    					sxtmp[ks] = 0.0;
    				}
                }
                for (is=0;is<nps;is++) fsx0(is, sx[is], data, im, isim, sensitivitySubset);
            }
			fsv(data, tstart, x, im, isim);
			dfxdp0(data, tstart, x, returndfdp0, im, isim);
		}
	}
    return 1;
}

/* Allocate memory used by the SUNDIALS solver */
int allocateSimMemoryCVODES( SimMemory sim_mem, int neq, int np, int sensi, int npSensi )
{
    int is, js, ks;
    realtype *atolV_tmp;
    
    sim_mem->neq = neq;
    sim_mem->np = np;
    sim_mem->npSensi = npSensi;
    sim_mem->sensi = sensi;
    
    /* Allocate userdata */
    sim_mem->data = (UserData) malloc(sizeof *sim_mem->data);
    if (sim_mem->data == NULL) { terminate_x_calc( sim_mem, 1 ); return 0; }
           
    /* Allocate event structure */
    sim_mem->event_data = (EventData) malloc(sizeof *sim_mem->event_data);
    if (sim_mem->event_data == NULL) { terminate_x_calc( sim_mem, 1 ); return 0; }    
    
    if ( neq > 0 ) {
        /* Create CVODES object */
        sim_mem->cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
        if (sim_mem->cvode_mem == NULL) { terminate_x_calc( sim_mem, 3 ); return 0; }        
        
        (sim_mem->x) = N_VNew_Serial(neq);
        if (sim_mem->x == NULL) {terminate_x_calc( sim_mem, 1 ); return 0; }

        /* Use private function to compute error weights */
        sim_mem->atolV = N_VNew_Serial(neq);
        if (sim_mem->atolV == NULL) { terminate_x_calc( sim_mem, 2 ); return 0; }
        
        for (is=0; is<neq; is++) 
            Ith(sim_mem->atolV, is+1) = 0.0;
        
        if (sensi == 1) {
            (sim_mem->sx) = N_VCloneVectorArray_Serial(npSensi, sim_mem->x);
            if (sim_mem->sx == NULL) { terminate_x_calc( sim_mem, 2 ); return 0; }

            sim_mem->atols_ss = N_VNew_Serial(npSensi);
            if (sim_mem->atols_ss == NULL) { terminate_x_calc( sim_mem, 2 ); return 0; }
            
            sim_mem->atolV_ss = N_VCloneVectorArray_Serial(npSensi, sim_mem->x);
            if (sim_mem->atolV_ss == NULL) { terminate_x_calc( sim_mem, 9 ); return 0; }
                    
            for(js=0; js < npSensi; js++) {
                atolV_tmp = NV_DATA_S(sim_mem->atolV_ss[js]);
                for(ks=0; ks < neq; ks++) {
                    atolV_tmp[ks] = 0.0;
                }
            }
        }
    }
    return 1;
}

/* Allocate memory for the states and sensitivities */
int allocateSimMemorySSA( SimMemory sim_mem, int nx )
{
    sim_mem->neq = nx;
    
    /* Allocate userdata */
    sim_mem->data = (UserData) malloc(sizeof *sim_mem->data);
    if (sim_mem->data == NULL) { terminate_x_calc( sim_mem, 1 ); return 0; }    
    
    if ( nx > 0 ) {
        sim_mem->x = N_VNew_Serial(nx);
        if (sim_mem->x == NULL) { terminate_x_calc( sim_mem, 2 ); return 0; }
        sim_mem->x_lb = N_VNew_Serial(nx);
        if (sim_mem->x_lb == NULL) { terminate_x_calc( sim_mem, 2 ); return 0; }
        sim_mem->x_ub = N_VNew_Serial(nx);
        if (sim_mem->x_ub == NULL) { terminate_x_calc( sim_mem, 2 ); return 0; }    
    }
    
    return 1;
}

/* Initialize the UserData structure for use with CVodes */
void initializeDataCVODES( SimMemory sim_mem, double tstart, int *abortSignal, mxArray *arcondition, double *qpositivex, int ic, int nsplines, int sensitivitySubset )
{
    int j;
    UserData data = sim_mem->data;
    data->splines = NULL;
    data->nsplines = 0;
    
    if ( nsplines > 0 ) {
        data->splines = (double**) malloc(nsplines * sizeof(double*));
        data->splineIndices = (int *) malloc(nsplines * sizeof(int));
        data->nsplines = nsplines;
        
        /* Initialize the pointers */
        for ( j = 0; j < nsplines; j++ )
            data->splines[j] = NULL;
        
        if (data->splines == NULL) { terminate_x_calc( sim_mem, 1 ); return; }   
    }
    
    if ( sensitivitySubset == 1 )
        data->sensIndices = (int32_T *) mxGetData(mxGetField(arcondition, ic, "sensIndices"));
	else
        data->sensIndices = NULL;
    
	data->abort = abortSignal;
	data->t = tstart;

	data->qpositivex = qpositivex;
	data->u = mxGetData(mxGetField(arcondition, ic, "uNum"));
	data->p = mxGetData(mxGetField(arcondition, ic, "pNum"));
	data->v = mxGetData(mxGetField(arcondition, ic, "vNum"));
	data->dvdx = mxGetData(mxGetField(arcondition, ic, "dvdxNum"));
	data->dvdu = mxGetData(mxGetField(arcondition, ic, "dvduNum"));
	data->dvdp = mxGetData(mxGetField(arcondition, ic, "dvdpNum"));

	if ( sim_mem->sensi == 1 ) {
        data->su = mxGetData(mxGetField(arcondition, ic, "suNum"));
        data->sv = mxGetData(mxGetField(arcondition, ic, "svNum"));
    }
}

int initializeEvents( SimMemory sim_mem, mxArray *arcondition, int ic, double tstart )
{
    int flag;
    int np = sim_mem->np;
    int neq = sim_mem->neq;
    int qEvents;
    EventData event_data = sim_mem->event_data;
    
	/* Initialize event list (points where solver needs to be reinitialized) */
	qEvents = init_list(arcondition, ic, tstart, &(event_data->n), &(event_data->t), &(event_data->i), "qEvents", "tEvents");        

	/* Allow state values and sensitivity values to be overwritten at events */
	event_data->overrides = 1;

	/* Grab additional data required for assignment operations */
	/* Assignment operations are of the form Ax+B where X is the state variable */
	flag = fetch_vector( arcondition, ic, &(event_data->value_A), "modx_A", neq*event_data->n );
	if ( flag < 0 ) { event_data->overrides = 0; };
	flag = fetch_vector( arcondition, ic, &(event_data->value_B), "modx_B", neq*event_data->n );
	if ( flag < 0 ) { event_data->overrides = 0; };
	flag = fetch_vector( arcondition, ic, &(event_data->sensValue_A), "modsx_A", neq*np*event_data->n );
	if ( flag < 0 ) { event_data->overrides = 0; };
	flag = fetch_vector( arcondition, ic, &(event_data->sensValue_B), "modsx_B", neq*np*event_data->n );
    if ( flag < 0 ) { event_data->overrides = 0; };
    
    return qEvents;
}

/* This function loads a vector/matrix from MATLAB and checks it against desired length */
int fetch_vector( mxArray* arcondition, int ic, double **vector, const char* fieldname, int desiredLength ) {
    
    mxArray *field;
    int nPoints;      
    
    field = mxGetField(arcondition, ic, fieldname);

    if ( field != NULL )
    {
        /* Check whether vector has the desired length */
        nPoints = (int) mxGetNumberOfElements( field );
        
        /* No points => silently drop */
        if (nPoints == 0)
           return -1;
        
        /* Wrong length => warn! */
        if (nPoints != desiredLength)
        {
           printf( "Warning, mod vector has incorrect size -> overrides disabled\n" );
           return -1;
        }

        /* Grab the data */
        *vector = (double*) mxGetData( field );
        return 1;
    } else {
        return -1;
    }
}

/* Free memory taken up by x_calc */
void terminate_x_calc( SimMemory sim_mem, double status )
{
    /* Something is seriously wrong */
	if ( sim_mem == NULL )
    {
        mexPrintf( "FATAL ERROR: Simulation memory is null upon terminate_x_calc!" );
		return;
    }

	/* Report status to user */
	sim_mem->status[0] = status;
  
	/* Make sure the thread terminates */
	sim_mem->threadStatus[0] = 1;
    
    /* Free the memory that was allocated */
	simFree( sim_mem );
}

/* This function initializes time point lists */
int init_list( mxArray* arcondition, int ic, double tstart, int* nPoints, double** timePoints, int* currentIndex, const char* flagFieldName, const char* timePointFieldName ) {
    int ID, flag;
    double *time;
          
    flag = (int) mxGetScalar(mxGetField(arcondition, ic, flagFieldName));
    if (flag==1) {
        mxArray *timePointField = mxGetField(arcondition, ic, timePointFieldName);

        if ( timePointField != NULL ) {
             time = (double*) mxGetData( timePointField );
             *nPoints = (int) mxGetNumberOfElements( timePointField );

             /* Move past pre-simulation points */
             ID = 0;
             while((time[ID] < tstart)&&(ID<(*nPoints)))
                 ID++;

             /* Set output values */
             *currentIndex = ID;
             *timePoints = time;
        } else { 
             thr_error( "Cannot find time point field. Not using events.\n" );
             return 0;
        }
    }

	 return flag;
}

/* calculate derived variables */
void z_calc(int im, int ic, int isim, mxArray *arcondition, int sensi) {
    
    /* printf("computing model #%i, condition #%i, derived variables\n", im, ic); */
    
  int nt, np, nx;
    int it;
            
    double *t;
    
    double *p;
    double *u;
    double *x;
    double *z;
    double *su;
    double *sx;
    double *sz;
    double *dzdx;
    
    /* MATLAB values */
    if(fine == 1){
        t = mxGetData(mxGetField(arcondition, ic, "tFine"));
        nt = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "tFine"));
        
        u = mxGetData(mxGetField(arcondition, ic, "uFineSimu"));
        x = mxGetData(mxGetField(arcondition, ic, "xFineSimu"));
        z = mxGetData(mxGetField(arcondition, ic, "zFineSimu"));
        if (sensi == 1) {
            su = mxGetData(mxGetField(arcondition, ic, "suFineSimu"));
            sx = mxGetData(mxGetField(arcondition, ic, "sxFineSimu"));
            sz = mxGetData(mxGetField(arcondition, ic, "szFineSimu"));
        }
    }
    else{
        t = mxGetData(mxGetField(arcondition, ic, "tExp"));
        nt = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "tExp"));
        
        u = mxGetData(mxGetField(arcondition, ic, "uExpSimu"));
        x = mxGetData(mxGetField(arcondition, ic, "xExpSimu"));
        z = mxGetData(mxGetField(arcondition, ic, "zExpSimu"));
	dzdx = mxGetData(mxGetField(arcondition, ic, "dzdx"));
        if (sensi == 1) {
            su = mxGetData(mxGetField(arcondition, ic, "suExpSimu"));
            sx = mxGetData(mxGetField(arcondition, ic, "sxExpSimu"));
            sz = mxGetData(mxGetField(arcondition, ic, "szExpSimu"));
        }
    }
    p = mxGetData(mxGetField(arcondition, ic, "pNum"));
    np = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "pNum"));
    nx = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "dxdt"));
    
    /* loop over output points */
    for (it=0; it < nt; it++) {
        /* printf("%f y-loop (im=%i id=%i)\n", t[it], im, id); */
        
        fz(t[it], nt, it, 0, 0, 0, 0, z, p, u, x, im, isim);
	if( fine == 0 ) {
	  dfzdx(t[it], nt, it, 0, nx, 0, 0, dzdx, z, p, u, x, im, isim);
	}
        if (sensi == 1) {
            fsz(t[it], nt, it, np, sz, p, u, x, z, su, sx, im, isim);
        }
    }
    
    /* printf("computing model #%i, condition #%i, derived variables (done)\n", im, ic); */
}

/* calculate observations */
void y_calc(int im, int id, mxArray *ardata, mxArray *arcondition, int sensi) {
    /* printf("computing model #%i, data #%i\n", im, id); */
    
    int nt, it, iy, ip, ntlink, itlink;
    int ny, np, ic;
    int has_yExp;
    
    double *t;
    double *tlink;
    
    double *y;
    double *ystd;
    double *yexp;
    
    double *sy;
    double *systd;
    
    double *qlogy;
    
    double *p;
    double *u;
    double *x;
    double *z;
    double *su;
    double *sx;
    double *sz;
    double *y_scale;
    double *dzdx;
    
    
    /* MATLAB values */
    ic = (int) mxGetScalar(mxGetField(ardata, id, "cLink")) - 1;
    has_yExp = safeGetToggle(ardata, id, "has_yExp");
    
    ny = (int) mxGetNumberOfElements(mxGetField(ardata, id, "y"));
    qlogy = mxGetData(mxGetField(ardata, id, "logfitting"));
    /*qlogp = mxGetData(mxGetField(ardata, id, "qLog10"));*/
    p = mxGetData(mxGetField(ardata, id, "pNum"));
    np = (int) mxGetNumberOfElements(mxGetField(ardata, id, "pNum"));
    
    if(fine == 1){
        t = mxGetData(mxGetField(ardata, id, "tFine"));
        nt = (int) mxGetNumberOfElements(mxGetField(ardata, id, "tFine"));
        tlink = mxGetData(mxGetField(ardata, id, "tLinkFine"));
        ntlink = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "tFine"));
        
        y = mxGetData(mxGetField(ardata, id, "yFineSimu"));
        ystd = mxGetData(mxGetField(ardata, id, "ystdFineSimu"));
        
        u = mxGetData(mxGetField(arcondition, ic, "uFineSimu"));
        x = mxGetData(mxGetField(arcondition, ic, "xFineSimu"));
        z = mxGetData(mxGetField(arcondition, ic, "zFineSimu"));
        
        if (sensi == 1) {
            sy = mxGetData(mxGetField(ardata, id, "syFineSimu"));
            systd = mxGetData(mxGetField(ardata, id, "systdFineSimu"));
            
            su = mxGetData(mxGetField(arcondition, ic, "suFineSimu"));
            sx = mxGetData(mxGetField(arcondition, ic, "sxFineSimu"));
            sz = mxGetData(mxGetField(arcondition, ic, "szFineSimu"));
        }
    }
    else{
        t = mxGetData(mxGetField(ardata, id, "tExp"));
        nt = (int) mxGetNumberOfElements(mxGetField(ardata, id, "tExp"));
        tlink = mxGetData(mxGetField(ardata, id, "tLinkExp"));
        ntlink = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "tExp"));
        
        y = mxGetData(mxGetField(ardata, id, "yExpSimu"));
        ystd = mxGetData(mxGetField(ardata, id, "ystdExpSimu"));
        
	y_scale = mxGetData(mxGetField(ardata, id, "y_scale"));
	dzdx = mxGetData(mxGetField(arcondition, ic, "dzdx"));
        u = mxGetData(mxGetField(arcondition, ic, "uExpSimu"));
        x = mxGetData(mxGetField(arcondition, ic, "xExpSimu"));
        z = mxGetData(mxGetField(arcondition, ic, "zExpSimu"));
        
        if (sensi == 1) {
            sy = mxGetData(mxGetField(ardata, id, "syExpSimu"));
            systd = mxGetData(mxGetField(ardata, id, "systdExpSimu"));
            
            su = mxGetData(mxGetField(arcondition, ic, "suExpSimu"));
            sx = mxGetData(mxGetField(arcondition, ic, "sxExpSimu"));
            sz = mxGetData(mxGetField(arcondition, ic, "szExpSimu"));
        }
        
        if (has_yExp == 1) {
            yexp = mxGetData(mxGetField(ardata, id, "yExp"));
        }
    }
    
    /* loop over output points */
    for (it=0; it < nt; it++) {
        /* printf("%f y-loop (im=%i id=%i)\n", t[it], im, id); */
        itlink = (int) tlink[it] - 1;
        
        fy(t[it], nt, it, ntlink, itlink, 0, 0, 0, 0, y, p, u, x, z, im, id);
        
	if(fine==0){
	  fy_scale(t[it], nt, it, ntlink, itlink, 0, 0, 0, 0, y_scale, p, u, x, z, dzdx, im, id);
	}

        /* log trafo of y */
        for (iy=0; iy<ny; iy++) {
            if(qlogy[iy] > 0.5){
                if(y[it + (iy*nt)]<-cvodes_atol) printf("WARNING, check for concentrations <= 0 in data %d and observable %d !!!\n", id+1, iy+1);
                if(fine==0)  y_scale[it + (iy*nt)] = y_scale[it + (iy*nt)] / y[it + (iy*nt)] / log(10.0);
                y[it + (iy*nt)] = log10(y[it + (iy*nt)]);
            }
        }
        
            fystd(t[it], nt, it, ntlink, itlink, ystd, y, p, u, x, z, im, id);

        if (sensi == 1) {
            fsy(t[it], nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz, im, id);
            
            /* log trafo of sy */
            for (iy=0; iy<ny; iy++) {
                if(qlogy[iy] > 0.5) {
                    for (ip=0; ip < np; ip++) {
                        sy[it + (iy*nt) + (ip*nt*ny)] =
                                sy[it + (iy*nt) + (ip*nt*ny)]
                                / pow(10.0, y[it + (iy*nt)])
                                / log(10.0);
                    }
                }
            }
            
            fsystd(t[it], nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz, im, id);
        }
        
        if ((has_yExp == 1) & (fine == 0)) {
            y_checkNaN(nt, ny, it, y, yexp, ystd);
        }
     
    }
    
    /* printf("computing model #%i, data #%i (done)\n", im, id); */
}

/* Cleaner alternative to debugPrint, but more overhead */
#ifndef MACRO_DEBUGPRINT
void debugPrint( int debugMode, int level, const char* format, ... )
{
    va_list args;

    if ( debugMode > level )
    {
        mexPrintf( "[D] " );
        va_start( args, format );
        mexPrintf( format, args );
        va_end( args );
    }    
}
#endif

/* if isnan(yexp) then set y and ystd also to NaN
 This part of the code was previously done in function fres */
void y_checkNaN(int nt, int ny, int it, double *y, double *yexp, double *ystd) {
    int iy;     
    for(iy=0; iy<ny; iy++){
        /* in case of missing data (nan) */
        if(mxIsNaN(yexp[it + (iy*nt)])) {
            y[it + (iy*nt)] = yexp[it + (iy*nt)];
            ystd[it + (iy*nt)] = yexp[it + (iy*nt)];
        }
    }
}

