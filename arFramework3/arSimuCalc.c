/*
 *  MATLAB usage: arSimuCalc(struct ar, int fine, int sensi)
 *
 *  Copyright Andreas Raue 2011 (andreas.raue@fdm.uni-freiburg.de)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

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
int    fiterrors;
int    cvodes_maxsteps;
double cvodes_maxstepsize;
int    cvodes_atolV;
int    cvodes_atolV_Sens;
double cvodes_rtol;
double cvodes_atol;
double fiterrors_correction;
int useFitErrorMatrix;
double *fiterrors_matrix;
mwSize nrows_fiterrors_matrix;

/* Name of the substructure with conditions we're currently evaluating */
char condition_name[MXSTRING];

/* Name of the threads substructure we are currently using */
char threads_name[MXSTRING];

/* Equilibration variables (Used when time point Inf is encountered) */
int    max_eq_steps;          /* Maximal equilibration steps */
double init_eq_step;          /* Initial equilibration stepsize attempt */
double eq_step_factor;        /* Factor with which to increase the time at each equilibration step */
double eq_tol;                /* Absolute tolerance for equilibration */

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
void x_calc(int im, int ic, int sensi, int setSparse, int *threadStatus, int *abortSignal);
void z_calc(int im, int ic, mxArray *arcondition, int sensi);
void y_calc(int im, int id, mxArray *ardata, mxArray *arcondition, int sensi);

void fres(int nt, int ny, int it, double *res, double *y, double *yexp, double *ystd, double *chi2, double fiterrors_correction_factor);
void fsres(int nt, int ny, int np, int it, double *sres, double *sy, double *yexp, double *ystd, double fiterrors_correction_factor);
void fres_error(int nt, int ny, int it, double *reserr, double *res, double *y, double *yexp, double *ystd, double *chi2);
void fsres_error(int nt, int ny, int np, int it, double *sres, double *sreserr, double *sy, double *systd, double *y, double *yexp, double *ystd, double *res, double *reserr);

int ewt(N_Vector y, N_Vector w, void *user_data);
void thr_error( const char* msg );
int fetch_vector( mxArray* arcondition, int ic, double **vector, const char* fieldname, int desiredLength );
int init_list( mxArray* arcondition, int ic, double tstart, int* nPoints, double** timePoints, int* currentIndex, const char* flagFieldName, const char* timePointFieldName );

/* user functions */
#include "arSimuCalcFunctions.c"

int handle_event(void *cvode_mem, EventData event_data, UserData user_data, N_Vector x, N_Vector *sx, int nps, int neq, int sensi, int sensi_meth );
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

    /* get ar.config */
    arconfig = mxGetField(prhs[0], 0, "config");
    parallel = (int) mxGetScalar(mxGetField(arconfig, 0, "useParallel"));
    jacobian = (int) mxGetScalar(mxGetField(arconfig, 0, "useJacobian"));
    setSparse = (int) mxGetScalar(mxGetField(arconfig, 0, "useSparseJac"));
    sensirhs = (int) mxGetScalar(mxGetField(arconfig, 0, "useSensiRHS"));
    cvodes_atolV = (int) mxGetScalar(mxGetField(arconfig, 0, "atolV"));
    cvodes_atolV_Sens = (int) mxGetScalar(mxGetField(arconfig, 0, "atolV_Sens"));
    cvodes_rtol = mxGetScalar(mxGetField(arconfig, 0, "rtol"));
    cvodes_atol = mxGetScalar(mxGetField(arconfig, 0, "atol"));
    cvodes_maxsteps = (int) mxGetScalar(mxGetField(arconfig, 0, "maxsteps"));
    cvodes_maxstepsize = mxGetScalar(mxGetField(arconfig, 0, "maxstepsize"));
    fiterrors = (int) mxGetScalar(mxGetField(arconfig, 0, "fiterrors"));
    fiterrors_correction = (double) mxGetScalar(mxGetField(arconfig, 0, "fiterrors_correction"));
    useFitErrorMatrix = (int) mxGetScalar(mxGetField(arconfig, 0, "useFitErrorMatrix"));
    if(useFitErrorMatrix==1){
        fiterrors_matrix = (double *) mxGetData(mxGetField(arconfig, 0, "fiterrors_matrix"));
        nrows_fiterrors_matrix = mxGetM(mxGetField(arconfig, 0, "fiterrors_matrix"));
    }
    mintau = mxGetScalar(mxGetField(arconfig, 0, "ssa_min_tau"));
    nruns = (int) mxGetScalar(mxGetField(arconfig, 0, "ssa_runs"));
    ms = (int) mxGetScalar(mxGetField(arconfig, 0, "useMS"));
    events = (int) mxGetScalar(mxGetField(arconfig, 0, "useEvents"));
    max_eq_steps = (int) mxGetScalar(mxGetField(arconfig, 0, "max_eq_steps"));
    init_eq_step = (double) mxGetScalar(mxGetField(arconfig, 0, "init_eq_step"));
    eq_step_factor = (double) mxGetScalar(mxGetField(arconfig, 0, "eq_step_factor"));
    eq_tol = (double) mxGetScalar(mxGetField(arconfig, 0, "eq_tol"));

    if ( ms == 1 ) events = 1;
        
    /* threads */
    arthread = mxGetField(arconfig, 0, threads_name);
    /*nthreads = (int) mxGetScalar(mxGetField(arconfig, 0, "nThreads"));*/
    nthreads = (int) mxGetNumberOfElements(arthread);
    
#ifdef HAS_PTHREAD
    if(NMAXTHREADS<nthreads) mexErrMsgTxt("ERROR at NMAXTHREADS < nthreads");
#endif
    
    /* printf("%i threads (%i fine, %i sensi, %i jacobian, %g rtol, %g atol, %i maxsteps)\n", nthreads, fine,
            sensi, jacobian, cvodes_rtol, cvodes_atol, cvodes_maxsteps); */
    
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
    
    for(in=0; in<n; ++in){
        /* printf("computing thread #%i, task %i/%i (m=%i, c=%i)\n", id, in, n, ms[in], cs[in]); */
        x_calc(ms[in], cs[in], globalsensi, setSparse, &threadStatus[id], &threadAbortSignal[id]);
    }
    
    /* printf("computing thread #%i(done)\n", id); */
    
#ifdef HAS_PTHREAD
    if(parallel==1) {pthread_exit(NULL);}
    return NULL;
#endif
}

/* calculate dynamics */
void x_calc(int im, int ic, int sensi, int setSparse, int *threadStatus, int *abortSignal) {
    mxArray    *arcondition;
    mxArray    *ardata;
    
    mxArray *dLink;
    double *dLinkints;      
    
    int nm, nc, id, nd, has_tExp, has_yExp;
    int flag;
    int is, js, ks, ids;
    int nout, neq, nyout;
    int nu, np, nps, nv, ny, nnz;
    
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
    N_Vector x;
    N_Vector atolV;
    N_Vector atols_ss;
    N_Vector *atolV_ss;
    realtype *atolV_tmp;
    N_Vector *sx;
    realtype *sxtmp;
    
    /* SSA variables */
    N_Vector x_lb;
    N_Vector x_ub;
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
    double *returnddxdtdp;
    double *equilibrated;

    double *y;
    double *yExp;
    double *yStd;
    double *y_scale;
    double *y_scale_S;
    double *y_max_scale;
    double *y_max_scale_S;

    struct timeval t2;
    struct timeval t3;
    struct timeval t4;
    struct timeval tdiff;
    double *ticks_start;
    double *ticks_stop_data;
    double *ticks_stop;
    
    int sensi_meth = CV_SIMULTANEOUS; /* CV_SIMULTANEOUS or CV_STAGGERED */
    bool error_corr = TRUE;
    
    /* printf("computing model #%i, condition #%i\n", im, ic); */
    only_sim = 0;
    
    /* Grab value of infinity (for steady state simulations) */
    inf = mxGetInf();    

    /* check if im in range */
    nm = (int) mxGetNumberOfElements(armodel);
    if(nm<=im) {
        thr_error("im > length(ar.model)\n");
        return;
    }
    
    /* get ar.model(im).condition */
    arcondition = mxGetField(armodel, im, condition_name);
           
    if(arcondition==NULL){
        return;
    }
    
    /* check if ic in range */
    nc = (int) mxGetNumberOfElements(arcondition);
    if(nc<=ic) {
        thr_error("ic > length(ar.model.condition)\n");
        return;
    }
    
    /* Get double handle to store equilibrium value */
    teq = mxGetData(mxGetField(arcondition, ic, "tEq"));
    
    has_tExp = (int) mxGetScalar(mxGetField(arcondition, ic, "has_tExp"));
    if(has_tExp == 0 && fine == 0) return;
    
    /* get ar.model(im).data */
    ardata = mxGetField(armodel, im, "data");

    ticks_start = mxGetData(mxGetField(arcondition, ic, "start"));
    ticks_stop = mxGetData(mxGetField(arcondition, ic, "stop"));
    ticks_stop_data = mxGetData(mxGetField(arcondition, ic, "stop_data"));

    gettimeofday(&t2, NULL);
    
    if(dynamics == 1) {
        if(ssa == 0) {
            /**** begin of CVODES ****/
            
            /* get MATLAB values */
            qpositivex = mxGetData(mxGetField(armodel, im, "qPositiveX"));
            status = mxGetData(mxGetField(arcondition, ic, "status"));
            tstart = mxGetScalar(mxGetField(arcondition, ic, "tstart"));
            neq = (int) mxGetNumberOfElements(mxGetField(armodel, im, "xs"));
            nnz = (int) mxGetScalar(mxGetField(armodel, im, "nnz"));
     
            if(fine == 1){
                ts = mxGetData(mxGetField(arcondition, ic, "tFine"));
                nout = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "tFine"));
                
                returnu = mxGetData(mxGetField(arcondition, ic, "uFineSimu"));
                returnv = mxGetData(mxGetField(arcondition, ic, "vFineSimu"));
                returnx = mxGetData(mxGetField(arcondition, ic, "xFineSimu"));
		y_max_scale = mxGetData(mxGetField(arcondition, ic, "y_atol"));
        y_max_scale_S = mxGetData(mxGetField(arcondition, ic, "y_atolS"));
                if (sensi == 1) {
                    returnsu = mxGetData(mxGetField(arcondition, ic, "suFineSimu"));
                    returnsv = mxGetData(mxGetField(arcondition, ic, "svFineSimu"));
                    returnsx = mxGetData(mxGetField(arcondition, ic, "sxFineSimu"));
                }
            }
            else{
                ts = mxGetData(mxGetField(arcondition, ic, "tExp"));
                nout = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "tExp"));
                
                returnu = mxGetData(mxGetField(arcondition, ic, "uExpSimu"));
                returnv = mxGetData(mxGetField(arcondition, ic, "vExpSimu"));
                returnx = mxGetData(mxGetField(arcondition, ic, "xExpSimu"));
		y_max_scale = mxGetData(mxGetField(arcondition, ic, "y_atol"));
        y_max_scale_S = mxGetData(mxGetField(arcondition, ic, "y_atolS"));
		/* Scaling part, take Residuals and y_scale from last iter */
	      if(ardata!=NULL && (cvodes_atolV ==1 || cvodes_atolV_Sens==1) && neq>0){  
		dLink = mxGetField(arcondition, ic, "dLink");
		dLinkints = mxGetData(dLink);
		nd = (int) mxGetNumberOfElements(dLink);
		/* loop over data */
		for(ids=0; ids<nd; ++ids){
		  id = ((int) dLinkints[ids]) - 1;
		  has_yExp = (int) mxGetScalar(mxGetField(ardata, id, "has_yExp"));
		  if(has_yExp == 1) {
		    y = mxGetData(mxGetField(ardata, id, "yExpSimu"));
		    ny = (int) mxGetNumberOfElements(mxGetField(ardata, id, "y"));
		    yExp = mxGetData(mxGetField(ardata, id, "yExp"));
		    yStd = mxGetData(mxGetField(ardata, id, "ystdExpSimu"));
		    nyout = (int) mxGetNumberOfElements(mxGetField(ardata, id, "tExp"));

            if( (useFitErrorMatrix == 0 && fiterrors == -1) || (useFitErrorMatrix == 1 && fiterrors_matrix[id*nrows_fiterrors_matrix+im] == -1) ) {
                yStd = mxGetData(mxGetField(ardata, id, "yExpStd"));
            }

		    y_scale = mxGetData(mxGetField(ardata, id, "y_scale"));
		    y_scale_S = mxGetData(mxGetField(ardata, id, "y_scale_S"));

		    for(is=0; is<neq; is++){
		      for(js=0; js<nyout; js++){
			for(ks=0; ks<ny; ks++){
			   
			  if(!mxIsNaN(yExp[js + (ks*nyout)]) && !mxIsNaN(y[js + (ks*nyout)]) && !mxIsNaN(yStd[js + (ks*nyout)]) && yStd[js + (ks*nyout)]>0.) {
			    if(useFitErrorMatrix == 1 && fiterrors_matrix[id*nrows_fiterrors_matrix+im] != 1) {
			      y_scale_S[js+ks*nyout+is*nyout*ny] = y_scale[js+ks*nyout+is*nyout*ny] * 2* fabs(yExp[js + (ks*nyout)] - y[js + (ks*nyout)]) / pow(yStd[js + (ks*nyout)],2);
			    } else {
			      y_scale_S[js+ks*nyout+is*nyout*ny] = y_scale[js+ks*nyout+is*nyout*ny] * 2* fabs(yExp[js + (ks*nyout)] - y[js + (ks*nyout)]) / pow(yStd[js + (ks*nyout)],2) * sqrt(fiterrors_correction);
			    }
			  }

			  if(fabs(y_scale[js+ks*nyout+is*nyout*ny])>y_max_scale[is] && !mxIsNaN(y_scale[js+ks*nyout+is*nyout*ny]))
			    y_max_scale[is] = fabs(y_scale[js+ks*nyout+is*nyout*ny]);
              
              if(fabs(y_scale_S[js+ks*nyout+is*nyout*ny])>y_max_scale_S[is] && !mxIsNaN(y_scale_S[js+ks*nyout+is*nyout*ny]))
			    y_max_scale_S[is] = fabs(y_scale_S[js+ks*nyout+is*nyout*ny]);
			  /*printf("y_scale old = %f and scale for neq %i, t %i, y %i, thus %i is = %f \n", y_max_scale[is], is, js, ks, js+ks*nout+is*nout*ny, y_scale[js+ks*nout+is*nout*ny]); */
			  
			}
		      }
		    }
		    
		  }
		  
		}
	      }
                if (sensi == 1) {
                    returnsu = mxGetData(mxGetField(arcondition, ic, "suExpSimu"));
                    returnsv = mxGetData(mxGetField(arcondition, ic, "svExpSimu"));
                    returnsx = mxGetData(mxGetField(arcondition, ic, "sxExpSimu"));
                }
            }
            
            returndxdt = mxGetData(mxGetField(arcondition, ic, "dxdt"));
                      
            if (sensi == 1) {
                returnddxdtdp = mxGetData(mxGetField(arcondition, ic, "ddxdtdp"));
            }

            /* User data structure */
            data = (UserData) malloc(sizeof *data);
            if (data == NULL) {status[0] = 1; return;}
            data->abort = abortSignal;
            data->t = tstart;
            
            /* Event structure */
            event_data = (EventData) malloc(sizeof *event_data);
            if (event_data == NULL) {status[0] = 18; return;}

            /* Initialize multiple shooting list */
            if (ms==1) 
               qMS = init_list(arcondition, ic, tstart, &(event_data->nMS), &(event_data->tMS), &(event_data->iMS), "qMS", "tMS");
            
            /* Initialize userdata and derivatives */
            data->qpositivex = qpositivex;
            data->u = mxGetData(mxGetField(arcondition, ic, "uNum"));
            nu = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "uNum"));
            
            data->p = mxGetData(mxGetField(arcondition, ic, "pNum"));
            np = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "pNum"));
            nps = np;
            
            /* If there are no parameters, do not compute sensitivities; otherwise failure at N_VCloneVectorArray_Serial */
            if (nps==0) sensi = 0;
            
            data->v = mxGetData(mxGetField(arcondition, ic, "vNum"));
            nv = (int) mxGetNumberOfElements(mxGetField(arcondition, ic, "vNum"));
            data->dvdx = mxGetData(mxGetField(arcondition, ic, "dvdxNum"));
            data->dvdu = mxGetData(mxGetField(arcondition, ic, "dvduNum"));
            data->dvdp = mxGetData(mxGetField(arcondition, ic, "dvdpNum"));

            /* Initialize event list (points where solver needs to be reinitialized) */
            qEvents = 0;
            if (events==1)
            {
                qEvents = init_list(arcondition, ic, tstart, &(event_data->n), &(event_data->t), &(event_data->i), "qEvents", "tEvents");        

                /* Allow state values and sensitivity values to be overwritten at events */
                event_data->overrides = 1;

                /* Grab additional data required for assignment operations */
                /* Assignment operations are of the form Ax+B where X is the state variable */
                flag = fetch_vector( arcondition, ic, &(event_data->value_A), "modx_A", neq*event_data->n );
                if ( flag < 0 ) { event_data->overrides = 0; };
                flag = fetch_vector( arcondition, ic, &(event_data->value_B), "modx_B", neq*event_data->n );
                if ( flag < 0 ) { event_data->overrides = 0; };
                flag = fetch_vector( arcondition, ic, &(event_data->sensValue_A), "modsx_A", neq*nps*event_data->n );
                if ( flag < 0 ) { event_data->overrides = 0; };
                flag = fetch_vector( arcondition, ic, &(event_data->sensValue_B), "modsx_B", neq*nps*event_data->n );
                if ( flag < 0 ) { event_data->overrides = 0; };
            }
            
            /* Override which condition to simulate */
            /* This is used for equilibration purposes */
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
            
            /* fill for t=0 */
            fu(data, tstart, im, isim);
            
            if(neq>0){
                /* Initial conditions */
                x = N_VNew_Serial(neq);
                if (x == NULL) {status[0] = 2; return;}
                for (is=0; is<neq; is++) Ith(x, is+1) = 0.0;
                fx0(x, data, im, isim);
                fv(data, tstart, x, im, isim);
                fx(tstart, x, returndxdt, data, im, isim);

                /* Create CVODES object */
                cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
                if (cvode_mem == NULL) {status[0] = 3; return;}
                
                /* Allocate space for CVODES */
                flag = AR_CVodeInit(cvode_mem, x, tstart, im, isim);
                if (flag < 0) {status[0] = 4; return;}
                
                /* Number of maximal internal steps */
                flag = CVodeSetMaxNumSteps(cvode_mem, cvodes_maxsteps);
                if(flag < 0) {status[0] = 15; return;}
                
                /* Maximal internal step size */
                flag = CVodeSetMaxStep(cvode_mem, cvodes_maxstepsize);
                if(flag < 0) {status[0] = 19; return;}
                
                /* Use private function to compute error weights */
		atolV = N_VNew_Serial(neq);
		if (atolV == NULL) {status[0] = 2; return;}
		for (is=0; is<neq; is++) Ith(atolV, is+1) = 0.0;
		 
		if(cvodes_atolV==1)   { 		
            double tmp_tol = 1.;
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
		}else{                
		  flag = CVodeSStolerances(cvode_mem, RCONST(cvodes_rtol), RCONST(cvodes_atol));
		}
                if (flag < 0) {status[0] = 5; return;}
                
                /* Attach user data */
                flag = CVodeSetUserData(cvode_mem, data);
                if (flag < 0) {status[0] = 6; return;}
                
                /* Attach linear solver */
                if(setSparse == 0){
                    /* Dense solver */
                    flag = CVDense(cvode_mem, neq);
                }else{              
                    /* sparse linear solver KLU */
                    flag = CVKLU(cvode_mem, neq, nnz);
                }
                if (flag < 0) {status[0] = 7; return;}
                
                /* Jacobian-related settings */
                if (jacobian == 1) {
                    flag = AR_CVDlsSetDenseJacFn(cvode_mem, im, isim, setSparse);
                    if (flag < 0) {status[0] = 8; return;}
                }
                
                /* custom error weight function */
                /*
            flag = CVodeWFtolerances(cvode_mem, ewt);
            if (flag < 0) return;
                 */
            }
            
            /* Sensitivity-related settings */
            if (sensi == 1) {
                /* User data structure */
                data->su = mxGetData(mxGetField(arcondition, ic, "suNum"));
                data->sv = mxGetData(mxGetField(arcondition, ic, "svNum"));
                
                /* fill inputs */
                fsu(data, tstart, im, isim);
                  
                if(neq>0){
                    /* Load sensitivity initial conditions */
                    sx = N_VCloneVectorArray_Serial(nps, x);
                    if (sx == NULL) {status[0] = 9; return;}
                    for(js=0; js < nps; js++) {
                        sxtmp = NV_DATA_S(sx[js]);
                        for(ks=0; ks < neq; ks++) {
                            sxtmp[ks] = 0.0;
                        }
                    }
                    for (is=0;is<nps;is++) fsx0(is, sx[is], data, im, isim);
                    fsv(data, tstart, x, im, isim);
                    dfxdp(data, tstart, x, returnddxdtdp, im, isim);
                    
                    flag = AR_CVodeSensInit1(cvode_mem, nps, sensi_meth, sensirhs, sx, im, isim);
                    if(flag < 0) {status[0] = 10; return;}
                    
                    /*
                flag = CVodeSensEEtolerances(cvode_mem);
                if(flag < 0) {status[0] = 11; return;}
                     */
                    
                    flag = CVodeSetSensParams(cvode_mem, data->p, NULL, NULL);
                    if (flag < 0) {status[0] = 13; return;}
                    
                    atols_ss = N_VNew_Serial(np);
                    if (atols_ss == NULL) {return;}
                    for (is=0; is<np; is++) Ith(atols_ss, is+1) = cvodes_atol;
                    
                    atolV_ss = N_VCloneVectorArray_Serial(nps, x);
                    if (atolV_ss == NULL) {status[0] = 9; return;}
                    
                    for(js=0; js < nps; js++) {
                        atolV_tmp = NV_DATA_S(atolV_ss[js]);
                        for(ks=0; ks < neq; ks++) {
                            atolV_tmp[ks] = 0.0;
                        }
                    }

                    if(cvodes_atolV_Sens==1)
                    { 
                        for(js=0; js < nps; js++) 
                        {
                            atolV_tmp = NV_DATA_S(atolV_ss[js]);
                            for(ks=0; ks < neq; ks++)
                            {
                                if(y_max_scale_S[ks]==0. || cvodes_atol/y_max_scale_S[ks]>1)
                                {
                                    atolV_tmp[ks] = 1;
                                } else if (cvodes_atol/y_max_scale_S[ks]<1e-8)
                                {   
                                /* && Ith(atolV, ks+1)==1.e-8){*/
                                /*printf("atolVS for neq=%i is %d \n", ks+1, atolV_tmp[ks]);*/
                			    atolV_tmp[ks] = 1e-8;			  
                                }else if(cvodes_atol/y_max_scale_S[ks]>1e-8 && cvodes_atol/y_max_scale_S[ks]<1) 
                                {
                                    atolV_tmp[ks] = cvodes_atol/y_max_scale_S[ks];
                                    /*if(atolV_tmp[ks] < Ith(atolV, ks+1)){
                                         atolV_tmp[ks] = Ith(atolV, ks+1);
                                    }*/
                                }else
                                {
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
                    
                    if(flag < 0) {status[0] = 11; return;}
                    
                    flag = CVodeSetSensErrCon(cvode_mem, error_corr);
                    if(flag < 0) {status[0] = 13; return;}
                }
            }

            /* Do we have a startup event? */
            if ( qEvents == 1 ) {
                if ( event_data->t[event_data->i] == tstart ) {
                    flag = handle_event( cvode_mem, event_data, data, x, sx, nps, neq, sensi, sensi_meth );
                    (event_data->i)++;

                    if (flag < 0) {status[0] = 16; thr_error("Failed to reinitialize solver at event"); return;}
                }
            }
            
            /* loop over output points */
            for (is=0; is < nout; is++) {
                /* printf("%f x-loop (im=%i ic=%i)\n", ts[is], im, ic); */
                
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
                                flag = equilibrate(cvode_mem, data, x, t, equilibrated, returndxdt, teq, neq, im, isim, abortSignal);
                            } else {
                                /* Simulate up to the next time point */
                                flag = CVode(cvode_mem, RCONST(ts[is]), x, &t, CV_NORMAL);
                                data->t = ts[is];
                            }
                            
                            /* Found an event */
                            if ((qEvents==1) && (event_data->i < event_data->n) && (ts[is]==event_data->t[event_data->i])) /*flag==CV_TSTOP_RETURN*/
                            {
                              qEvents = 2;    /* qEvents=2 denotes that an event just happened */
                              flag = 0;     /* Re-set the flag for legacy error-checking reasons */
                            }
                            
                            if ( flag==CV_TSTOP_RETURN )
                            {
                              thr_error( "Error in the event system. Did the model link properly?" );
                            }

                            status[0] = flag;
                        }
                    }
                }

                /* Store time step results */
                if(status[0] == 0.0) {
                    fu(data, ts[is], im, isim);
                    fv(data, ts[is], x, im, isim);
                    
                    for(js=0; js < nu; js++) returnu[js*nout+is] = data->u[js];
                    for(js=0; js < nv; js++) returnv[js*nout+is] = data->v[js];
                    for(js=0; js < neq; js++) {
                        returnx[js*nout+is] = Ith(x, js+1);
                        /* set negative values to zeros */
                        if(qpositivex[js]>0.5 && returnx[js*nout+is]<0.0) returnx[js*nout+is] = 0.0;
                    }
                } else {
                    for(js=0; js < nu; js++) returnu[js*nout+is] = 0.0;
                    for(js=0; js < nv; js++) returnv[js*nout+is] = 0.0;
                    for(js=0; js < neq; js++) returnx[js*nout+is] = 0.0;
                }
                
                /* Store output sensitivities */
                if(status[0] == 0.0) {
                    if (sensi == 1) {
                        if(ts[is] > tstart) {
                            if(neq>0) {
                                flag = CVodeGetSens(cvode_mem, &t, sx);
                                if (flag < 0) {status[0] = 14; return;}
                            }
                        }
                        fsu(data, ts[is], im, isim);
                        fsv(data, ts[is], x, im, isim);
                        
                        for(js=0; js < nps; js++) {
                            if(neq>0) {
                                /* Output state sensitivities */
                                sxtmp = NV_DATA_S(sx[js]);
                                for(ks=0; ks < neq; ks++) {
                                    returnsx[(js*neq+ks)*nout + is] = sxtmp[ks];
                                }
                                
                                /* Output flux sensitivities */
                                csv(ts[is], x, js, sx[js], data, im, ic);
                                for(ks=0; ks < nv; ks++) {
                                    returnsv[(js*nv+ks)*nout + is] = data->sv[ks];
                                }      
                            }
                            
                            /* Output input sensitivities */
                            for(ks=0; ks < nu; ks++) {
                                returnsu[(js*nu+ks)*nout + is] = data->su[(js*nu)+ks];
                            }
                        }
                    }
                } else {
                    if (sensi == 1) {
                        for(js=0; js < nps; js++) {
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
                    flag = handle_event(cvode_mem, event_data, data, x, sx, nps, neq, sensi, sensi_meth );
                    if (flag < 0) {status[0] = 16; thr_error("Failed to reinitialize solver at event"); return;}
                    
                    qEvents = 1;
                    (event_data->i)++;
                }
            } /* End of simulation loop */
            
            /* Free memory */
            if(neq>0) {
                N_VDestroy_Serial(x);
                N_VDestroy_Serial(atolV);
                if (sensi == 1) {
                    N_VDestroyVectorArray_Serial(sx, nps);
                    N_VDestroy_Serial(atols_ss);
                    N_VDestroyVectorArray_Serial(atolV_ss, nps);
                }
               CVodeFree(&cvode_mem);              
            }
            free(data);
            free(event_data);
            
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
            
            status = mxGetData(mxGetField(arcondition, ic, "status"));
            
            /* User data structure */
            data = (UserData) malloc(sizeof *data);
            if (data == NULL) {status[0] = 1; return;}
            
            data->abort = abortSignal;
            data->u = mxGetData(mxGetField(arcondition, ic, "uNum"));
            data->p = mxGetData(mxGetField(arcondition, ic, "pNum"));
            data->v = mxGetData(mxGetField(arcondition, ic, "vNum"));
            
            /* State vectors */
            x = N_VNew_Serial(nx);
            if (x == NULL) {status[0] = 2; return;}
            x_lb = N_VNew_Serial(nx);
            if (x_lb == NULL) {status[0] = 2; return;}
            x_ub = N_VNew_Serial(nx);
            if (x_ub == NULL) {status[0] = 2; return;}
            
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
            
            /* Free memory */
            N_VDestroy_Serial(x);
            N_VDestroy_Serial(x_lb);
            N_VDestroy_Serial(x_ub);
            free(data);
            
            /**** end of SSA ****/
        }
        
        /* call z_calc */
        if ( only_sim == 0 )
            z_calc(im, ic, arcondition, sensi);
    }

    gettimeofday(&t3, NULL);
    
    /* printf("computing model #%i, condition #%i (done)\n", im, ic); */
    
    /* call y_calc */
    if(ardata!=NULL){
        dLink = mxGetField(arcondition, ic, "dLink");
        dLinkints = mxGetData(dLink);
        
        nd = (int) mxGetNumberOfElements(dLink);
        
        /* loop over data */
        for(ids=0; ids<nd; ++ids){
            id = ((int) dLinkints[ids]) - 1;
            has_tExp = (int) mxGetScalar(mxGetField(ardata, id, "has_tExp"));
            
            if((has_tExp == 1) | (fine == 1)) {
                y_calc(im, id, ardata, arcondition, sensi);
            }
        }
    }

    gettimeofday(&t4, NULL);
    timersub(&t2, &t1, &tdiff);
    ticks_start[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
    timersub(&t3, &t1, &tdiff);
    ticks_stop_data[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
    timersub(&t4, &t1, &tdiff);
    ticks_stop[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
    
    /* Signal finished simulation */
    *threadStatus = 1;    
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
    while( !converged )
    {        
        time = time + current_stepsize;

        /* Simulate up to next checkpoint */
        CVodeSetStopTime(cvode_mem, RCONST(time));
        flag = CVode(cvode_mem, RCONST(time), x, &t, CV_NORMAL);

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
                converged = (converged && ( (equilibrated[i] < 0.1) || ( fabs(returndxdt[i])<eq_tol ) ) );
        } else {
            for (i=0; i<neq; i++)
                converged = (converged && fabs(returndxdt[i])<eq_tol );
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
    #ifdef HAS_PTHREAD
        if(parallel==1) {pthread_exit(NULL);}
    #endif
}

/* Event handler */
/* Put functions that are supposed to be evaluated on events here */
int handle_event( void* cvode_mem, EventData event_data, UserData user_data, N_Vector x, N_Vector *sx, int nps, int neq, int sensi, int sensi_meth )
{
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
        }
	}

    /* Reinitialize the solver */
    flag = CVodeReInit(cvode_mem, RCONST(event_data->t[event_data->i]), x);
    if (flag>=0) {
        if (sensi==1) {
            flag = CVodeSensReInit(cvode_mem, sensi_meth, sx);
        }
    }
        
    return flag;
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

/* This function initializes time point lists */
int init_list( mxArray* arcondition, int ic, double tstart, int* nPoints, double** timePoints, int* currentIndex, const char* flagFieldName, const char* timePointFieldName ) {
    int ID, flag;
    double *time;
          
    flag = (int) mxGetScalar(mxGetField(arcondition, ic, flagFieldName));
    if (flag==1) {
        mxArray *timePointField = mxGetField(arcondition, ic, timePointFieldName);

        if ( timePointField != NULL ) {
             time = (double*) mxGetData( timePointField );
             *nPoints     = (int) mxGetNumberOfElements( timePointField );

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
void z_calc(int im, int ic, mxArray *arcondition, int sensi) {
    
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
        
        fz(t[it], nt, it, 0, 0, 0, 0, z, p, u, x, im, ic);
	if( fine == 0 ) {
	  dfzdx(t[it], nt, it, 0, nx, 0, 0, dzdx, z, p, u, x, im, ic);
	}
        if (sensi == 1) {
            fsz(t[it], nt, it, np, sz, p, u, x, z, su, sx, im, ic);
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
    double *res;
    double *reserr;
    
    double *sy;
    double *systd;
    double *sres;
    double *sreserr;
    
    double *qlogy;
    double *qlogp;
    
    double *p;
    double *u;
    double *x;
    double *z;
    double *su;
    double *sx;
    double *sz;
    double *y_scale;
    double *dzdx;
    
    double *chi2;
    double *chi2err;
    
    double fiterrors_correction_factor;
    
    if(useFitErrorMatrix == 1 && fiterrors_matrix[id*nrows_fiterrors_matrix+im] != 1) {
        fiterrors_correction_factor = 1;
    } else {
        fiterrors_correction_factor = fiterrors_correction;
    }
    
    /* MATLAB values */
    ic = (int) mxGetScalar(mxGetField(ardata, id, "cLink")) - 1;
    has_yExp = (int) mxGetScalar(mxGetField(ardata, id, "has_yExp"));
    
    ny = (int) mxGetNumberOfElements(mxGetField(ardata, id, "y"));
    qlogy = mxGetData(mxGetField(ardata, id, "logfitting"));
    qlogp = mxGetData(mxGetField(ardata, id, "qLog10"));
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
            if( (useFitErrorMatrix == 0 && fiterrors == -1) || (useFitErrorMatrix == 1 && fiterrors_matrix[id*nrows_fiterrors_matrix+im] == -1) ) {
                ystd = mxGetData(mxGetField(ardata, id, "yExpStd"));
             }
            res = mxGetData(mxGetField(ardata, id, "res"));
            reserr = mxGetData(mxGetField(ardata, id, "reserr"));
            if (sensi == 1) {
                sres = mxGetData(mxGetField(ardata, id, "sres"));
                sreserr = mxGetData(mxGetField(ardata, id, "sreserr"));
            }
            chi2 = mxGetData(mxGetField(ardata, id, "chi2"));
            chi2err = mxGetData(mxGetField(ardata, id, "chi2err"));
            for(iy=0; iy<ny; iy++) {
                chi2[iy] = 0.0;
                if( (useFitErrorMatrix == 0 && fiterrors == 1) || (useFitErrorMatrix == 1 && fiterrors_matrix[id*nrows_fiterrors_matrix+im] == 1) ) {
                    chi2err[iy] = 0.0;
                }
            }
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
                if(y[it + (iy*nt)]<-cvodes_atol) printf("WARNING, check for concentrations <= 0 !!!\n");
                if(fine==0)  y_scale[it + (iy*nt)] = y_scale[it + (iy*nt)] / y[it + (iy*nt)] / log(10.0);
                y[it + (iy*nt)] = log10(y[it + (iy*nt)]);
            }
        }
        
        if( (useFitErrorMatrix == 0 && fiterrors != -1) || (useFitErrorMatrix == 1 && fiterrors_matrix[id*nrows_fiterrors_matrix+im] != -1) ) {
            fystd(t[it], nt, it, ntlink, itlink, ystd, y, p, u, x, z, im, id);
        }

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
            
            if( (useFitErrorMatrix == 0 && fiterrors != -1) || (useFitErrorMatrix == 1 && fiterrors_matrix[id*nrows_fiterrors_matrix+im] != -1) ) {
                fsystd(t[it], nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz, im, id);
            }
        }
        
        if ((has_yExp == 1) & (fine == 0)) {
            fres(nt, ny, it, res, y, yexp, ystd, chi2, fiterrors_correction_factor);
            if(sensi == 1) fsres(nt, ny, np, it, sres, sy, yexp, ystd, fiterrors_correction_factor);
            
            if( (useFitErrorMatrix == 0 && fiterrors == 1) || (useFitErrorMatrix == 1 && fiterrors_matrix[id*nrows_fiterrors_matrix+im] == 1) ) {
                fres_error(nt, ny, it, reserr, res, y, yexp, ystd, chi2err);
                if (sensi == 1) fsres_error(nt, ny, np, it, sres, sreserr, sy, systd, y, yexp, ystd, res, reserr);
            }
        }
        
        /* log trafo of parameters */
        if ((sensi == 1) & (has_yExp == 1) & (fine == 0)) {
            for (ip=0; ip < np; ip++) {
                if (qlogp[ip] > 0.5) {
                    for (iy=0; iy<ny; iy++) {
                        sres[it + (iy*nt) + (ip*nt*ny)] *= p[ip] * log(10.0);
                        if( (useFitErrorMatrix == 0 && fiterrors == 1) || (useFitErrorMatrix == 1 && fiterrors_matrix[id*nrows_fiterrors_matrix+im] == 1) ) {
                            sreserr[it + (iy*nt) + (ip*nt*ny)] *= p[ip] * log(10.0);
                        }
                    }
                }
            }
        }
    }
    
    /* printf("computing model #%i, data #%i (done)\n", im, id); */
}

/* standard least squares */
void fres(int nt, int ny, int it, double *res, double *y, double *yexp, double *ystd, double *chi2, double fiterrors_correction_factor) {
    int iy;
    
    for(iy=0; iy<ny; iy++){
        res[it + (iy*nt)] = (yexp[it + (iy*nt)] - y[it + (iy*nt)]) / ystd[it + (iy*nt)] * sqrt(fiterrors_correction_factor);
        /* in case of missing data (nan) */
        if(mxIsNaN(yexp[it + (iy*nt)])) {
            res[it + (iy*nt)] = 0.0;
            y[it + (iy*nt)] = yexp[it + (iy*nt)];
            ystd[it + (iy*nt)] = yexp[it + (iy*nt)];
        }
        /* in case of Inf data after log10(0) */
        if(mxIsInf(yexp[it + (iy*nt)])) {
            res[it + (iy*nt)] = 0.0;
        }
        chi2[iy] += pow(res[it + (iy*nt)], 2);
    }
}
void fsres(int nt, int ny, int np, int it, double *sres, double *sy, double *yexp, double *ystd, double fiterrors_correction_factor) {
    int iy, ip;
    
    for(iy=0; iy<ny; iy++){
        for(ip=0; ip<np; ip++){
            sres[it + (iy*nt) + (ip*nt*ny)] = - sy[it + (iy*nt) + (ip*nt*ny)] / ystd[it + (iy*nt)] * sqrt(fiterrors_correction_factor);
            /* in case of missing data (nan) */
            if(mxIsNaN(yexp[it + (iy*nt)])) {
                sres[it + (iy*nt) + (ip*nt*ny)] = 0.0;
            }
            /* in case of Inf data after log10(0) */
            if(mxIsInf(yexp[it + (iy*nt)])) {
                sres[it + (iy*nt) + (ip*nt*ny)] = 0.0;
            }
        }
    }
}

/* least squares for error model fitting */
void fres_error(int nt, int ny, int it, double *reserr, double *res, double *y, double *yexp, double *ystd, double *chi2err) {
    int iy;
    
    double add_c = 50.0;
    
    for(iy=0; iy<ny; iy++){
        reserr[it + (iy*nt)] = 2.0*log(ystd[it + (iy*nt)]);
        if(mxIsNaN(yexp[it + (iy*nt)])) {
            reserr[it + (iy*nt)] = 0.0;
            y[it + (iy*nt)] = yexp[it + (iy*nt)];
            ystd[it + (iy*nt)] = yexp[it + (iy*nt)];
        } else {
            reserr[it + (iy*nt)] += add_c;
            /* 2*log(ystd) + add_c > 0 */
            if(reserr[it + (iy*nt)] < 0) {
                printf("ERROR error model < 1e-10 not allowed\n");
                return;
            }
            reserr[it + (iy*nt)] = sqrt(reserr[it + (iy*nt)]);
            chi2err[iy] += pow(reserr[it + (iy*nt)], 2) - add_c;
        }
    }
}
void fsres_error(int nt, int ny, int np, int it, double *sres, double *sreserr, double *sy, double *systd, double *y, double *yexp, double *ystd, double *res, double *reserr) {
    int iy, ip;
    
    for(iy=0; iy<ny; iy++){
        for(ip=0; ip<np; ip++){
            sres[it + (iy*nt) + (ip*nt*ny)] -= systd[it + (iy*nt) + (ip*nt*ny)] * res[it + (iy*nt)] / ystd[it + (iy*nt)];
            sreserr[it + (iy*nt) + (ip*nt*ny)] = systd[it + (iy*nt) + (ip*nt*ny)] / (reserr[it + (iy*nt)] * ystd[it + (iy*nt)]);
            if(mxIsNaN(yexp[it + (iy*nt)])) {
                sres[it + (iy*nt) + (ip*nt*ny)] = 0.0;
                sreserr[it + (iy*nt) + (ip*nt*ny)] = 0.0;
            }
        }
    }
}


/* custom error weight function */
/*
int ewt(N_Vector y, N_Vector w, void *user_data)
{
  int i;
  realtype yy, ww;
  
  for (i=1; i<=NV_LENGTH_S(y); i++) {
    yy = Ith(y,i);
    ww = cvodes_rtol * ABS(yy) + cvodes_atol;  
    if (ww <= 0.0) return (-1);
    Ith(w,i) = 1.0/ww;
    printf("%e ", ww);
  }
  printf("\n");
  return(0);
} 
*/
