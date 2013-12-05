/*
 *  MATLAB usage: arSimuCalc(struct ar, int fine, int sensi)
 *
 *  (adaptation from Scott D. Cohen, Alan C. Hindmarsh, and Radu Serban @ LLNL)
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

#ifdef HAS_SYSTIME
#include <sys/time.h>
#endif

#include <cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDENSE fcts. and constants */
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

/* Accessor macros */
#define Ith(v, i)     NV_Ith_S(v, i-1)        /* i-th vector component i=1..neq */
#define IJth(A, i, j) DENSE_ELEM(A, i-1, j-1) /* (i,j)-th matrix component i,j=1..neq */

#define MXNCF        20
#define MXNEF        20

#ifdef HAS_PTHREAD
struct thread_data_x {
    int	id;
};
#endif

mxArray *armodel;
mxArray *arthread;

int    fine;
int    sensi;
int    dynamics;
int    jacobian;
int    parallel;
int    fiterrors;
int    cvodes_maxsteps;
double cvodes_rtol;
double cvodes_atol;
double fiterrors_correction;

#ifdef HAS_SYSTIME
struct timeval t1;
#endif

/* Prototypes of private functions */
#ifdef HAS_PTHREAD
void *thread_calc(void *threadarg);
#else
void thread_calc(int id);
#endif
void x_calc(int im, int ic);
void y_calc(int im, int id, mxArray *ardata, mxArray *arcondition);

void fres(int nt, int ny, int it, double *res, double *y, double *yexp, double *ystd, double *chi2);
void fsres(int nt, int ny, int np, int it, double *sres, double *sy, double *yexp, double *ystd);
void fres_error(int nt, int ny, int it, double *reserr, double *res, double *y, double *yexp, double *ystd, double *chi2);
void fsres_error(int nt, int ny, int np, int it, double *sres, double *sreserr, double *sy, double *systd, double *y, double *yexp, double *ystd, double *res, double *reserr);

int ewt(N_Vector y, N_Vector w, void *user_data);

/* user functions */
#include "arSimuCalcFunctions.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int nthreads, ithreads, tid;
#ifdef HAS_PTHREAD
    int rc;
    pthread_t threads_x[NMAXTHREADS];
    struct thread_data_x thread_data_x_array[NMAXTHREADS];
#endif

    mxArray    *arconfig;
    
#ifdef HAS_SYSTIME
    struct timeval t2, tdiff;
    double *ticks_stop;
    ticks_stop = mxGetData(mxGetField(prhs[0], 0, "stop"));
    gettimeofday(&t1, NULL);
#endif
    
    /* get ar.model */
    armodel = mxGetField(prhs[0], 0, "model");
    if(armodel==NULL){
        mexErrMsgTxt("field ar.model not existing");
    }
    
    fine = (int) mxGetScalar(prhs[1]);
    sensi = (int) mxGetScalar(prhs[2]);
    dynamics = (int) mxGetScalar(prhs[3]);
    
    /* get ar.config */
    arconfig = mxGetField(prhs[0], 0, "config");
    parallel = (int) mxGetScalar(mxGetField(arconfig, 0, "useParallel"));
    jacobian = (int) mxGetScalar(mxGetField(arconfig, 0, "useJacobian"));
    cvodes_rtol = mxGetScalar(mxGetField(arconfig, 0, "rtol"));
    cvodes_atol = mxGetScalar(mxGetField(arconfig, 0, "atol"));
    cvodes_maxsteps = (int) mxGetScalar(mxGetField(arconfig, 0, "maxsteps"));
    fiterrors = (int) mxGetScalar(mxGetField(arconfig, 0, "fiterrors"));
    fiterrors_correction = (double) mxGetScalar(mxGetField(arconfig, 0, "fiterrors_correction"));
            
    /* threads */
    arthread = mxGetField(arconfig, 0, "threads");
    nthreads = (int) mxGetScalar(mxGetField(arconfig, 0, "nThreads"));
    
#ifdef HAS_PTHREAD
    if(NMAXTHREADS<nthreads) mexErrMsgTxt("ERROR at NMAXTHREADS < nthreads");
#endif
    
    /* printf("%i threads (%i fine, %i sensi, %i jacobian, %g rtol, %g atol, %i maxsteps)\n", nthreads, fine,
            sensi, jacobian, cvodes_rtol, cvodes_atol, cvodes_maxsteps); */
    
#ifdef HAS_PTHREAD
    /* loop over threads parallel */
    for(ithreads=0; ithreads<nthreads; ++ithreads){
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
    
    /* wait for termination of condition threads */
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
    
#ifdef HAS_SYSTIME
    gettimeofday(&t2, NULL);
    timersub(&t2, &t1, &tdiff);
    ticks_stop[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
#endif
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
        x_calc(ms[in], cs[in]);
    }
    
    /* printf("computing thread #%i(done)\n", id); */
    
#ifdef HAS_PTHREAD
    if(parallel==1) {pthread_exit(NULL);}
#endif
}

/* calculate dynamics by CVODES */
void x_calc(int im, int ic) {
    mxArray    *arcondition;
    mxArray    *ardata;
    
    mxArray *dLink;
    double *dLinkints;      
    
    int nm, nc, id, nd, has_tExp;
    int flag;
    int is, js, ks, ids;
    int nout, neq;
    int nu, np, nps, nv;
    
    void *cvode_mem;
    UserData data;
    
    realtype t;
    double tstart;
    N_Vector x;
    N_Vector atols_ss;
    N_Vector *sx;
    realtype *sxtmp;
    
    double *qpositivex;
    double *status;
    double *ts;
    double *returnu;
    double *returnsu;
    double *returnv;
    double *returnsv;
    double *returnx;
    double *returnsx;
    double *returndxdt;
    double *returnddxdtdp;
    
    int sensi_meth = CV_SIMULTANEOUS; /* CV_SIMULTANEOUS or CV_STAGGERED */
    bool error_corr = TRUE;
        
    /* printf("computing model #%i, condition #%i\n", im, ic); */
            
    /* check if im in range */
    nm = mxGetNumberOfElements(armodel);
    if(nm<=im) {
        printf("im > length(ar.model)\n");
#ifdef HAS_PTHREAD
        if(parallel==1) {pthread_exit(NULL);}
#endif
        return;
    }
        
    /* get ar.model(im).condition */
    arcondition = mxGetField(armodel, im, "condition");
    if(arcondition==NULL){
        printf("field ar.model.condition not existing\n");
#ifdef HAS_PTHREAD
        if(parallel==1) {pthread_exit(NULL);}
#endif
        return;
    }
    
    /* check if ic in range */
    nc = mxGetNumberOfElements(arcondition);
    if(nc<=ic) {
        printf("ic > length(ar.model.condition)\n");
#ifdef HAS_PTHREAD
        if(parallel==1) {pthread_exit(NULL);}
#endif
        return;
    }
    
    has_tExp = (int) mxGetScalar(mxGetField(arcondition, ic, "has_tExp"));
    if(has_tExp == 0 && fine == 0) return;
    
    /* get ar.model(im).data */
    ardata = mxGetField(armodel, im, "data");
     
#ifdef HAS_SYSTIME
    struct timeval t2, t3, t4, tdiff;
    double *ticks_start, *ticks_stop_data, *ticks_stop;
    ticks_start = mxGetData(mxGetField(arcondition, ic, "start"));
    ticks_stop = mxGetData(mxGetField(arcondition, ic, "stop"));
    ticks_stop_data = mxGetData(mxGetField(arcondition, ic, "stop_data"));

    gettimeofday(&t2, NULL);
#endif
    
    if(dynamics == 1) { /* begin of CVODES */
        
        /* get MATLAB values */
        qpositivex = mxGetData(mxGetField(armodel, im, "qPositiveX"));
        status = mxGetData(mxGetField(arcondition, ic, "status"));
        tstart = mxGetScalar(mxGetField(arcondition, ic, "tstart"));
        neq = mxGetNumberOfElements(mxGetField(armodel, im, "xs"));
        
        if(fine == 1){
            ts = mxGetData(mxGetField(arcondition, ic, "tFine"));
            nout = mxGetNumberOfElements(mxGetField(arcondition, ic, "tFine"));
            
            returnu = mxGetData(mxGetField(arcondition, ic, "uFineSimu"));
            returnv = mxGetData(mxGetField(arcondition, ic, "vFineSimu"));
            returnx = mxGetData(mxGetField(arcondition, ic, "xFineSimu"));
            if (sensi == 1) {
                returnsu = mxGetData(mxGetField(arcondition, ic, "suFineSimu"));
                returnsv = mxGetData(mxGetField(arcondition, ic, "svFineSimu"));
                returnsx = mxGetData(mxGetField(arcondition, ic, "sxFineSimu"));
            }
        }
        else{
            ts = mxGetData(mxGetField(arcondition, ic, "tExp"));
            nout = mxGetNumberOfElements(mxGetField(arcondition, ic, "tExp"));
            
            returnu = mxGetData(mxGetField(arcondition, ic, "uExpSimu"));
            returnv = mxGetData(mxGetField(arcondition, ic, "vExpSimu"));
            returnx = mxGetData(mxGetField(arcondition, ic, "xExpSimu"));
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
        
        data->qpositivex = qpositivex;
        data->u = mxGetData(mxGetField(arcondition, ic, "uNum"));
        nu = mxGetNumberOfElements(mxGetField(arcondition, ic, "uNum"));
        
        data->p = mxGetData(mxGetField(arcondition, ic, "pNum"));
        np = mxGetNumberOfElements(mxGetField(arcondition, ic, "pNum"));
        nps = np;
        
        data->v = mxGetData(mxGetField(arcondition, ic, "vNum"));
        nv = mxGetNumberOfElements(mxGetField(arcondition, ic, "vNum"));
        data->dvdx = mxGetData(mxGetField(arcondition, ic, "dvdxNum"));
        data->dvdu = mxGetData(mxGetField(arcondition, ic, "dvduNum"));
        data->dvdp = mxGetData(mxGetField(arcondition, ic, "dvdpNum"));
        
        /* fill for t=0 */
        fu(data, tstart, im, ic);
        
        if(neq>0){
            /* Initial conditions */
            x = N_VNew_Serial(neq);
            if (x == NULL) {status[0] = 2; return;}
            for (is=0; is<neq; is++) Ith(x, is+1) = 0.0;
            fx0(x, data, im, ic);
            fv(data, tstart, x, im, ic);
            fx(tstart, x, returndxdt, data, im, ic);
            
            /* Create CVODES object */
            cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
            if (cvode_mem == NULL) {status[0] = 3; return;}
            
            /* Allocate space for CVODES */
            flag = AR_CVodeInit(cvode_mem, x, tstart, im, ic);
            if (flag < 0) {status[0] = 4; return;}
            
            /* Number of maximal internal steps */
            flag = CVodeSetMaxNumSteps(cvode_mem, cvodes_maxsteps);
            if(flag < 0) {status[0] = 15; return;}
            
            /* Use private function to compute error weights */
            flag = CVodeSStolerances(cvode_mem, RCONST(cvodes_rtol), RCONST(cvodes_atol));
            if (flag < 0) {status[0] = 5; return;}
            
            /* Attach user data */
            flag = CVodeSetUserData(cvode_mem, data);
            if (flag < 0) {status[0] = 6; return;}
            
            /* Attach linear solver */
            flag = CVDense(cvode_mem, neq);
            if (flag < 0) {status[0] = 7; return;}
            
            /* Jacobian-related settings */
            if (jacobian == 1) {
                flag = AR_CVDlsSetDenseJacFn(cvode_mem, im, ic);
                if (flag < 0) {status[0] = 8; return;}
            }
            
            /* custom error weight function */
            /*
        flag = CVodeWFtolerances(cvode_mem, ewt);
        if (flag < 0) return;}
             */
        }
        
        /* Sensitivity-related settings */
        if (sensi == 1) {
            /* User data structure */
            data->su = mxGetData(mxGetField(arcondition, ic, "suNum"));
            data->sv = mxGetData(mxGetField(arcondition, ic, "svNum"));
            
            /* fill inputs */
            fsu(data, tstart, im, ic);
            
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
                for (is=0;is<nps;is++) fsx0(is, sx[is], data, im, ic);
                fsv(data, tstart, x, im, ic);
                dfxdp(data, tstart, x, returnddxdtdp, im, ic);
                
                flag = AR_CVodeSensInit1(cvode_mem, nps, sensi_meth, sx, im, ic);
                if(flag < 0) {status[0] = 10; return;}
                
                /*
            flag = CVodeSensEEtolerances(cvode_mem);
            if(flag < 0) {status[0] = 11; return;}
             
            flag = CVodeSetSensParams(cvode_mem, data->p, NULL, NULL);
            if (flag < 0) {status[0] = 13; return;}
                 */
                
                atols_ss = N_VNew_Serial(np);
                if (atols_ss == NULL) {return;}
                for (is=0; is<np; is++) Ith(atols_ss, is+1) = cvodes_atol;
                
                flag = CVodeSensSStolerances(cvode_mem, RCONST(cvodes_rtol), N_VGetArrayPointer(atols_ss));
                if(flag < 0) {status[0] = 11; return;}
                
                flag = CVodeSetSensErrCon(cvode_mem, error_corr);
                if(flag < 0) {status[0] = 13; return;}
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
                        flag = CVode(cvode_mem, RCONST(ts[is]), x, &t, CV_NORMAL);
                        status[0] = flag;
                        /*
                    if(flag==-1) printf("CVODES stoped at t=%f, TOO_MUCH_WORK, did not reach output time after %i steps (m=%i, c=%i).\n", t, cvodes_maxsteps, im, ic);
                    if(flag<-1) printf("CVODES stoped at t=%f (m=%i, c=%i).\n", t, im, ic);
                         */
                    }
                }
                fu(data, ts[is], im, ic);
                fv(data, ts[is], x, im, ic);
                
                for(js=0; js < nu; js++) returnu[js*nout+is] = data->u[js];
                for(js=0; js < nv; js++) {
                    returnv[js*nout+is] = data->v[js];
                }
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
            
            /* only set output sensitivities if no errors occured */
            if(status[0] == 0.0) {
                if (sensi == 1) {
                    if(ts[is] > tstart) {
                        if(neq>0) {
                            flag = CVodeGetSens(cvode_mem, &t, sx);
                            if (flag < 0) {status[0] = 14; return;}
                        }
                    }
                    fsu(data, ts[is], im, ic);
                    fsv(data, ts[is], x, im, ic);
                    
                    for(js=0; js < nps; js++) {
                        if(neq>0) {
                            sxtmp = NV_DATA_S(sx[js]);
                            for(ks=0; ks < neq; ks++) {
                                returnsx[js*neq*nout + ks*nout + is] = sxtmp[ks];
                            }
                        }
                        for(ks=0; ks < nu; ks++) {
                            returnsu[js*nu*nout + ks*nout + is] = data->su[(js*nu)+ks];
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
                    }
                }
            }
        }
        
        /* Free memory */
        if(neq>0) {
            N_VDestroy_Serial(x);
            if (sensi == 1) {
                N_VDestroyVectorArray_Serial(sx, nps);
            }
            CVodeFree(&cvode_mem);
        }
        free(data);
        
    } /* end of CVODES */

#ifdef HAS_SYSTIME
    gettimeofday(&t3, NULL);
#endif
    
    /* printf("computing model #%i, condition #%i (done)\n", im, ic); */
    
    /* call y_calc */
    if(ardata!=NULL){
        dLink = mxGetField(arcondition, ic, "dLink");
        dLinkints = mxGetData(dLink);
        
        nd = mxGetNumberOfElements(dLink);
        
        /* loop over data */
        for(ids=0; ids<nd; ++ids){
            id = ((int) dLinkints[ids]) - 1;
            has_tExp = (int) mxGetScalar(mxGetField(ardata, id, "has_tExp"));
            
            if(has_tExp == 1 | fine == 1) {
                y_calc(im, id, ardata, arcondition);
            }
        }
    }
    
#ifdef HAS_SYSTIME
    gettimeofday(&t4, NULL);
    timersub(&t2, &t1, &tdiff);
    ticks_start[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
    timersub(&t3, &t1, &tdiff);
    ticks_stop_data[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
    timersub(&t4, &t1, &tdiff);
    ticks_stop[0] = ((double) tdiff.tv_usec) + ((double) tdiff.tv_sec * 1e6);
#endif
}



/* calculate observations */
void y_calc(int im, int id, mxArray *ardata, mxArray *arcondition) {
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
    double *su;
    double *sx;
    
    double *chi2;
    double *chi2err;
    
    /* MATLAB values */
    ic = (int) mxGetScalar(mxGetField(ardata, id, "cLink")) - 1;
    has_yExp = (int) mxGetScalar(mxGetField(ardata, id, "has_yExp"));
    
    ny = mxGetNumberOfElements(mxGetField(ardata, id, "y"));
    qlogy = mxGetData(mxGetField(ardata, id, "logfitting"));
    qlogp = mxGetData(mxGetField(ardata, id, "qLog10"));
    p = mxGetData(mxGetField(ardata, id, "pNum"));
    np = mxGetNumberOfElements(mxGetField(ardata, id, "pNum"));
    
    if(fine == 1){
        t = mxGetData(mxGetField(ardata, id, "tFine"));
        nt = mxGetNumberOfElements(mxGetField(ardata, id, "tFine"));
        tlink = mxGetData(mxGetField(ardata, id, "tLinkFine"));
        ntlink = mxGetNumberOfElements(mxGetField(arcondition, ic, "tFine"));
        
        y = mxGetData(mxGetField(ardata, id, "yFineSimu"));
        ystd = mxGetData(mxGetField(ardata, id, "ystdFineSimu"));
        
        u = mxGetData(mxGetField(arcondition, ic, "uFineSimu"));
        x = mxGetData(mxGetField(arcondition, ic, "xFineSimu"));
        
        if (sensi == 1) {
            sy = mxGetData(mxGetField(ardata, id, "syFineSimu"));
            systd = mxGetData(mxGetField(ardata, id, "systdFineSimu"));
            
            su = mxGetData(mxGetField(arcondition, ic, "suFineSimu"));
            sx = mxGetData(mxGetField(arcondition, ic, "sxFineSimu"));
        }
    }
    else{
        t = mxGetData(mxGetField(ardata, id, "tExp"));
        nt = mxGetNumberOfElements(mxGetField(ardata, id, "tExp"));
        tlink = mxGetData(mxGetField(ardata, id, "tLinkExp"));
        ntlink = mxGetNumberOfElements(mxGetField(arcondition, ic, "tExp"));
        
        y = mxGetData(mxGetField(ardata, id, "yExpSimu"));
        ystd = mxGetData(mxGetField(ardata, id, "ystdExpSimu"));
        
        u = mxGetData(mxGetField(arcondition, ic, "uExpSimu"));
        x = mxGetData(mxGetField(arcondition, ic, "xExpSimu"));
        
        if (sensi == 1) {
            sy = mxGetData(mxGetField(ardata, id, "syExpSimu"));
            systd = mxGetData(mxGetField(ardata, id, "systdExpSimu"));
            
            su = mxGetData(mxGetField(arcondition, ic, "suExpSimu"));
            sx = mxGetData(mxGetField(arcondition, ic, "sxExpSimu"));
        }
        
        if (has_yExp == 1) {
            yexp = mxGetData(mxGetField(ardata, id, "yExp"));
            if(fiterrors==-1) ystd = mxGetData(mxGetField(ardata, id, "yExpStd"));
            
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
                if(fiterrors==1) chi2err[iy] = 0.0;
            }
        }
    }
    
    /* loop over output points */
    for (it=0; it < nt; it++) {
        /* printf("%f y-loop (im=%i id=%i)\n", t[it], im, id); */
        itlink = (int) tlink[it] - 1;
        
        fy(t[it], nt, it, ntlink, itlink, 0, 0, 0, y, p, u, x, im, id);
        
        /* log trafo of y */
        for (iy=0; iy<ny; iy++) {
            if(qlogy[iy] > 0.5){
                if(y[it + (iy*nt)]<0.0) printf("WARNING, check for concentrations <= 0 !!!\n");
                y[it + (iy*nt)] = log10(y[it + (iy*nt)]);
            }
        }
        
        if(fiterrors!=-1) fystd(t[it], nt, it, ntlink, itlink, ystd, y, p, u, x, im, id);
        
        if (sensi == 1) {
            fsy(t[it], nt, it, ntlink, itlink, sy, p, u, x, su, sx, im, id);
            
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
            
            if(fiterrors!=-1) fsystd(t[it], nt, it, ntlink, itlink, systd, p, y, u, x, sy, su, sx, im, id);
        }
        
        if ((has_yExp == 1) & (fine == 0)) {
            fres(nt, ny, it, res, y, yexp, ystd, chi2);
            if(sensi == 1) fsres(nt, ny, np, it, sres, sy, yexp, ystd);
            
            if(fiterrors==1) {
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
                        if(fiterrors==1) sreserr[it + (iy*nt) + (ip*nt*ny)] *= p[ip] * log(10.0);
                    }
                }
            }
        }
    }
    
    /* printf("computing model #%i, data #%i (done)\n", im, id); */
}

/* standard least squares */
void fres(int nt, int ny, int it, double *res, double *y, double *yexp, double *ystd, double *chi2) {
    int iy;
    
    for(iy=0; iy<ny; iy++){
        res[it + (iy*nt)] = (yexp[it + (iy*nt)] - y[it + (iy*nt)]) / ystd[it + (iy*nt)] * sqrt(fiterrors_correction);
        if(mxIsNaN(yexp[it + (iy*nt)])) {
            res[it + (iy*nt)] = 0.0;
            y[it + (iy*nt)] = yexp[it + (iy*nt)];
            ystd[it + (iy*nt)] = yexp[it + (iy*nt)];
        }
        chi2[iy] += pow(res[it + (iy*nt)], 2);
    }
}
void fsres(int nt, int ny, int np, int it, double *sres, double *sy, double *yexp, double *ystd) {
    int iy, ip;
    
    for(iy=0; iy<ny; iy++){
        for(ip=0; ip<np; ip++){
            sres[it + (iy*nt) + (ip*nt*ny)] = - sy[it + (iy*nt) + (ip*nt*ny)] / ystd[it + (iy*nt)] * sqrt(fiterrors_correction);
            if(mxIsNaN(yexp[it + (iy*nt)])) {
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