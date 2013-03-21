/*
 *  MATLAB usage: arSSACalc(struct ar)
 *
 *  Copyright Andreas Raue 2011 (andreas.raue@fdm.uni-freiburg.de)
 *
 * (a) generate two randon numbers r1 and r2 uniformly distributed in (0,1)
 * (b) Compute alpha0 = \sum_{i=1}^q alpha_i(t),
 *     where alpha_i(t) = propensity function of i-th reaction
 * (c) The next reaction take place at time t+tau, where
 *     tau = 1/alpha0 * log(1/r1)   ... natural logrithm
 * (d) Determine which reaction occurs. Find j such that
 *     r2 >= 1/alpha0 \sum_{i=1}^{j-1} alpha_i(t)
 *     r2 <  1/alpha0 \sum_{i=1}^j     alpha_i(t)
 *     Update the number of reactants and products of the j-th reaction
 * (e) Go to (a) with t = t + tau
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <pthread.h>

#include <cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDENSE fcts. and constants */
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

/* Accessor macros */
#define Ith(v, i)    NV_Ith_S(v, i-1)       /* i-th vector component i=1..neq */
#define IJth(A, i, j) DENSE_ELEM(A, i-1, j-1) /* (i,j)-th matrix component i,j=1..neq */

/* user variables */
#include "arSimuCalcVariables.c"

struct thread_data_x {
    int	im;
    int ic;
    mxArray *arcondition;
};

struct thread_data_y {
    int	im;
    int id;
    mxArray *ardata;
    mxArray *arcondition;
    int dt;
};

typedef struct {
    double *u;
    double *su;
    double *p;
    double *v;
    double *dvdx;
    double *dvdu;
    double *dvdp;
    double *sv;
} *UserData;

mxArray *armodel;

int     parallel;
double  mintau;
int     nruns;

/* Prototypes of private functions */
void *x_calc(void *threadarg);
void *y_calc(void *threadarg);

static int check_flag(void *flagvalue, char *funcname, int opt);

/* user functions */
#include "arSimuCalcFunctions.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    srand( (unsigned)time(NULL) );
    
    int nthreads_x, nthreads_y, nm, im, nc, ic, tid, dtid, nd, id, rc, has_tExp;
    
    mxArray    *arconfig;
    mxArray    *arcondition;
    mxArray    *ardata;
    
    /* get ar.model */
    armodel = mxGetField(prhs[0], 0, "model");
    if(armodel==NULL){
        mexErrMsgTxt("field ar.model not existing");
    }
    
    parallel = 1;
    
    /* get ar.config */
    arconfig = mxGetField(prhs[0], 0, "config");
    mintau = mxGetScalar(mxGetField(arconfig, 0, "ssa_min_tau"));
    nruns = (int) mxGetScalar(mxGetField(arconfig, 0, "ssa_runs"));
    
    /* threads */
    nthreads_x = mxGetScalar(mxGetField(arconfig, 0, "nthreads_x"));
    nthreads_y = mxGetScalar(mxGetField(arconfig, 0, "nthreads_y"));
    pthread_t threads_x[nthreads_x];
    pthread_t threads_y[nthreads_y];
    struct thread_data_x thread_data_x_array[nthreads_x];
    struct thread_data_y thread_data_y_array[nthreads_y];
    
/*    printf("%i x-threads (%i runs, %g min-tau)\n", nthreads_x, nruns, mintau);
    printf("%i y-threads\n", nthreads_y); */
    
    nm = mxGetNumberOfElements(armodel);
    /* loop over models */
    for(im=0; im<nm; ++im){
        
        /* get ar.model(im).condition */
        arcondition = mxGetField(armodel, im, "condition");
        if(arcondition==NULL){
            mexErrMsgTxt("field ar.model.condition not existing");
        }
        
        nc = mxGetNumberOfElements(arcondition);
        /* loop over conditions */
        for(ic=0; ic<nc; ++ic){
            tid = mxGetScalar(mxGetField(arcondition, ic, "thread_id"));
            
            thread_data_x_array[tid].im = im;
            thread_data_x_array[tid].ic = ic;
            thread_data_x_array[tid].arcondition = arcondition;
            
            if(parallel==1){
                /* printf("creating condition thread %i, m=%i, c=%i\n", tid, im, ic); */
                rc = pthread_create(&threads_x[tid], NULL, x_calc, (void *) &thread_data_x_array[tid]);
                if (rc){
                    mexErrMsgTxt("ERROR at pthread_create x");
                }
            } else {
                x_calc(&thread_data_x_array[tid]);
            }
        }
    }
    
    /* wait for termination of condition threads = pthread_exit(NULL);*/
    if(parallel==1){
        /* loop over models */
        for(im=0; im<nm; ++im){
            
            /* get ar.model(im).condition */
            arcondition = mxGetField(armodel, im, "condition");
            if(arcondition==NULL){
                mexErrMsgTxt("field ar.model.condition not existing");
            }
            
            nc = mxGetNumberOfElements(arcondition);
            /* loop over conditions */
            for(ic=0; ic<nc; ++ic){    
                tid = mxGetScalar(mxGetField(arcondition, ic, "thread_id"));
                
                rc = pthread_join(threads_x[tid], NULL);
                if (rc){
                    mexErrMsgTxt("ERROR at pthread_join x");
                }
            }
        }
    }
    
    /* loop over models */
    for(im=0; im<nm; ++im){
        /* get ar.model(im).data */
        ardata = mxGetField(armodel, im, "data");
        
        if(ardata!=NULL){
            arcondition = mxGetField(armodel, im, "condition");
            
            nd = mxGetNumberOfElements(ardata);
            /* loop over data */
            for(id=0; id<nd; ++id){
                tid = mxGetScalar(mxGetField(ardata, id, "thread_id"));
                
                thread_data_y_array[tid].im = im;
                thread_data_y_array[tid].id = id;
                thread_data_y_array[tid].ardata = ardata;
                thread_data_y_array[tid].arcondition = arcondition;
                
                dtid = mxGetScalar(mxGetField(ardata, id, "dthread_id"));
                thread_data_y_array[tid].dt = dtid;
                
                if(parallel==1){
                    /* printf("creating data thread %i, m=%i, d=%i (wait call for condition #%i)\n", tid, im, id, dtid); */
                    rc = pthread_create(&threads_y[tid], NULL, y_calc, (void *) &thread_data_y_array[tid]);
                    if (rc){
                        mexErrMsgTxt("ERROR at pthread_create y");
                    }
                } else {
                    y_calc(&thread_data_y_array[tid]);
                }
            }
        }
    }
    
    /* wait for termination of data threads = pthread_exit(NULL);*/
    if(parallel==1){
        /* loop over models */
        for(im=0; im<nm; ++im){
            /* get ar.model(im).data */
            ardata = mxGetField(armodel, im, "data");
            
            if(ardata!=NULL){
                arcondition = mxGetField(armodel, im, "condition");
                
                nd = mxGetNumberOfElements(ardata);
                /* loop over data */
                for(id=0; id<nd; ++id){
                    tid = mxGetScalar(mxGetField(ardata, id, "thread_id"));
                    
                    rc = pthread_join(threads_y[tid], NULL);
                    if (rc){
                        mexErrMsgTxt("ERROR at pthread_join y");
                    }
                }
            }
            
        }
    }
}



/* calculate stochastic dynamics by SSA */
void *x_calc(void *threadarg) {
    struct thread_data_x *my_data = (struct thread_data_x *) threadarg;
    
    int im = my_data->im;
    int ic = my_data->ic;
    mxArray *arcondition = my_data->arcondition;
    
    /* printf("computing model #%i, condition #%i\n", im, ic); */
    
    /* begin of SSA */
    
    UserData data;
    N_Vector x;
    N_Vector x_lb;
    N_Vector x_ub;
    
    double t, tfin, tau, meantau;
    double r1, r2;
    double alpha0, sumalpha;
    int iruns, it, itexp, ix, iv;
            
    int has_texp = (int) mxGetScalar(mxGetField(arcondition, ic, "has_tExp"));
    double lasttau[] = {1,1,1,1,1,1,1,1,1,1};
    int ilasttau = 0;
            
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
    int nx = mxGetNumberOfElements(mxGetField(armodel, im, "xs"));
    int nt = mxGetNumberOfElements(mxGetField(arcondition, ic, "tFine"));
    int nv = mxGetNumberOfElements(mxGetField(arcondition, ic, "vNum"));
    
    double *texp;
    int ntexp;
    double *xssaexp;
    if(has_texp==1) {
        texp = mxGetData(mxGetField(arcondition, ic, "tExp"));
        ntexp = mxGetNumberOfElements(mxGetField(arcondition, ic, "tExp"));
        xssaexp = mxGetData(mxGetField(arcondition, ic, "xExpSSA"));
    }
    
    /* User data structure */
    data = (UserData) malloc(sizeof *data);
    if (check_flag((void *)data, "malloc", 2)) {return;}
    
    data->u = mxGetData(mxGetField(arcondition, ic, "uNum"));
    data->p = mxGetData(mxGetField(arcondition, ic, "pNum"));
    data->v = mxGetData(mxGetField(arcondition, ic, "vNum"));
    
    /* State vectors */
    x = N_VNew_Serial(nx);
    if (check_flag((void *)x, "N_VNew_Serial", 0)) {return;}
    x_lb = N_VNew_Serial(nx);
    if (check_flag((void *)x, "N_VNew_Serial", 0)) {return;}
    x_ub = N_VNew_Serial(nx);
    if (check_flag((void *)x, "N_VNew_Serial", 0)) {return;}
    
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
        while(t<=tfin & it<nt){
            
            /* (a) */
            r1 = ((double)rand()/(double)RAND_MAX);
            r2 = ((double)rand()/(double)RAND_MAX);
            
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
            while(tfine[it]<t+tau & it<nt){
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
            if(has_texp==1 & ntexp>0) {
                while(itexp<ntexp & texp[itexp]<t+tau){
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
    }
//     printf("model #%i, condition #%i: done\n", im+1, ic+1);
    
    N_VDestroy_Serial(x);
    N_VDestroy_Serial(x_lb);
    N_VDestroy_Serial(x_ub);
    
    free(data);
    
    /* end of SSA */
    
    /* printf("computing model #%i, condition #%i (done)\n", im, ic); */
    
    if(parallel==1) {pthread_exit(NULL);}
}



/* calculate observations */
void *y_calc(void *threadarg) {
    struct thread_data_y *my_data = (struct thread_data_y *) threadarg;
    
    int im = my_data->im;
    int id = my_data->id;
    int dt = my_data->dt;
    int rc;   
    mxArray *ardata = my_data->ardata;
    mxArray *arcondition = my_data->arcondition;
    
    /* wait for dependend x thread */
/*    printf("computing model #%i, data #%i (waiting for %i)\n", im, id, dt);
    if(parallel==1) {
        rc = pthread_join(threads_x[dt], NULL);
        if (rc){
            mexErrMsgTxt("ERROR at pthread_join");
        }
    }*/
    
    /* printf("computing model #%i, data #%i\n", im, id); */
    
    int ic, ny, nx, iy, iruns, has_texp;
    int nt, it, ntlink, itlink;
    int ntexp, itexp, ntexplink, itexplink;
    
    double *qlogy;
    double *p;
    
    double *t;
    double *tlink;
    
    double *y;
    double *y_lb;
    double *y_ub;
    
    double *u;
    
    double *x;
    double *x_lb;
    double *x_ub;
    
    double *texp;
    double *texplink;
    
    double *yexp;
    double *uexp;
    double *xexp;
    
    /* MATLAB values */
    ic = (int) mxGetScalar(mxGetField(ardata, id, "cLink")) - 1;
    has_texp = (int) mxGetScalar(mxGetField(ardata, id, "has_tExp"));
    
    nx = mxGetNumberOfElements(mxGetField(armodel, im, "x"));
    ny = mxGetNumberOfElements(mxGetField(ardata, id, "y"));
    qlogy = mxGetData(mxGetField(ardata, id, "logfitting"));
    p = mxGetData(mxGetField(ardata, id, "pNum"));
    
    t = mxGetData(mxGetField(ardata, id, "tFine"));
    nt = mxGetNumberOfElements(mxGetField(ardata, id, "tFine"));
    tlink = mxGetData(mxGetField(ardata, id, "tLinkFine"));
    ntlink = mxGetNumberOfElements(mxGetField(arcondition, ic, "tFine"));
    
    y = mxGetData(mxGetField(ardata, id, "yFineSSA"));
    y_lb = mxGetData(mxGetField(ardata, id, "yFineSSA_lb"));
    y_ub = mxGetData(mxGetField(ardata, id, "yFineSSA_ub"));
    
    u = mxGetData(mxGetField(arcondition, ic, "uFineSimu"));
    
    x = mxGetData(mxGetField(arcondition, ic, "xFineSSA"));
    x_lb = mxGetData(mxGetField(arcondition, ic, "xFineSSA_lb"));
    x_ub = mxGetData(mxGetField(arcondition, ic, "xFineSSA_ub"));
    
    if(has_texp == 1) {
        texp = mxGetData(mxGetField(ardata, id, "tExp"));
        ntexp = mxGetNumberOfElements(mxGetField(ardata, id, "tExp"));
        texplink = mxGetData(mxGetField(ardata, id, "tLinkExp"));
        ntexplink = mxGetNumberOfElements(mxGetField(arcondition, ic, "tExp"));
        
        yexp = mxGetData(mxGetField(ardata, id, "yExpSSA"));
        uexp = mxGetData(mxGetField(arcondition, ic, "uExpSimu"));
        xexp = mxGetData(mxGetField(arcondition, ic, "xExpSSA"));
    }
    
    /* nruns loop */
    for (iruns=0; iruns<nruns; iruns++) {
        
        /* loop over fine time points */
        for (it=0; it < nt; it++) {
            /*         printf("%f y-loop (im=%i id=%i)\n", t[it], im, id);*/
            itlink = (int) tlink[it] - 1;
            
            fy(t[it], nt, it, ntlink, itlink, ny, nx, iruns, y, p, u, x, im, id);
            fy(t[it], nt, it, ntlink, itlink, ny, nx, iruns, y_lb, p, u, x_lb, im, id);
            fy(t[it], nt, it, ntlink, itlink, ny, nx, iruns, y_ub, p, u, x_ub, im, id);
            
            /* log trafo of y */
            for (iy=0; iy<ny; iy++) {
                if(qlogy[iy] > 0.5){
                    y[it + (iy*nt) + (iruns*nt*ny)] = log10(y[it + (iy*nt) + (iruns*nt*ny)]);
                    y_lb[it + (iy*nt) + (iruns*nt*ny)] = log10(y_lb[it + (iy*nt) + (iruns*nt*ny)]);
                    y_ub[it + (iy*nt) + (iruns*nt*ny)] = log10(y_ub[it + (iy*nt) + (iruns*nt*ny)]);
                }
            }
        }
        
        /* loop over experimental time points */
        if(has_texp == 1) {
            for (itexp=0; itexp < ntexp; itexp++) {
                /*         printf("%f y-loop (im=%i id=%i)\n", texp[itexp], im, id);*/
                itexplink = (int) texplink[itexp] - 1;
                
                fy(texp[itexp], ntexp, itexp, ntexplink, itexplink, ny, nx, iruns, yexp, p, u, xexp, im, id);
                
                /* log trafo of y */
                for (iy=0; iy<ny; iy++) {
                    if(qlogy[iy] > 0.5){
                        yexp[itexp + (iy*ntexp) + (iruns*ntexp*ny)] = log10(yexp[itexp + (iy*ntexp) + (iruns*ntexp*ny)]);
                    }
                }
            }
        }
    }

    /* printf("computing model #%i, data #%i (done)\n", im, id); */
    
    if(parallel==1) pthread_exit(NULL);
}



/*
 * Check function return value of CVODES.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a flag so check if
 *             flag >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer
 */

static int check_flag(void *flagvalue, char *funcname, int opt) {
    int *errflag;
    
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
        /* printf("\nSUNDIALS ERROR: %s() failed - returned NULL pointer\n\n", funcname); */
        return(1);
    }
    
    /* Check if flag < 0 */
    else if (opt == 1) {
        errflag = (int *) flagvalue;
        if (*errflag < 0) {
            /* printf("\nSUNDIALS ERROR: %s() failed with flag = %d\n\n", funcname, *errflag); */
            return(1);
        }
    }
    
    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
        /* printf("\nMEMORY ERROR: %s() failed - returned NULL pointer\n\n", funcname); */
        return(1);
    }
    
    return(0);
}

