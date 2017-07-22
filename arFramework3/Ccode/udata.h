#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */

#ifndef _MY_UDATA
#define _MY_UDATA

typedef struct {
	double *qpositivex;
	double *u;
	double *su;
	double *p;
	double *v;
	double *sv;
	double *dvdx;
	double *dvdu;
	double *dvdp;
    int    nsplines;
    double **splines;
	double t;
    int    *abort;
	} *UserData;

    
typedef struct {
   double* t;				/* Time */
   int     i;				/* Current index/iterator */
   int     n;				/* Length */

   /* Used for explicit changes of state in the event handler */
   /* Can be used for reassignments of the form aX+b          */
   int     overrides;
   double* value_A;
   double* value_B;
   double* sensValue_A;
   double* sensValue_B;

   /* Multiple shooting (currently not implemented) */
   double* tMS;
   int     nMS;
   int     iMS;
   
   } *EventData;

/* Global memory structure */
typedef struct {
    /* State vector */
    N_Vector    x;
    
	/* ODE integration via CVODES */
	N_Vector    atolV;
	N_Vector    *atolV_ss;
	N_Vector    atols_ss;
	N_Vector    *sx;
	EventData   event_data;
	void        *cvode_mem;
	UserData    data;
    
    int         sensi;
    int         neq;
    int         np;
    
	/* SSA integration */
	N_Vector    x_lb;
	N_Vector    x_ub;    
    
    /* Logging purposes */
    int         *threadStatus;
    double      *status;
} *SimMemory;

/* Functions defined in udata.h */
SimMemory simCreate( int *threadStatus, double* status );
void simFree( SimMemory sim_mem );

#endif /* _MY_UDATA */
