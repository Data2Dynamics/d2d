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

#endif /* _MY_UDATA */
