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
	} *UserData;
#endif /* _MY_UDATA */
