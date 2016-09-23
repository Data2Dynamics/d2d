#include <stdio.h>
#include <stdlib.h>
#include "udata.h"

SimMemory simCreate( int *threadStatus, double* status )
{
    /* If object creation fails, abort! */
	SimMemory sim_mem = (SimMemory) malloc(sizeof *sim_mem);
	if ( sim_mem == NULL )
		return NULL;

	/* Initialize all to NULL to facilitate easier cleanup */
    
    /* CVODES */
	sim_mem->cvode_mem   	= NULL;
    sim_mem->data           = NULL;
	sim_mem->event_data  	= NULL;
	sim_mem->x           	= NULL;
	sim_mem->atolV       	= NULL;
    sim_mem->sx          	= NULL;
	sim_mem->atols_ss    	= NULL;
	sim_mem->atolV_ss    	= NULL;
    
    /* SSA */
    sim_mem->x_lb           = NULL;
    sim_mem->x_ub           = NULL;
    
    sim_mem->neq            = 0;
    sim_mem->np             = 0;
    sim_mem->sensi          = 0;

	/* Used for storing errors */
	sim_mem->status 		= status;
    
    /* Used for terminating the thread */
	sim_mem->threadStatus   = threadStatus;
}

void simFree( SimMemory sim_mem )
{
    /* Something is seriously wrong. Anything we free would cause a segfault */
	if ( sim_mem == NULL )
		return;
    
	/* CVODES memory */
	if ( sim_mem->data )
		free(sim_mem->data);
	if ( sim_mem->event_data )
		free(sim_mem->event_data);
	if ( sim_mem->cvode_mem ) 
		CVodeFree(&(sim_mem->cvode_mem));
	if ( sim_mem->x )
		N_VDestroy_Serial(sim_mem->x);	
	if ( sim_mem->atolV )
		N_VDestroy_Serial(sim_mem->atolV);
	if ( sim_mem->sx )
		N_VDestroyVectorArray_Serial(sim_mem->sx, sim_mem->np);
	if ( sim_mem->atols_ss )
		N_VDestroy_Serial(sim_mem->atols_ss);
	if ( sim_mem->atolV_ss )
		N_VDestroyVectorArray_Serial(sim_mem->atolV_ss, sim_mem->np);

	/* SSA memory */
	if ( sim_mem->x_lb )
		N_VDestroy_Serial(sim_mem->x_lb);
	if ( sim_mem->x_ub )
		N_VDestroy_Serial(sim_mem->x_ub);
    
    free( sim_mem );
}