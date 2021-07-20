function Results = CeVNS_thread(par_struct,Nruns)
	
    Results=[];
		
	%Reset the random seed.
	randstate = 1e3*Nruns+ sum(100*clock);
	rand('state',randstate);
	randn('state',randstate);
	
	 res=rvnds_hamming(...
        par_struct(Nruns).problem,...
        par_struct(Nruns).opts);
			
	Results.f=res.func;
	Results.time=res.time;
	Results.neval=res.neval;
	Results.numeval=res.numeval;
	Results.cpu_time=res.cpu_time;
	Results.fbest=res.fbest;
	Results.xbest=res.xbest;
		
end