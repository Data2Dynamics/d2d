function Results = CeSS_thread(par_struct,Nruns)

    %Pass ess options to ess_kernel:
    ess_options = []; 
    ess_options.ndiverse   = par_struct(Nruns).opts.ndiverse;
    ess_options.maxtime    = par_struct(Nruns).opts.maxtime_cess;
    ess_options.maxeval    = par_struct(Nruns).opts.maxeval;
    ess_options.log_var    = par_struct(Nruns).opts.log_var;
    ess_options.local      = par_struct(Nruns).opts.local;
    ess_options.dim_refset = par_struct(Nruns).opts.dim_refset;
    
    %Run optimization:
    res = ess_kernel(par_struct(Nruns).problem,ess_options);

    %Store results:
    Results.x        = res.x;
    Results.f        = res.f;
    Results.refset_x = res.Refset.x;
    Results.refset_f = res.Refset.f;
    Results.neval    = res.neval;
    Results.numeval  = res.numeval;
    Results.cpu_time = res.cpu_time;
    Results.fbest    = res.fbest;
    Results.xbest    = res.xbest;
    Results.time     = res.time;
		
end