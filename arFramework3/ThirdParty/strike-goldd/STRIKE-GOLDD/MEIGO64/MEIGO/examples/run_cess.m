nthreads=8;
npars=1000;
n_iter=4;
is_parallel=false;
%THis is too short for this problem, only for testing purposes
maxtime_per_iteration=30;

par_struct=get_CeSS_options(nthreads,npars,maxtime_per_iteration);
x_L=-100*ones(1,npars);
x_U= 100*ones(1,npars);

for i=1:nthreads
		
	par_struct(i).problem.f='f17';
	par_struct(i).problem.x_L=x_L;
	par_struct(i).problem.x_U=x_U;
	%DHC takes as bit to run, uncoment the next two lines for tests
	par_struct(i).opts.local.solver  = 0; 
	par_struct(i).opts.local.finish  = 0;
	
end

Results = CeSS(par_struct,n_iter,is_parallel)

plot(Results.time,log10(Results.f))

