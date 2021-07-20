nthreads=2;
n_iter=2;
is_parallel=true;
npars=50;
maxtime_per_iteration=10;

x_L=-5*ones(1,npars);
x_U= 5*ones(1,npars);
x_0= round(rand(1,npars).*(x_U-x_L+1)+x_L-0.5)

par_struct=get_CeVNS_options(nthreads,npars,maxtime_per_iteration);

for i=1:nthreads
		
	par_struct(i).problem.f='rosen10';
	par_struct(i).problem.x_L=x_L;
	par_struct(i).problem.x_U=x_U;
	par_struct(i).problem.x_0=x_0;
	
end

Results = CeVNS(par_struct,n_iter,is_parallel)

