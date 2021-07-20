% Optimization settings:
nthreads    = 2;            % number of threads
n_iter      = 2;            % number of cooperative iterations
is_parallel = true;         % parallel (true) or sequential (false)
maxtime_per_iteration = 10; % time limit for each iteration

% Number of parameters; bounds; and initial point:
npars = 50;
x_L   =-5*ones(1,npars);
x_U   = 5*ones(1,npars);
x_0   = round(rand(1,npars).*(x_U-x_L+1)+x_L-0.5);

% Read default optimization settings:
par_struct=get_CeVNS_options(nthreads,npars,maxtime_per_iteration);

% Overwrite the following default options in par_struct:
for i=1:nthreads		
	par_struct(i).problem.f   = 'rosen10';
	par_struct(i).problem.x_L = x_L;
	par_struct(i).problem.x_U = x_U;
	par_struct(i).problem.x_0 = x_0;	
end

% Run CeVNS:
Results = CeVNS(par_struct,n_iter,is_parallel)
