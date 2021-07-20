% Here, we consider Ackley's optimization test function, which is defined
% for arbitrary dimension and possesses many local minima with one steep
% global minimum fval = 0 at x = [0,...,0].

% consider the problem in dimension 100

dim = 100;

% create problem struct

problem = struct();
problem.f = 'ex6';
problem.x_L = -20*ones(dim,1);
problem.x_U = 20*ones(dim,1);

% create options struct

options = struct();
options.maxeval = 2000;

% as local solver, we consider ydhc. Alternatively, try e.g. dhc and
% fmincon, or any other solver.
options.local.solver = 'ydhc';
options.local.iterprint = 0;

% run optimization

results = MEIGO(problem, options, 'ESS');
