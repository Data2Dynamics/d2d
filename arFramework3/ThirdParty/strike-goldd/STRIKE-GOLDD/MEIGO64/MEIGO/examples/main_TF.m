% In this example, we use problems from the TF class, which basically
% contains diverse optimization test problems.

% Feel free to try out also other problems and settings, e.g.
% non-differentiable problems and the impact of noise.

% test problem
test_problem = TF.hyperellipse;

% dimension
dim = 20;

% get values for dimension
[lb,ub] = TF.f_getVectors(test_problem, dim);

% set up MEIGO

problem = struct();
problem.f = test_problem.fun; % here we pass a function handle
problem.x_L = lb;
problem.x_U = ub;

options = struct();
options.maxeval = 200*dim;
options.local.solver = 'fmincon';
options.local.finish = 'fmincon';
options.local.iterprint = 1;

% run optimization

results = MEIGO(problem, options, 'ESS');
