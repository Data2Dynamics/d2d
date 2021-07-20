
function Results_multistart = ssm_multistart(problem,opts,varargin)

% Function   : ssm_multistart
% Written by : Process Engineering Group IIM-CSIC (jegea@iim.csic.es)
% Created on : 01/07/2005
% Last Update: 05/11/2014 by A.Gabor
%
% Script to do multistart optimization (e.g. Multiple local optimization
% starting from different initial points uniformly distributed withing the
% bounds)
%
%
%          Results_multistart=ssm_multistart(problem,opts,p1,p2,...,pn)
%
% Input Parameters:
%         problem   - Structure containing problem
%
%               problem.f               = Name of the file containing the objective function
%               problem.x_L             = Lower bounds of decision variables
%               problem.x_U             = Upper bounds of decision variables
%               problem.x_0             = Initial point(s) (optional)
%
%               Additionally, fill the following fields if your problem has
%               non-linear constraints
%
%               problem.neq             = Number of equality constraints (do not define it if there are no equality constraints)
%               problem.c_L             = Lower bounds of nonlinear constraints
%               problem.c_U             = Upper bounds of nonlinear constraints
%               problem.int_var         = Number of integer variables
%               problem.bin_var         = Number of binary variables
%               NOTE: The order of decision variables is x=[cont int bin]
%
%         opts      - Structure containing options
%
%           opts.ndiverse: Number of uniformly distributed points inside the bounds to do the multistart optimization (default 100)
%
%           Local options
%               opts.local.solver             = Choose local solver ( default 'fmincon')
%                                                     0: Local search deactivated
%                                                     'fmincon' (Mathworks)
%                                                     'clssolve' (Tomlab)
%                                                     'snopt' (Tomlab)
%                                                     'nomad' (Abramson)
%                                                     'npsol' (Tomlab)
%                                                     'solnp' (Ye)
%                                                     'n2fb' (For parameter stimation problems. Single precision)
%                                                     'dn2fb' (For parameter stimation problems. Double precision)
%                                                     'dhc' (Yuret) Direct search
%                                                     'fsqp'
%                                                     'ipopt'
%                                                     'misqp' (Exler) For mixed integer problems
%               opts.local.tol                = Level of tolerance in local search
%               opts.local.iterprint          = Print each iteration of local solver on screen
%
%   p1,p2... :  optional input parameters to be passed to the objective
%   function
%
%
% Output Parameters:
%         Results_multistart      -Structure containing results
%
%               Results_multistart.fbest       :Best objective function value found after the multistart optimization
%               Results_multistart.xbest       :Vector providing the best function value
%               Results_multistart.x0          :Array containing the vectors used for the multistart optimization (in rows)
%               Results_multistart.f0          :Vector containing the objective function values of the initial solutions used in the multistart optimization
%               Results_multistart.func        :Vector containing the objective function values obtained after every local search
%               Results_multistart.xxx         :Array containing the vectors provided by the local optimization (in rows)
%               Results_multistart.no_conv     :Array containing the initial points that did not provide any solution when the local search was applied (in rows)
%               Results_multistart.nfuneval    :Vector containing the number of function evaluations in every optimisation
%               Results_multistart.time        :Total CPU time to carry out the multistart optimization

%               NOTE: To plot an histogram of the results: hist(Results_multistart.func)

global input_par
cpu_time = cputime;
fin = 0;

%Extra input parameters for fsqp and n2fb
if nargin>2
    if  strcmp(opts.local.solver,'dn2fb') | strcmp(opts.local.solver,'n2fb') |...
            strcmp(opts.local.solver,'fsqp') | strcmp(opts.local.solver,'nomad') | strcmp(opts.local.solver,'nl2sol')
        input_par = varargin;
    end
end

x_U = problem.x_U;
x_L = problem.x_L;

if not(isfield(problem,'neq')) | isempty(problem.neq)
    neq = 0;
else
    neq = problem.neq;
end

if not(isfield(problem,'x_0'))
    x_0 = [];
else
    x_0 = problem.x_0;
end

if not(isfield(problem,'c_U'))
    c_U = [];
    c_L = [];
else
    c_U = problem.c_U;
    c_L = problem.c_L;
end

if not(isfield(problem,'int_var')) | isempty(problem.int_var)
    int_var = 0;
else
    int_var = problem.int_var;
end

if not(isfield(problem,'bin_var')) | isempty(problem.bin_var)
    bin_var = 0;
else
    bin_var = problem.bin_var;
end

if not(isfield(problem,'ndata'))
    ndata = [];
else
    ndata = problem.ndata;
end

%Load default values & set options:
default      = ssm_defaults;
opts         = ssm_optset(default,opts);
weight       = opts.weight;
tolc         = opts.tolc;
local_solver = opts.local.solver;
% maxtime      = opts.maxtime;
% maxeval      = opts.maxeval;

%Check if bounds have the same dimension
if length(x_U) ~= length(x_L)
    disp('Upper and lower bounds have different dimension!!!!')
    disp('EXITING')
    Results_multistart = [];
    return
else
    %Number of decision variables
    nvar = length(x_L);
end

%Transformacion de variables de cadena en funciones para que la llamada con
%feval sea mas rapida
fobj = str2func(problem.f);

isFjacDefined = 0;
if isfield(problem,'fjac') & ~isempty(problem.fjac)
    isFjacDefined = 1;
    if ischar(problem.fjac)
        fjac = str2func(problem.fjac);
    elseif isa(problem.fjac,'function_handle')
        fjac = problem.fjac;
    else
        error('fjac is neither a string nor a function handle.')
    end
end

if (int_var + bin_var) & local_solver & not(strcmp(local_solver,'misqp'))
    fprintf('For problems with integer and/or binary variables you must use MISQP as a local solver \n');
    fprintf('EXITING \n');
    Results_multistart = [];
    return
end

if local_solver & strcmp(local_solver,'n2fb') | strcmp(local_solver,'dn2fb') | strcmp(local_solver,'nl2sol')
    n_out_f = nargout(problem.f);
    if n_out_f<3
        fprintf('%s requires 3 output arguments \n',local_solver);
        fprintf('EXITING \n');
        Results_multistart = [];
        return
    else
        %Generate a random point within the bounds
        randx   = rand(1,nvar).*(x_U-x_L)+x_L;
        [f g R] = feval(fobj,randx,varargin{:});
        ndata   = length(R);
    end
else
    ndata=[];
end

%If there are equality constraints
if neq
    %Set the new bounds for all the constraints
    c_L = [-tolc*ones(1,neq) c_L];
    c_U = [tolc*ones(1,neq) c_U];
end
nconst = length(c_U);
if nconst
    n_out_f = nargout(problem.f);
    if n_out_f < 2
        fprintf('For constrained problems the objective function must have at least 2 output arguments \n');
        fprintf('EXITING \n');
        Results_multistart=[];
        return
    end
end

% Load aux files:
if isFjacDefined
    ssm_aux_local(problem.f,problem.fjac,x_L,x_U,c_L,c_U,neq,local_solver,nvar,varargin{:});
else
    ssm_aux_local(problem.f,[],x_L,x_U,c_L,c_U,neq,local_solver,nvar,varargin{:});
end

if opts.ndiverse
    multix = rand(opts.ndiverse,nvar);    
    a      = repmat(x_U-x_L,[opts.ndiverse,1]);
    b      = repmat(x_L,[opts.ndiverse,1]);
    multix = multix.*a+b; % to put the variables inside the bounds
else
    multix = [];
end
multix = [x_0; multix];
if (int_var) || (bin_var)
    multix = ssm_round_int(multix,int_var+bin_var,x_L,x_U);
end
nrun     = size(multix,1);
f0       = zeros(nrun,1);
tic
f0(1)    = feval(problem.f,multix(1,:),varargin{:}); 
tf0      = toc;
func     = zeros(nrun,1);
xxx      = zeros(nrun,nvar);
time     = zeros(nrun,1);
nfuneval = zeros(nrun,1);
no_conv  = [];
x0       = [];

for ind_search = 1:nrun
    x0 = multix(ind_search,:);  
    if strcmp(local_solver,'mix')
        [this_local_solver,obj_fun] = ssm_mix(ind_search,nrun);
        this_f                      = obj_fun;
        fobj                        = str2func(this_f);        
        %-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        % If the solver from the 'mix' is NL2SOL:
        %-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        if nargin > 2 
            if  strcmp(this_local_solver,'nl2sol'), input_par = varargin; end
        end
        if strcmp(this_local_solver,'nl2sol')
            opts.local.use_gradient_for_finish = 0;
            opts.local.check_gradient_for_finish = 0;
            n_out_f = nargout(this_f);
            if n_out_f<3
                fprintf('%s requires 3 output arguments \n',this_local_solver);
                fprintf('EXITING \n');
                Results_multistart = [];
                return
            else
                %Generate a random point within the bounds
                randx   = rand(1,nvar).*(x_U-x_L)+x_L;
                [f g R] = feval(fobj,randx,varargin{:});
                ndata   = length(R);
            end
        else
            ndata=[];
        end
        % If we are using NL2SOL and the Jacobian is provided:
        isFjacDefined = 0;
        if isfield(problem,'fjac') & ~isempty(problem.fjac)
            isFjacDefined = 1;
            if ischar(problem.fjac)
                fjac = str2func(problem.fjac);
            elseif isa(problem.fjac,'function_handle')
                fjac = problem.fjac;
            else
                error('fjac is neither a string nor a function handle.')
            end
        end
        if isFjacDefined
            ssm_aux_local(this_f,problem.fjac,x_L,x_U,c_L,c_U,neq,this_local_solver,nvar,varargin{:});
        else
            ssm_aux_local(this_f,[],x_L,x_U,c_L,c_U,neq,this_local_solver,nvar,varargin{:});
        end
        %-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        % If the solver from the 'mix' is DHC:
        %-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        opts.local.use_gradient_for_finish   = 0; % provide gradient?
        opts.local.check_gradient_for_finish = 0; % check gradient?      
        %-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        %-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    else
        this_f = problem.f;
        this_local_solver = opts.local.solver;
    end   
    f0(ind_search) = feval(this_f,x0,varargin{:});
    if opts.local.iterprint
        fprintf('\n');
        fprintf('Local search number: %i \n',ind_search);
        fprintf('Call local solver: %s \n', upper(this_local_solver))
        fprintf('Initial point function value: %f \n',f0(ind_search));
    end
    tic
    %[x,fval,exitflag,numeval]=ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,[],...
    %    this_local_solver,opts.local.iterprint,opts.local.tol,weight,nconst,tolc,opts.local,varargin{:});
%%% AFV start handling errors: %%%
  try
    if isFjacDefined
        [x,fval,exitflag,numeval]=ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,fjac,...
            this_local_solver,opts.local.iterprint,opts.local.tol,weight,nconst,tolc,opts.local,varargin{:});
    else
        [x,fval,exitflag,numeval]=ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,[],...
            this_local_solver,opts.local.iterprint,opts.local.tol,weight,nconst,tolc,opts.local,varargin{:});
    end
  catch
    fprintf('This local search produced an error \n');
    x        = x0;
    fval     = inf;
    exitflag = -1;
    numeval  = 1;
  end
%%% AFV end of error handling code %%%
    time(ind_search) = toc;
    if opts.local.iterprint
        fprintf('Local solution function value: %f \n',fval);        
        if exitflag<0
            fprintf('Warning!!! This could be an infeasible solution: \n');
        end        
        fprintf('Number of function evaluations in the local search: %i \n',numeval);
        fprintf('CPU Time of the local search: %f  seconds \n\n',time(ind_search));
    end
    
    if exitflag<0 | isnan(fval) | isinf(fval)
        no_conv = [no_conv;x0];
    end
    func(ind_search) = fval;
    xxx(ind_search,:) = x;
    nfuneval(ind_search) = numeval;
    
%     if sum(nfuneval) > maxeval
%         fprintf('Multistart stopped after %d starts. Maximum number of function evaluation (%d) is reached.\n',ind_search,maxeval);
%         fin = 1;
%         break;
%     elseif sum(time) > maxtime
%         fprintf('Multistart stopped after %d starts. Maximum CPU computation time (%d) is reached.\n',ind_search,maxtime);
%         fin = 2;
%         break;
%     end
end

[aaa,iii]                = min(func);
Results_multistart.x0    = multix;
Results_multistart.f0    = f0;
Results_multistart.func  = func;
Results_multistart.xxx   = xxx;
Results_multistart.fbest = aaa;
if ~isempty(xxx)
    Results_multistart.xbest = xxx(iii,:);
elseif ~isinf(f0)
    disp('No feasible improvement of the cost function was achieved.')
    Results_multistart.xbest = x0;
else
    Results_multistart.xbest = [];
end
Results_multistart.no_conv  = no_conv;
Results_multistart.nfuneval = nfuneval;
Results_multistart.time     = cputime-cpu_time;  % total CPU time
Results_multistart.itime    = time;           % CPU time of each search.
Results_multistart.end_crit = fin;
% Compute the  convergence curve
% initial cost and time:
% [fconv(1) I]= min(f0);
fconv(1) = f0(1);
tconv(1) = tf0;
com_time = tf0; % commulative time
bestit   = multix(1,:);
% add new point if the cost function is improved:
for i = 1:length(func)
    com_time = com_time + time(i); % add the CPU time.    
    if func(i) < fconv(end)
        fconv = [fconv;func(i)];
        bestit = [bestit;xxx(i,:)];
        tconv = [tconv;com_time];
    elseif i == length(func)
        % if not better, but the last point is processed
        % duplicate the last point at the last time
        fconv = [fconv;fconv(end)];
        tconv = [tconv;com_time];
        bestit = [bestit;bestit(end,:)];
    end
end

Results_multistart.convcurve = [tconv, fconv];
Results_multistart.bestit = bestit;
% save ssm_multistart_report problem opts Results_multistart
% ssm_delete_files(local_solver,c_U);

return
