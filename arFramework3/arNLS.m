% [p,resnorm,res,exitflag,output] = arNLS(fun,p,lb,ub,options,method)
%
% use like LSQNONLIN
% 
% exitflag:
%  0  Too many function evaluations or iterations.
%  1  Converged to a solution.
%  2  Change in p too small.
%  3  Change in resnorm too small.
%  4  Computed search direction too small.
%
% method:   
%  0 = trust region (based on modified trust.m)
%  1 = Levenberg-Marquardt
%  2 = Newton (with maximal step length mu)
%  3 = gradient descent (with steplength mu)
%  4 = gradient descent (to cauchy point with steplength mu)
%  5 = dogleg
%  6 = generalized trust region (based on modified trust.m)
%  7 = MATLABs trdog
%  8 = Newton pcgr (with maximal step length mu)
%  9 = trdog pcgr (with maximal step length mu)
% 10 = dogleg Newton pcgr
% 11 = dogleg trdog pcgr
% 12 = trdog pcgr (no DM)
% 13 = trdog pcgr (no DG)
% 14 = trdog pcgr Levenberg-Marquardt
% 15 = trdog pcgr 2D subspace 
% 16 = trdog pcgr (no DM) 2D subspace 
% 17 = trdog pcgr (no DM) 3D subspace 

function [p,resnorm,res,exitflag,output,lambda,jac] = arNLS(fun,p,lb,ub,options,method)

if(nargin==0)
    p = arNLSstep;
    return;
end

if(~exist('method','var'))
    method = 0;
end

% check bounds
if (sum(ub<=lb)>0)
    error('arNLS: bounds are inconsistent');
end
if(sum(p>ub | p<lb)>0)
    error('arNLS: parameters inconsistent with bounds');
end

% not used
lambda = [];
jac = [];

exitflag = 1;

if(nargin>4)
    options = optimset(optimset('lsqnonlin'), options);
else
    options = optimset('lsqnonlin');
end

% output level
switch(options.Display)
    case('off')
        debug = 0;
    case('notify')
        debug = 1;
    case('final')
        debug = 2;
    case('iter')
        debug = 3;
end

% inertia effect using memory
% 0 = no
% useInertia has to be < 1
useInertia = 0;
dpmem = [];

% initial trust region size
if(isempty(options.InitTrustRegionRadius))
    mu = 1;
else
    mu = options.InitTrustRegionRadius;
end
if(method==6)
    mu = eye(length(p))*mu;
end

% trust region size scale factor
mu_fac = 2;

% counters
iter = 0;
funevals = 0;
solver_calls = 0;

% initial function evaluation
funevals = funevals + 1;
if(nargout(fun)==2 || nargout(fun)==-1)
    [res, sres] = feval(fun, p);
    H = 2*(sres'*sres);             % Hessian matrix approximation
elseif(nargout(fun)==3)
    [res, sres, H] = feval(fun, p); % user Hessian matrix
else
    error('nargout(fun)==2 or 3');
end
resnorm = sum(res.^2);
resnorm_start = resnorm;
llh = sum(res.^2);      % objective function
g = -2*res*sres;        % gradient

% calculate first order optimality criterion
onbound = [p==ub; p==lb];
exbounds = [g>0; g<0];
firstorderopt = norm(g(sum(onbound & exbounds,1)==0));

% output
if(debug>2)
    fprintf('%3i/%3i  resnorm=%-8.2g   norm(g)=%-8.2g\n', iter, options.MaxIter, resnorm, firstorderopt);
end
if(~isempty(options.OutputFcn))
    feval(options.OutputFcn,p,[],'iter');
end

q_converged = false;
while(iter < options.MaxIter && ~q_converged)
    iter = iter + 1;
    
    % solve subproblem - get trial point
    [dp, solver_calls, qred, grad_dir_frac, resnorm_expect, normdpmu_type] = ...
        arNLSstep(llh, g, H, sres, mu, p, lb, ub, solver_calls, dpmem, method);
    
    pt = p + dp;
    
    % ensure strict feasibility - the hard way
    pt(pt<lb) = lb(pt<lb);
    pt(pt>ub) = ub(pt>ub);
    
    % evaluate trial point
    try
        if(nargout(fun)==2 || nargout(fun)==-1)
            [rest, srest] = feval(fun, pt);
        elseif(nargout(fun)==3 && nargin(fun)==1)
            [rest, srest, Ht] = feval(fun, pt);
        elseif(nargout(fun)==3 && nargin(fun)==5)
            [rest, srest, Ht] = feval(fun, pt, p, res, sres, H);
        end
        resnormt = sum(rest.^2);
    catch
        resnormt = Inf;
    end
    funevals = funevals + 1;
    
    % fit improvement statistics
    dresnorm = resnormt - resnorm; 
    dresnorm_expect = resnorm_expect - resnorm;
    
    % approximation quality
    approx_qual = dresnorm/dresnorm_expect;
    q_approx_qual = approx_qual > 0.75;
    
    % calculate first order optimality criterion
    onbound = [p==ub; p==lb];
    exbounds = [g>0; g<0];
    firstorderopt = norm(g(sum(onbound & exbounds,1)==0));
    
    % reduction achieved ?
    q_reduction = dresnorm<0;
    
    % accept step ?
%     q_accept_step = q_reduction;
    q_accept_step = q_reduction && q_approx_qual;
    
    % output
    if(debug>2)
        printiter(iter, options.MaxIter, resnorm, mu, norm(dp), normdpmu_type, dresnorm, ...
            firstorderopt, find(qred), grad_dir_frac, approx_qual, q_accept_step, cond(H));
    end
    
    % call output function
    if(~isempty(options.OutputFcn))
        feval(options.OutputFcn,p,[],'iter');
    end
    
    % call plot function
    if(~isempty(options.PlotFcns))
        feval(options.PlotFcns, p, H, mu);
    end
    
    % update if step was accepted
    if(q_accept_step)    
        p = pt;
        
        res = rest;
        sres = srest;
        resnorm = resnormt;
        
        llh = sum(res.^2);      % objective function
        g = -2*res*sres;        % gradient
        if(nargout(fun)==2 || nargout(fun)==-1)
            H = 2*(sres'*sres); % Hessian matrix approximation
        elseif(nargout(fun)==3)
            H = Ht;             % user Hessian matrix
        end
        
        % inertial effect using memory
        if(~isempty(dpmem))
            dpmem = useInertia*dpmem + (1-useInertia)*dp;
        else
            dpmem = dp;
        end
    end
    
    % update schedule for trust region
    if(normdpmu_type == -1)
        q_enlarge = q_reduction && q_approx_qual && (norm(dp)/mu) > 0.9;
    elseif(normdpmu_type > 0)
        q_enlarge = q_reduction && q_approx_qual && normdpmu_type > 0.9;
    else
        q_enlarge = q_reduction && q_approx_qual;
    end
    q_shrink = ~q_reduction || ~q_approx_qual;
    
    dmu = 0;
    if(q_enlarge) % enlarge trust region
        mu = arNLSTrustTrafo(mu, mu_fac, dp, false);
        dmu = 1;
    elseif(q_shrink) % shrink trust region
        mu_red_fac = 1/mu_fac;
        if(isscalar(mu) && normdpmu_type==-1)
            mu_red_fac = min([mu_red_fac mu_red_fac*norm(dp)/mu]);
        end
        mu = arNLSTrustTrafo(mu, mu_red_fac, dp, false);
        dmu = -1;
    end
    
    % check convergence
    if(isscalar(mu))
        q_converged = mu < options.TolX || firstorderopt < 1e-6;
    else
        q_converged = norm(dp) < options.TolX || firstorderopt < 1e-6;
    end
end

% exit condition
if(iter==options.MaxIter)
    exitflag = 0;
end
if(mu < options.TolX || norm(dp) < options.TolX)
    exitflag = 2;
end
if(firstorderopt < 1e-6)
    exitflag = 1;
end

% assign output structure
if (nargout>4)
    output = struct([]);
    output(1).iterations = iter;
    output.funcCount = funevals;
    output.solverCount = solver_calls;
    output.algorithm = 'arNLS';
    output.firstorderopt = firstorderopt;
    switch(exitflag)
        case(0)
            output.message = 'Too many function evaluations or iterations.';
        case(1)
            output.message = 'Converged to a solution.';
        case(2)
            output.message = 'Change in p too small.';
    end
    output.totalimprove = resnorm_start - resnorm;
end

% output
if(debug>1 || (debug>0 && exitflag==0))
    switch(exitflag)
        case(0)
            fprintf('NLS_TRUST: Too many function evaluations or iterations.'); 
        case(1)
            fprintf('NLS_TRUST: Converged to a solution.');
        case(2)
            fprintf('NLS_TRUST: Change in p too small.');
    end
    fprintf(' Achieved improvement %g in %i/%i iterations.\n', resnorm_start - resnorm, iter, options.MaxIter);
end




% print iteration
function printiter(iter, maxIter, resnorm, mu, norm_dp, normdpmu, dresnorm, ...
    norm_gred, dim_red, grad_dir_frac, approx_qual, step_accept, condition_number)

if(~step_accept)
    outstream = 2;
else
    outstream = 1;
end

fprintf(outstream, '%3i/%3i  resnorm=%-8.2g  ', iter, maxIter, resnorm);
if(~isscalar(mu))
    figure(1)
    subplot(3,1,1)
    plot(log10(real(eig(mu))),'*-');
    subplot(3,1,2:3)
    maxmu = max(abs(mu(:)));
    imagesc(mu, [-maxmu maxmu]);
    colorbar;
    drawnow;
     
    fprintf(outstream, 'mu=%-8.2g (det=%-8.2g cond=%-8.2g maxeig=%-8.2g)  ', normdpmu, det(mu), cond(mu), max(eig(mu)));
else
    fprintf(outstream, 'mu=%-8.2g ', mu);
    if(normdpmu>0)
        fprintf(outstream, '(%-5.2f) ', normdpmu);
    end
end
fprintf(outstream, 'norm(dp)=%-8.2g  dresnorm=%-8.2g  ', norm_dp, dresnorm);
fprintf(outstream, 'approx_qual=%-8.2g  norm(g)=%-8.2g  grad_dir=%3.1f  cond(H)=%-8.2g  dim_red: ', ...
    approx_qual, norm_gred, grad_dir_frac, condition_number);
fprintf(outstream, '%i ', dim_red);
fprintf(outstream, '\n');
