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
%  3 = gradient descent (up to cauchy point)
%  4 = dogleg
%  5 = generalized trust region (based on modified trust.m)
%  6 = trdog
%  7 = trdog pcgr

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
if(method==5)
    mu = eye(length(p))*mu;
end

% trust region size scale factor
mu_fac = 2;

% counters
iter = 0;
funevals = 0;
solver_calls = 0;

% initial function evaluation
[res, sres] = feval(fun, p);
resnorm = sum(res.^2);
resnorm_start = resnorm;
funevals = funevals + 1;

llh = sum(res.^2);      % objective function
g = -2*res*sres;        % gradient
H = 2*(sres'*sres);     % Hessian matrix

% output
if(debug>2)
    fprintf('%3i/%3i  resnorm=%-8.2g\n', iter, options.MaxIter, resnorm);
end
optimValues = struct([]);
optimValues(1).iteration = 0;
optimValues(1).mu = mu;
optimValues(1).normdp = nan;
if(~isempty(options.OutputFcn))
    feval(options.OutputFcn,[],optimValues,'iter');
end

q_converged = false;
while(iter < options.MaxIter && ~q_converged)
    iter = iter + 1;
    
    % solve subproblem - get trial point
    [dp, solver_calls, qred, dpmem, grad_dir_frac, resnorm_expect, normdpmu] = ...
        arNLSstep(llh, g, H, sres, mu, p, lb, ub, solver_calls, dpmem, useInertia, method);
    pt = p + dp;
    
    % ensure strict feasibility - the hard way
    pt(pt<lb) = lb(pt<lb);
    pt(pt>ub) = ub(pt>ub);
    
    % evaluate trial point
    [rest, srest] = feval(fun, pt);
    resnormt = sum(rest.^2);
    funevals = funevals + 1;
    
    % fit improvement statistics
    dresnorm = resnormt - resnorm; 
    dresnorm_expect = resnorm_expect - resnorm;
    
    % approximation quality
    approx_qual = dresnorm/dresnorm_expect;
    q_approx_qual = approx_qual > 0.75;
    
    % reduction achieved ?
    q_reduction = dresnorm<0;
    
%     q_accept_step = q_reduction;
    q_accept_step = q_reduction && q_approx_qual;
    % update if step was accepted
    if(q_accept_step)
        firstorderopt = norm(g(~qred));
        
        p = pt;
        
        res = rest;
        sres = srest;
        resnorm = resnormt;
        
        llh = sum(res.^2);      % objective function
        g = -2*res*sres;        % gradient
        H = 2*(sres'*sres);     % Hessian matrix
        
        % call output function
        if(~isempty(options.OutputFcn))
            optimValues(1).iteration = optimValues(1).iteration + 1;
            optimValues(1).mu = mu;
            optimValues(1).normdp = norm(dp);
            feval(options.OutputFcn,[],optimValues,'iter');
        end
    else
        firstorderopt = nan;
    end
    
    % update schedule for trust region
    q_enlarge = q_reduction && q_approx_qual && (normdpmu > 0.9 || isnan(normdpmu));
    q_shrink = ~q_reduction || ~q_approx_qual;
    
    mu_old = mu;
    dmu = 0;
    if(q_enlarge) % enlarge trust region
        mu = arNLSTrustTrafo(mu, mu_fac, dp, false);
        dmu = 1;
    elseif(q_shrink) % shrink trust region
        mu_red_fac = 1/mu_fac;
        if(isscalar(mu) && ~isnan(normdpmu))
            mu_red_fac = min([mu_red_fac mu_red_fac*normdpmu]);
        end
        mu = arNLSTrustTrafo(mu, mu_red_fac, dp, false);
        dmu = -1;
    end
   
    % output
    if(debug>2)
        printiter(iter, options.MaxIter, resnorm, mu_old, dmu, norm(dp), normdpmu, dresnorm, ...
            firstorderopt, sum(qred), grad_dir_frac, approx_qual, q_accept_step);
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
function printiter(iter, maxIter, resnorm, mu, dmu, norm_dp, normdpmu, dresnorm, ...
    norm_gred, dim_red, grad_dir_frac, approx_qual, step_accept)

if(~step_accept)
    outstream = 2;
else
    outstream = 1;
end

if(dmu==-1)
    dmu = '-';
elseif(dmu==+1)
    dmu = '+';
else
    dmu = '0';
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
     
    fprintf(outstream, 'mu=%-8.2g %s (det=%-8.2g cond=%-8.2g maxeig=%-8.2g)  ', normdpmu, dmu, det(mu), cond(mu), max(eig(mu)));
else
    fprintf(outstream, 'mu=%-8.2g %s ', mu, dmu);
end
fprintf(outstream, 'norm(dp)=%-8.2g  dresnorm=%-8.2g  ', norm_dp, dresnorm);
fprintf(outstream, 'approx_qual=%-5.2f  norm(g)=%-8.2g  dim_red=%i  grad_dir=%3.1f\n', ...
    approx_qual, norm_gred, dim_red, grad_dir_frac);
