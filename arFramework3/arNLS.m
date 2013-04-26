% [p,resnorm,res,exitflag,output] = arNLS(fun,p,lb,ub,options,varargin)
%
% use like LSQNONLIN
% 
% exitflag:
%  0  Too many function evaluations or iterations.
%  1  Converged to a solution.
%  2  Change in p too small.
%  3  Change in resnorm too small.
%  4  Computed search direction too small.

function [p,resnorm,res,exitflag,output,lambda,jac] = arNLS(fun,p,lb,ub,options)

% check bounds
if (sum(ub<=lb)>0)
    error('nls_trust: bounds are inconsistent');
end
if(sum(p>ub | p<lb)>0)
    error('nls_trust: parameters inconsistent with bounds');
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

% method:   0 = trust region (based on modified trust.m)
%           1 = Levenberg-Marquardt
%           2 = Newton (with maximal step length mu)
%           3 = gradient descent (up to cauchy point)
%           4 = dogleg
%           5 = generalized trust region (based on modified trust.m)
method = 5;

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

% maximum trust region size
mu_max = Inf;

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
if(~isempty(options.OutputFcn))
    optimValues(1).iteration = 0;
    optimValues(1).mu = nan;
    optimValues(1).normdp = nan;
    feval(options.OutputFcn,[],optimValues,'iter');
end

dp = 1;
dresnorm = -1;
while(iter < options.MaxIter && dresnorm < 0 && doContinue(mu, dresnorm, norm(dp), options))
    iter = iter + 1;
    
    % solve subproblem
    [dp, solver_calls, qred, dpmem, grad_dir_frac, resnorm_expect, mudp] = ...
        arNLSstep(llh, g, H, mu, p, lb, ub, solver_calls, dpmem, useInertia, method);
    pt = p + dp;
    
    % ensure strict feasibility
    pt(pt<lb) = lb(pt<lb);
    pt(pt>ub) = ub(pt>ub);
    
    % function evaluation
    [rest, srest] = feval(fun, pt);
    resnormt = sum(rest.^2);
    funevals = funevals + 1;
    
    % fit improve
    dresnorm = resnormt - resnorm; 
    dresnorm_expect = resnorm_expect - resnorm;
    
    % output
    if(debug>2)
        printiter(iter, options.MaxIter, resnorm, mudp, mu, norm(dp), dresnorm, ...
            0, sum(qred), grad_dir_frac, dresnorm_expect);
    end
    
    if(doContinue(mu, dresnorm, norm(dp), options))
        % adjust mu - shrinc trust region if approximation is bad
        did_shric = false;
        while(doContinueApprox(dresnorm, dresnorm_expect))
            % output
            if(debug>2)
                fprintf('  -\n');
            end
            
            % shrinc trust region
            mu = arNLSTrustTrafo(mu, mu_fac, dp, true);
            did_shric = true;
            
            % solve subproblem
            [dp, solver_calls, qred, dpmem, grad_dir_frac, resnorm_expect, mudp] = ...
                arNLSstep(llh, g, H, mu, p, lb, ub, solver_calls, dpmem, useInertia, method);
            pt = p + dp;
            
            % ensure strict feasibility
            pt(pt<lb) = lb(pt<lb);
            pt(pt>ub) = ub(pt>ub);
            
            % function evaluation
            [rest, srest] = feval(fun, pt);
            resnormt = sum(rest.^2);
            funevals = funevals + 1;
            
            % fit improve
            dresnorm = resnormt - resnorm;
            dresnorm_expect = resnorm_expect - resnorm;
            
            % output
            if(debug>2)
                printiter(iter, options.MaxIter, resnorm, mudp, mu, norm(dp), dresnorm, ...
                    0, sum(qred), grad_dir_frac, dresnorm_expect);
            end
            
            if(~doContinue(mu, dresnorm, norm(dp), options))
                break;
            end
        end
    end
    
    % update if step lead to reduction
    if(dresnorm<0)
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
            optimValues(1).mu = mudp;
            optimValues(1).normdp = norm(dp);
            feval(options.OutputFcn,[],optimValues,'iter');
        end
    end
    
    % output
    if(debug>2)
        if(dresnorm<0)
            fprintf('  *\n');
        else
            fprintf('  -\n');
        end
    end
    
    % adjust mu - enlarge trust region if not shricted before
    if(doContinue(mu, dresnorm, norm(dp), options) && ~did_shric)
        if(isscalar(mu))
            if(norm(dp)>0.5*mu)
                mu = mu * mu_fac;
            end
            if(mu>mu_max)
                mu = mu_max;
            end
        else
            if(norm(dp)>0.5*mudp)
                mu = arNLSTrustTrafo(mu, mu_fac, dp, false);
            end
        end
    end
end

% exit condition
if(iter==options.MaxIter)
    exitflag = 0;
end
if(~doContinue(mu, dresnorm, norm(dp), options))
    exitflag = 2;
end

% assign output structure
if (nargout>4)
    output = struct([]);
    output(1).iterations = iter;
    output.funcCount = funevals;
    output.solverCount = solver_calls;
    output.algorithm = 'nls_trust';
    output.firstorderopt = 0;
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



% check if convergence
function q = doContinue(mu, dresnorm, dpnorm, options)
if(isscalar(mu))
    q = mu >= options.TolX;
else
    q = dresnorm<-options.TolFun || dpnorm>=options.TolX;
end



% check if trust region shrinking necessary
function q = doContinueApprox(dresnorm, dresnorm_expect)
q = dresnorm>=0;
% q = dresnorm>=0 || abs((dresnorm / dresnorm_expect)-1)>0.25;



% print interation
function printiter(iter, maxIter, resnorm, mudp, mu, norm_dp, dresnorm, ...
    norm_gred, dim_red, grad_dir_frac, dresnorm_expect)
if(isscalar(mu))
    fprintf('%3i/%3i  resnorm=%-8.2g  mu=%-8.2g  norm(dp)=%-8.2g  ', iter, maxIter, resnorm, mudp, norm_dp);
else
    % plot(eig(real(mu)),'*-'); drawnow;
    fprintf('%3i/%3i  resnorm=%-8.2g  mu=%-8.2g(det=%-8.2g cond=%-8.2g maxeig=%-8.2g)  norm(dp)=%-8.2g  ', ...
        iter, maxIter, resnorm, mudp, det(mu), cond(mu), max(eig(mu)), norm_dp);
end

if(dresnorm<0)
    fprintf('dresnorm=%-8.2g  ', dresnorm);
else
    fprintf('dresnorm=%8s  ', '');
end
fprintf('approx_qual=%-5.2f  norm(g)=%-8.2g  dim_red=%i  grad_dir=%3.1f', ...
    dresnorm/dresnorm_expect, norm_gred, dim_red, grad_dir_frac);
