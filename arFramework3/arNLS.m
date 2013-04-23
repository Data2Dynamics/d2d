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
method = 2;

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
optimValues(1).iteration = 0;
optimValues(1).mu = mu;
optimValues(1).normdp = nan;

dp = nan;
dresnorm = -1;
while(iter < options.MaxIter && dresnorm < 0 && mu >= options.TolX)
    iter = iter + 1;
    
    % call output function
    if(~isempty(options.OutputFcn))
        feval(options.OutputFcn,[],optimValues,'iter');
        optimValues(1).iteration = optimValues(1).iteration + 1;
        optimValues(1).mu = mu;
        optimValues(1).normdp = norm(dp);
    end
    
    % solve subproblem
    [dp, solver_calls, qred, dpmem, grad_dir_frac, resnorm_expect] = ...
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
        printiter(iter, options.MaxIter, resnorm, mu, norm(dp), dresnorm, ...
            0, sum(qred), grad_dir_frac, dresnorm_expect);
    end
    
    if(mu >= options.TolX)
        % adjust mu
        did_shric = false;
        while(dresnorm>=0)
            % output
            if(debug>2)
                fprintf('  -\n');
            end
            
            mu = mu / mu_fac; % shrinc trust region
            did_shric = true;
            
            % solve subproblem
            [dp, solver_calls, qred, dpmem, grad_dir_frac, resnorm_expect] = ...
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
                printiter(iter, options.MaxIter, resnorm, mu, norm(dp), dresnorm, ...
                    0, sum(qred), grad_dir_frac, dresnorm_expect);
            end
            
            if(mu < options.TolX)
                break;
            end
        end
        
        if(mu > options.TolX)
%             %  shrinc trust region if approximation is bad
%             dresnorm_rel = dresnorm / dresnorm_expect;
%             if(dresnorm_rel<0 || abs(dresnorm_rel-1)>0.25)
%                 mu = mu / mu_fac;
%                 did_shric = true;
%             end
            
            % enlarge mu if not shricted before
            if(~did_shric)
                if(norm(dp)>0.5*mu)
                    mu = mu * mu_fac;
                end
                if(mu>mu_max)
                    mu = mu_max;
                end
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
    end
    
    % output
    if(debug>2)
        if(dresnorm<0)
            fprintf('  *\n');
        else
            fprintf('  -\n');
        end
    end
end

% exit condition
if(iter==options.MaxIter)
    exitflag = 0;
end
if(mu<options.TolX)
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



function printiter(iter, maxIter, resnorm, mu, norm_dp, dresnorm, ...
    norm_gred, dim_red, grad_dir_frac, dresnorm_expect)
fprintf('%3i/%3i  resnorm=%-8.2g  mu=%-8.2g  norm(dp)=%-8.2g  ', iter, maxIter, resnorm, mu, norm_dp);

if(dresnorm<0)
    fprintf('dresnorm=%-8.2g  ', dresnorm);
else
    fprintf('dresnorm=%8s  ', '');
end
fprintf('approx_qual=%-5.2f  norm(g)=%-8.2g  dim_red=%i  grad_dir=%3.1f', ...
    dresnorm/dresnorm_expect, norm_gred, dim_red, grad_dir_frac);
