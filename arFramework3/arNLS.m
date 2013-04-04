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

% initial trust region size
if(isempty(options.InitTrustRegionRadius))
    mu = 1;         
else
    mu = options.InitTrustRegionRadius;         
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

% step inertia
dpmem = [];

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
    [dp, solver_calls, gred, dpmem] = ...
        getStep(g, H, mu, p, lb, ub, solver_calls, dpmem);
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
    
    % output
    if(debug>2)
        if(dresnorm<0)
            fprintf('%3i/%3i  resnorm=%-8.2g  mu=%-8.2g  norm(dp)=%-8.2g  dresnorm=%-8.2g  norm(g)=%-8.2g  sub-dim=%i', ...
                iter, options.MaxIter, resnorm, mu, norm(dp), dresnorm, norm(gred), length(gred));
        else
            fprintf('%3i/%3i  resnorm=%-8.2g  mu=%-8.2g  norm(dp)=%-8.2g  dresnorm=%8s  norm(g)=%-8.2g  sub-dim=%i', ...
                iter, options.MaxIter, resnorm, mu, norm(dp), '', norm(gred), length(gred));
        end
    end
    
    if(mu >= options.TolX)
        % adjust mu
        if(dresnorm>=0)
            while(dresnorm>=0)
                % output
                if(debug>2)
                    fprintf('  -\n');
                end
                
                mu = mu / mu_fac; % shrinc trust region
                
                % solve subproblem
                [dp, solver_calls, gred, dpmem] = ...
                    getStep(g, H, mu, p, lb, ub, solver_calls, dpmem);
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
                
                % output
                if(debug>2)
                    if(dresnorm<0)
                        fprintf('%3i/%3i  resnorm=%-8.2g  mu=%-8.2g  norm(dp)=%-8.2g  dresnorm=%-8.2g  norm(g)=%-8.2g  sub-dim=%i', ...
                            iter, options.MaxIter, resnorm, mu, norm(dp), dresnorm, norm(gred), length(gred));
                    else
                        fprintf('%3i/%3i  resnorm=%-8.2g  mu=%-8.2g  norm(dp)=%-8.2g  dresnorm=%8s  norm(g)=%-8.2g  sub-dim=%i', ...
                            iter, options.MaxIter, resnorm, mu, norm(dp), '', norm(gred), length(gred));
                    end
                end
                
                if(mu < options.TolX)
                    break;
                end
            end
        else
            if(norm(dp)>0.5*mu)
                mu = mu * mu_fac;
            end
        end
    end
    
    % update if step lead to reduction
    if(dresnorm<0)
        p = pt;
        
        res = rest;
        sres = srest;
        resnorm = resnormt;
        
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
    output.firstorderopt = norm(gred);
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


% generate trial steps
%
% g:    gradient
% H:    Hessian matrix
% mu:   trust region size
% p:    current parameters
% lb:   lower bounds
% ub:   upper bounds
function [dp, solver_calls, gred, dpmem, imem, dpmemsize] = getStep(g, H, mu, p, lb, ub, solver_calls, dpmem, imem, dpmemsize)

% solve subproblem
dp = getDP(g, H, mu);
solver_calls = solver_calls + 1;

distp = -[p-ub; -(p - lb)]; % distance to bounds (should be always positive)
onbound = distp==0; % which parameter is exactly on the bound ?
exbounds = [p+dp-ub; -(p+dp - lb)] > 0; % dp bejond bound ?

% if p on bounds and dp pointing bejond bound, reduce problem
gred = g;
qred = sum(onbound & exbounds,1)>0;
if(sum(~qred)==0)
    error('solution outside bounds, no further step possible');
end
qred_presist = qred;
while(sum(qred)>0)
    gred = g(~qred_presist);
    Hred = H(~qred_presist,~qred_presist);
    
    % solve reduced subproblem
    dp_red = getDP(gred, Hred, mu);
    solver_calls = solver_calls + 1;
    
    dp(:) = 0;
    dp(~qred_presist) = dp_red;
    
    exbounds = [p+dp-ub; -(p+dp - lb)] > 0; % dp bejond bound ?
    qred = sum(onbound & exbounds,1)>0;
    qred_presist = qred | qred_presist;
    
    if(sum(~qred_presist)==0)
        error('solution outside bounds, no further step possible');
    end
end

% if dp too long cut to bounds
dptmp = [dp; dp];
dpredfac = abs(min(distp(exbounds)./dptmp(exbounds)));
if(~isempty(dpredfac))
    if(dpredfac==0)
        error('zero step size');
    end
    dp = dp * dpredfac;
end

% memory
if(~isempty(dpmem))
%     dp = mean([dpmem;dp]);
end
dpmem = dp;



% solver function
function dp = getDP(g, H, mu)

% % trust region solution
% dp = trust(-g',H,mu)'; 
% % PROBLEM: ensuring norm(dp)<=mu in trust function
% if(norm(dp)>mu)
%     dp = dp/norm(dp)*mu;
% end

% trust region solution
dp = trust(-g',H,mu)';
lambda = 1e-6;
while(norm(dp)-mu > 1e-6)
%     fprintf('trust problem %g %g\n', norm(dp)-mu, lambda);
    lambda = lambda * 10;
    dp = trust(-g',H+lambda*eye(size(H)),mu)'; 
end

% % levenberg-marquardt
% dp = transpose(pinv(H + eye(size(H))/mu)*g');
