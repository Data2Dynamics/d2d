% generate trial steps within box bounds
%
% [dp, solver_calls, gred, dpmem, grad_dir_frac] = ...
%       arNLSstep(llh, g, H, mu, p, lb, ub, solver_calls, ...
%       dpmem, useInertia, method, nu)
%
% llh:  likelihood
% g:    gradient
% H:    Hessian matrix
% mu:   trust region size (for LM: lambda = 10/mu)
% p:    current parameters
% lb:   lower bounds
% ub:   upper bounds
%
% method:   0 = trust region (based on modified trust.m)
%           1 = Levenberg-Marquardt
%           2 = Newton (with maximal step length mu)
%           3 = gradient descent (up to cauchy point)
%           4 = dogleg
%           5 = generalized trust region (based on modified trust.m)
%           6 = trdog

function [dp, solver_calls, qred, dpmem, grad_dir_frac, llh_expect, mudp] = ...
    arNLSstep(llh, g, H, sres, mu, p, lb, ub, solver_calls, dpmem, useInertia, method)

qred = false(size(p));

% snls like scaling
[v, dv] = definev(g,p,lb,ub);
D = diag(v);

if(method == 6)
    g2 = 0.5*g';
    pcoptions = Inf;
    pcgtol = 0.1;
    kmax = max(1,floor(length(p)/2));
    gopt = v.*g2;
    optnrm = norm(gopt,inf);
    theta = max(.95,1-optnrm);  
    
    dp = trdog(p',g2,sres,D,mu,dv,@atamult,@aprecon,...
        pcoptions,pcgtol,kmax,theta,lb',ub',[],[],'jacobprecon')';
    
    solver_calls = solver_calls + 1;
    
    % proportion of step in direction of gradient
    grad_dir_frac = dp(~qred)*g(~qred)'/norm(dp(~qred))/norm(g(~qred));
    
    % expected likelihood
    llh_expect = llh - g*dp' + 0.5*dp*H*dp';
    
    % current maximum step size
    mudp = norm(dp) / sqrt((dp/mu)*(dp/mu)');
    
    return
end

% solve subproblem
[dp, solver_calls_tmp] = getDP(g, H, mu, method, qred, D);
solver_calls = solver_calls + solver_calls_tmp;

distp = -[p-ub; -(p - lb)]; % distance to bounds (should be always positive)
onbound = distp==0; % which parameter is exactly on the bound ?
exbounds = [p+dp-ub; -(p+dp - lb)] > 0; % dp bejond bound ?

% if p on bounds and dp pointing bejond bound, reduce problem
qred = sum(onbound & exbounds,1)>0;
if(sum(~qred)==0)
    error('solution outside bounds, no further step possible');
end
qred_tmp = qred;
while(sum(qred_tmp)>0)
    gred = g(~qred);
    Hred = H(~qred,~qred);
    
    % solve reduced subproblem
    [dp_red, solver_calls_tmp] = getDP(gred, Hred, mu, method, qred, D);
    solver_calls = solver_calls + solver_calls_tmp;
    
    dp(:) = 0;
    dp(~qred) = dp_red;
    
    exbounds = [p+dp-ub; -(p+dp - lb)] > 0; % dp bejond bound ?
    qred_tmp = sum(onbound & exbounds,1)>0;
    qred = qred_tmp | qred;
    
    if(sum(~qred)==0)
        error('solution in the conner, dp pointing outside, no further step possible');
    end
end

% if dp too long cut to bounds
dptmp = [dp; dp];
dpredfac = abs(min(distp(exbounds)./dptmp(exbounds)));
if(~isempty(dpredfac))
    if(dpredfac==0)
        error('zero step size error');
    end
    dp = dp * dpredfac;
end

% inertial effect using memory
if(useInertia>0 && ~isempty(dpmem))
    dp = useInertia*dpmem + (1-useInertia)*dp;
end
dpmem = dp;

% proportion of step in direction of gradient
grad_dir_frac = dp(~qred)*g(~qred)'/norm(dp(~qred))/norm(g(~qred));

% expected likelihood
llh_expect = llh - g*dp' + 0.5*dp*H*dp';

% current maximum step size
mudp = norm(dp) / sqrt((dp/mu)*(dp/mu)');


% generate trial steps within ball of size mu (for LM: lambda = 10/mu)
%
% method:   0 = trust region (based on modified trust.m)
%           1 = Levenberg-Marquardt
%           2 = Newton (with maximal step length mu)
%           3 = gradient descent (up to cauchy point)
%           4 = dogleg

function [dp, solver_calls_tmp] = getDP(g, H, mu, method, qred, D)

solver_calls_tmp = 1;

if(mu==0)
    dp = zeros(size(g));
    return;
end

% snls like scaling
% Dred = D(~qred,~qred);
% g = g*Dred;
% H = Dred'*H*Dred;

switch method
    case 0 % trust region solution
        dp = trust(-g',H,mu)';
        
        % the function trust has a bug, therefore, in some cases
        % the problem has to be regularized
        lambda = 1e-6;
        while(norm(dp)/mu > 1.01)
            fprintf('trust.m problem %g, regularizing with new lambda=%g\n', norm(dp)/mu, lambda);
            dp = trust(-g',H+lambda*eye(size(H)),mu)';
            solver_calls_tmp = solver_calls_tmp + 1;
            lambda = lambda * 10;
        end
        
    case 1 % levenberg-marquardt
        lambda = 1/mu;
        dp = transpose(pinv(H + lambda*eye(size(H)))*g');
        
    case 2 % newton
        dp = transpose(pinv(H)*g');
        if(norm(dp)>mu)
            dp = mu*dp/norm(dp);
        end
    case 3 % gradient up to cauchy point
%         dp = norm(g)^2/(g*H*g') * g;
        dp = g;
        if(norm(dp)>mu)
            dp = mu*dp/norm(dp);
        end
    case 4
        dpC = norm(g)^2/(g*H*g') * g;
        if(norm(dpC)>mu)
            % truncated gradient
            dp = mu*g/norm(g);
%             fprintf('truncated gradient - ');
        else
            dpN = transpose(pinv(H)*g');
            if(norm(dpN)<=mu)
                % full newton
                dp = dpN;
%                 fprintf('full newton - ');
            else
                % dogleg
                a = sum((dpN - dpC).^2);
                b = 2*sum(dpC.*(dpN - dpC));
                c = sum(dpC.^2) - mu^2;
                t = (-b + sqrt(b^2 - 4 * a * c))/(2*a);
                dp = dpC + t*(dpN - dpC);
%                 fprintf('dogleg (%4.2f) - ', norm(dpC)/mu);
            end
        end
        
    case 5 % generalized trust region
        mut = mu(~qred,~qred);
        gt = g*mut;
        Ht = mut'*H*mut;
        Ht = 0.5*(Ht + Ht');
        
        dp = trust(-gt',Ht,1)';
        
        % the function trust has a bug, therefore, in some cases
        % the problem has to be regularized
        lambda = 1e-6;
        while(norm(dp) > 1.01)
            fprintf('trust.m problem %g, regularizing with new lambda=%g\n', norm(dp), lambda);
            dp = trust(-gt',Ht+lambda*eye(size(Ht)),1)';
            solver_calls_tmp = solver_calls_tmp + 1;
            lambda = lambda * 10;
        end
        dp = dp * mut; 
end

% snls like scaling
% dp = dp * Dred; 




function [v,dv]= definev(g,x,l,u)
%DEFINEV Scaling vector and derivative
%
%	[v,dv]= DEFINEV(g,x,l,u) returns v, distances to the
%   bounds corresponding to the sign of the gradient g, where
%   l is the vector of lower bounds, u is the vector of upper 
%   bounds. Vector dv is 0-1 sign vector.
%

%   Copyright 1990-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/07/06 20:45:56 $

n = length(x); 
v = zeros(n,1); 
dv=zeros(n,1);
arg1 = (g < 0)  & (u <  inf ); 
arg2 = (g >= 0) & (l > -inf);
arg3 = (g < 0)  & (u == inf); 
arg4 = (g >= 0) & (l == -inf);
v(arg1)  = (x(arg1) - u(arg1)); 
dv(arg1) = 1;
v(arg2)  = (x(arg2) - l(arg2)); 
dv(arg2) = 1;
v(arg3)  = -1;
dv(arg3) = 0;
v(arg4)  = 1;
dv(arg4) = 0;

