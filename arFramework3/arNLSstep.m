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

function [dp, solver_calls, qred, dpmem, grad_dir_frac, llh_expect] = ...
    arNLSstep(llh, g, H, mu, p, lb, ub, solver_calls, dpmem, useInertia, method)
        
% solve subproblem
[dp, solver_calls_tmp] = getDP(g, H, mu, method);
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
    [dp_red, solver_calls_tmp] = getDP(gred, Hred, mu, method);
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


% generate trial steps within ball of size mu (for LM: lambda = 10/mu)
%
% method:   0 = trust region (based on modified trust.m)
%           1 = Levenberg-Marquardt
%           2 = Newton (with maximal step length mu)
%           3 = gradient descent (up to cauchy point)
%           4 = dogleg

function [dp, solver_calls_tmp] = getDP(g, H, mu, method)

solver_calls_tmp = 1;

if(mu==0)
    dp = zeros(size(g));
    return;
end

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
end

