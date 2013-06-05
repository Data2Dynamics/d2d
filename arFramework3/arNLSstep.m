% generate trial steps within box bounds
%
% [dp, solver_calls, qred, dpmem, grad_dir_frac, llh_expect, mudp] = ...
%     arNLSstep(llh, g, H, sres, mu, p, lb, ub, solver_calls, dpmem, useInertia, method)
%
% llh:  likelihood
% g:    gradient
% H:    Hessian matrix
% mu:   trust region size (for LM: lambda = 10/mu)
% p:    current parameters
% lb:   lower bounds
% ub:   upper bounds
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
% 15 = trdog pcgr 2D subspace (open: how to combine, how to implement updates)

function [dp, solver_calls, qred, dpmem, grad_dir_frac, llh_expect, normdpmu] = ...
    arNLSstep(llh, g, H, sres, mu, p, lb, ub, solver_calls, dpmem, useInertia, method)

if(nargin==0)
    dp = {};
    dp{1} = 'trustregion';
    dp{2} = 'Levenberg-Marquardt';
    dp{3} = 'Newton';
    dp{4} = 'gradient-descent';
    dp{5} = 'gradient-descent';
    dp{6} = 'dogleg';
    dp{7} = 'generalized-trustregion';
    dp{8} = 'trdog';
    dp{9} = 'Newton-pcgr';
    dp{10} = 'trdog-pcgr';
    dp{11} = 'dogleg-Newton-pcgr';
    dp{12} = 'dogleg-trdog-pcgr';
    dp{13} = 'trdog-pcgr-noDM';
    dp{14} = 'trdog-pcgr-noDG';
    dp{15} = 'trdog-pcgr-Levenberg-Marquardt';
    dp{16} = 'trdog-pcgr-2D-subspace';
    return;
end

qred = false(size(p));

if(isscalar(mu) && mu==0)
    dp = zeros(size(g));
    grad_dir_frac = 0;
    llh_expect = llh;
    normdpmu = 0;
    return
end

if(method==7) % use MATLABs trdog
    g_trdog = -0.5*g';
    [v, dv] = definev(g_trdog',p,lb,ub);
    dd = sqrt(abs(v)); 
    D = sparse(diag(dd));
    pcoptions = Inf;
    pcgtol = 0.1;
    kmax = max(1,floor(length(p)/2));
    gopt = v.*g_trdog;
    optnrm = norm(gopt,inf);
    theta = max(.95,1-optnrm);  
    
    dp = trdog(p',g_trdog,sres,D,mu,dv,@atamult,@aprecon,...
        pcoptions,pcgtol,kmax,theta,lb',ub',[],[],'jacobprecon')';
    
    solver_calls = solver_calls + 1;
    normdpmu = nan;
    
else % own step implementation
    % solve subproblem
    [dp, normdpmu, solver_calls_tmp] = getDP(p, lb, ub, g, H, sres, mu, method, qred);
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
        sresred = sres(:,~qred);
        
        % solve reduced subproblem
        [dp_red, normdpmu, solver_calls_tmp] = getDP(p(~qred), lb(~qred), ub(~qred), ...
            gred, Hred, sresred, mu, method, qred);
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
% normdpmu = 
%       norm(dp) in units of mu
%       nan if trust region scaling based on normdpmu not applicable
function [dp, normdpmu, solver_calls_tmp] = getDP(p, lb, ub, g, H, sres, mu, method, qred)

if(mu==0)
    dp = zeros(size(g));
    normdpmu = 0;
    solver_calls_tmp = 0;
    return;
else
    solver_calls_tmp = 1;
end

switch method
    case 0 % trust region solution
        dp = trust(-g',H,mu)';
        
        % the function trust has a bug, therefore, in some cases
        % the problem has to be regularized
        lambda = 1e-6;
        while(norm(dp)/mu > 1.01)
%             fprintf('trust.m problem %g, regularizing with new lambda=%g\n', norm(dp)/mu, lambda);
            dp = trust(-g',H+lambda*eye(size(H)),mu)';
            solver_calls_tmp = solver_calls_tmp + 1;
            lambda = lambda * 10;
        end
        normdpmu = norm(dp)/mu;
        
    case 1 % levenberg-marquardt
        lambda = 1/mu;
        dp = transpose(pinv(H + lambda*eye(size(H)))*g');
        normdpmu = nan;
        
    case 2 % newton with maximum steplength mu
        dp = transpose(pinv(H)*g');
        if(norm(dp)>mu)
            dp = mu*dp/norm(dp);
        end
        normdpmu = norm(dp)/mu;
        
    case 3 % gradient with steplength mu
        dp = mu*g/norm(g);
        normdpmu = 1;
       
    case 4 % gradient to cauchy point with steplength mu
        dp = norm(g)^2/(g*H*g') * g;
        if(norm(dp)>mu)
            dp = mu*dp/norm(dp);
        end
        normdpmu = norm(dp)/mu;
        
    case 5 % dogleg
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
        normdpmu = norm(dp)/mu;
        
    case 6 % generalized trust region
        mut = mu(~qred,~qred);
        gt = g*mut;
        Ht = mut'*H*mut;
        Ht = 0.5*(Ht + Ht');
        
        dp = trust(-gt',Ht,1)';
        
        % the function trust has a bug, therefore, in some cases
        % the problem has to be regularized
        lambda = 1e-6;
        while(norm(dp) > 1.01)
%             fprintf('trust.m problem %g, regularizing with new lambda=%g\n', norm(dp), lambda);
            dp = trust(-gt',Ht+lambda*eye(size(Ht)),1)';
            solver_calls_tmp = solver_calls_tmp + 1;
            lambda = lambda * 10;
        end
        dp = dp * mut;
        normdpmu = norm(dp) / sqrt((dp/mut)*(dp/mut)');
        
    case 8 % Newton pcgr
        g_pcgr = -0.5*g';
        pcoptions = Inf;
        pcgtol = 0.1;
        kmax = max(1,floor(length(g)/2));
        DM = eye(length(g));
        DG = DM;
        [R,permR] = feval(@aprecon,sres,pcoptions,DM,DG);
        dp = pcgr(DM,DG,g_pcgr,kmax,pcgtol,@atamult,sres,R,permR,'jacobprecon',pcoptions)';
        if(norm(dp)>mu)
            dp = mu*dp/norm(dp);
        end
        normdpmu = norm(dp)/mu;
        
    case 9 % trdog pcgr
        g_pcgr = -0.5*g';
        pcoptions = Inf;
        pcgtol = 0.1;
        kmax = max(1,floor(length(g)/2));
        [v, dv] = definev(g_pcgr',p,lb,ub);
        dd = sqrt(abs(v));
        DM = sparse(diag(dd));
%         DM = sparse(eye(size(DM)));
        DG = sparse(1:length(g),1:length(g),full(abs(g_pcgr).*dv));
%         DG = sparse(1:length(g),1:length(g),full(ones(size(g_pcgr)).*dv));
        [R,permR] = feval(@aprecon,sres,pcoptions,DM,DG);
        dp = pcgr(DM,DG,DM*g_pcgr,kmax,pcgtol,@atamult,sres,R,permR,'jacobprecon',pcoptions);
        dp = transpose(DM*dp);
        if(norm(dp)>mu)
            dp = mu*dp/norm(dp);
        end
        normdpmu = norm(dp)/mu;
   
    case 10 % dogleg Newton pcgr
        dpC = norm(g)^2/(g*H*g') * g;
        if(norm(dpC)>mu)
            % truncated gradient
            dp = mu*g/norm(g);
            fprintf('truncated gradient - ');
        else
            g_pcgr = -0.5*g';
            pcoptions = Inf;
            pcgtol = 0.1;
            kmax = max(1,floor(length(g)/2));
            DM = eye(length(g));
            DG = DM;
            [R,permR] = feval(@aprecon,sres,pcoptions,DM,DG);
            dpN = pcgr(DM,DG,g_pcgr,kmax,pcgtol,@atamult,sres,R,permR,'jacobprecon',pcoptions)';
            
            if(norm(dpN)<=mu)
                % full newton
                dp = dpN;
                fprintf('full newton - ');
            else
                % dogleg
                a = sum((dpN - dpC).^2);
                b = 2*sum(dpC.*(dpN - dpC));
                c = sum(dpC.^2) - mu^2;
                t = (-b + sqrt(b^2 - 4 * a * c))/(2*a);
                dp = dpC + t*(dpN - dpC);
                fprintf('dogleg (%4.2f) - ', norm(dpC)/mu);
            end
        end
        normdpmu = norm(dp)/mu;
        
    case 11 % dogleg trdog pcgr
        dpC = norm(g)^2/(g*H*g') * g;
        if(norm(dpC)>mu)
            % truncated gradient
            dp = mu*g/norm(g);
            fprintf('truncated gradient - ');
        else
            g_pcgr = -0.5*g';
            pcoptions = Inf;
            pcgtol = 0.1;
            kmax = max(1,floor(length(g)/2));
            [v, dv] = definev(g_pcgr',p,lb,ub);
            dd = sqrt(abs(v));
            DM = sparse(diag(dd));
%             DM = sparse(eye(size(DM)));
            DG = sparse(1:length(g),1:length(g),full(abs(g_pcgr).*dv));
%             DG = sparse(1:length(g),1:length(g),full(ones(size(g_pcgr)).*dv));
            [R,permR] = feval(@aprecon,sres,pcoptions,DM,DG);
            dpN = pcgr(DM,DG,DM*g_pcgr,kmax,pcgtol,@atamult,sres,R,permR,'jacobprecon',pcoptions);
            dpN = transpose(DM*dpN);
            
            if(norm(dpN)<=mu)
                % full newton
                dp = dpN;
                fprintf('full newton - ');
            else
                % dogleg
                a = sum((dpN - dpC).^2);
                b = 2*sum(dpC.*(dpN - dpC));
                c = sum(dpC.^2) - mu^2;
                t = (-b + sqrt(b^2 - 4 * a * c))/(2*a);
                dp = dpC + t*(dpN - dpC);
                fprintf('dogleg (%4.2f) - ', norm(dpC)/mu);
            end
        end
        normdpmu = norm(dp)/mu;
        
    case 12 % trdog pcgr no DM
        g_pcgr = -0.5*g';
        pcoptions = Inf;
        pcgtol = 0.1;
        kmax = max(1,floor(length(g)/2));
        [v, dv] = definev(g_pcgr',p,lb,ub);
        dd = sqrt(abs(v));
        DM = sparse(diag(dd));
        DM = sparse(eye(size(DM)));
        DG = sparse(1:length(g),1:length(g),full(abs(g_pcgr).*dv));
        [R,permR] = feval(@aprecon,sres,pcoptions,DM,DG);
        dp = pcgr(DM,DG,DM*g_pcgr,kmax,pcgtol,@atamult,sres,R,permR,'jacobprecon',pcoptions);
        dp = transpose(DM*dp);
        if(norm(dp)>mu)
            dp = mu*dp/norm(dp);
        end
        normdpmu = norm(dp)/mu;
        
    case 13 % trdog pcgr no DG
        g_pcgr = -0.5*g';
        pcoptions = Inf;
        pcgtol = 0.1;
        kmax = max(1,floor(length(g)/2));
        [v, dv] = definev(g_pcgr',p,lb,ub);
        dd = sqrt(abs(v));
        DM = sparse(diag(dd));
        DG = sparse(1:length(g),1:length(g),full(ones(size(g_pcgr)).*dv));
        [R,permR] = feval(@aprecon,sres,pcoptions,DM,DG);
        dp = pcgr(DM,DG,DM*g_pcgr,kmax,pcgtol,@atamult,sres,R,permR,'jacobprecon',pcoptions);
        dp = transpose(DM*dp);
        if(norm(dp)>mu)
            dp = mu*dp/norm(dp);
        end
        normdpmu = norm(dp)/mu;
        
    case 14 % trdog pcgr Levenberg-Marquardt
        lambda = 1/mu;
        sresLM = [sres;sqrt(lambda)*eye(size(sres,2))];
        
        g_pcgr = -0.5*g';
        pcoptions = Inf;
        pcgtol = 0.1;
        kmax = max(1,floor(length(g)/2));
        [v, dv] = definev(g_pcgr',p,lb,ub);
        dd = sqrt(abs(v));
        DM = sparse(diag(dd));
%         DM = sparse(eye(size(DM)));
        DG = sparse(1:length(g),1:length(g),full(abs(g_pcgr).*dv));
%         DG = sparse(1:length(g),1:length(g),full(ones(size(g_pcgr)).*dv));
        [R,permR] = feval(@aprecon,sresLM,pcoptions,DM,DG);
        dp = pcgr(DM,DG,DM*g_pcgr,kmax,pcgtol,@atamult,sresLM,R,permR,'jacobprecon',pcoptions);
        dp = transpose(DM*dp);
        
        normdpmu = nan;
        
    case 15 % trdog pcgr 2D subspace
        g_pcgr = -0.5*g';
        pcoptions = Inf;
        pcgtol = 0.1;
        kmax = max(1,floor(length(g)/2));
        [v, dv] = definev(g_pcgr',p,lb,ub);
        dd = sqrt(abs(v));
        DM = sparse(diag(dd));
%         DM = sparse(eye(size(DM)));
        DG = sparse(1:length(g),1:length(g),full(abs(g_pcgr).*dv));
%         DG = sparse(1:length(g),1:length(g),full(ones(size(g_pcgr)).*dv));
        grad = DM*g_pcgr;
        [R,permR] = feval(@aprecon,sres,pcoptions,DM,DG);
        [dp, posdef] = pcgr(DM,DG,grad,kmax,pcgtol,@atamult,sres,R,permR,'jacobprecon',pcoptions);
        
        % calculate subspace Z
        tol2 = sqrt(eps);
        if norm(dp) > 0
            v1 = dp/norm(dp);
        else
            v1 = dp;
        end
        Z(:,1) = v1;
        if length(p) > 1
            if (posdef < 1)
                v2 = D*sign(g_pcgr);
                if norm(v2) > 0
                    v2 = v2/norm(v2);
                end
            else
                if norm(g_pcgr) > 0
                    v2 = g_pcgr/norm(g_pcgr);
                else
                    v2 = g_pcgr;
                end
                
            end
            v2 = v2 - v1*(v1'*v2);
            nrmv2 = norm(v2);
            if nrmv2 > tol2
                v2 = v2/nrmv2;
                Z(:,2) = v2;
            end
        end
        
        % reduce to subspace
        W = DM*Z;        
        WW = feval(@atamult,sres,W,0);
        
        W = DM*WW;
        MM = full(Z'*W + Z'*DG*Z);
        rhs = full(Z'*grad);
        
        % determine 2D trust solution
        dp = trust(rhs,MM,mu);
        
        % the function trust has a bug, therefore, in some cases
        % the problem has to be regularized
        lambda = 1e-6;
        while(norm(dp)/mu > 1.01)
%             fprintf('trust.m problem %g, regularizing with new lambda=%g\n', norm(dp)/mu, lambda);
            dp = trust(rhs,MM+lambda*eye(size(MM)),mu);
            solver_calls_tmp = solver_calls_tmp + 1;
            lambda = lambda * 10;
        end
        
        % got back to original space
        dp = Z*dp;
        dp = transpose(DM*dp);
        
        normdpmu = norm(dp)/mu;
end




% Scaling vector and derivative
function [v,dv]= definev(g,x,l,u)
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

