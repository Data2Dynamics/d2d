% generate trial steps within box bounds
%
% [dp, solver_calls, qred, dpmem, grad_dir_frac, llh_expect, normdpmu_type] = ...
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
% 15 = trdog pcgr 2D subspace 
% 16 = trdog pcgr (no DM) 2D subspace 
% 17 = trdog pcgr (no DM) 3D subspace 
% 18 = trust region (with DG, based on modified trust.m)
% 19 = trust region solution (with EV)

function [dp, solver_calls, qred, grad_dir_frac, llh_expect, normdpmu_type] = ...
    arNLSstep(llh, g, H, sres, mu, p, lb, ub, solver_calls, dpmem, method)

if(nargin==0)
    dp = {};
    dp{1} = 'trustregion';
    dp{2} = 'Levenberg-Marquardt';
    dp{3} = 'Newton';
    dp{4} = 'gradient-descent';
    dp{5} = 'gradient-descent';
    dp{6} = 'dogleg';
    dp{7} = 'generalized-trustregion';
    dp{8} = 'MATLAB-trdog';
    dp{9} = 'Newton-pcgr';
    dp{10} = 'trdog-pcgr';
    dp{11} = 'dogleg-Newton-pcgr';
    dp{12} = 'dogleg-trdog-pcgr';
    dp{13} = 'trdog-pcgr-noDM';
    dp{14} = 'trdog-pcgr-noDG';
    dp{15} = 'trdog-pcgr-Levenberg-Marquardt';
    dp{16} = 'trdog-pcgr-2D-subspace';
    dp{17} = 'trdog-pcgr-noDM-2D-subspace';
    dp{18} = 'trdog-pcgr-noDM-3D-subspace';
    dp{19} = 'trustregion-withDG';
    dp{20} = 'trustregion-withEV';
    return;
end

qred = false(size(p));

if(isscalar(mu) && mu==0)
    dp = zeros(size(g));
    grad_dir_frac = 0;
    llh_expect = llh;
    normdpmu_type = nan;
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
    normdpmu_type = nan;
    
else % own step implementation
    % solve subproblem
    [dp, normdpmu_type, solver_calls_tmp] = getDP(p, lb, ub, g, H, sres, mu, dpmem, method, qred);
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
        [dp_red, normdpmu_type, solver_calls_tmp] = getDP(p(~qred), lb(~qred), ub(~qred), ...
            gred, Hred, sresred, mu, dpmem, method, qred);
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

% proportion of step in direction of gradient
grad_dir_frac = dp(~qred)*g(~qred)'/norm(dp(~qred))/norm(g(~qred));

% expected likelihood
llh_expect = llh - g*dp' + 0.5*dp*H*dp';



% generate trial steps within ball of size mu (for LM: lambda = 10/mu)
% normdpmu = 
%       norm(dp) in units of mu
%       nan if trust region scaling based on normdpmu not applicable
function [dp, normdpmu_type, solver_calls_tmp] = getDP(p, lb, ub, g, H, sres, mu, dpmem, method, qred)

if(mu==0)
    dp = zeros(size(g));
    normdpmu_type = nan;
    solver_calls_tmp = 0;
    return;
else
    solver_calls_tmp = 1;
end

if(~isempty(dpmem))
    dpmem = dpmem(~qred);
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
        normdpmu_type = -1;
        
    case 1 % levenberg-marquardt
        lambda = 1/mu;
        dp = transpose(pinv(H + lambda*eye(size(H)))*g');
        normdpmu_type = nan;
        
    case 2 % newton with maximum steplength mu
        dp = transpose(pinv(H)*g');
        if(norm(dp)>mu)
            dp = mu*dp/norm(dp);
        end
        normdpmu_type = -1;
        
    case 3 % gradient with steplength mu
        dp = mu*g/norm(g);
        normdpmu_type = -1;
       
    case 4 % gradient to cauchy point with steplength mu
        dp = norm(g)^2/(g*H*g') * g;
        if(norm(dp)>mu)
            dp = mu*dp/norm(dp);
        end
        normdpmu_type = -1;
        
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
        normdpmu_type = -1;
        
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
        normdpmu_type = norm(dp) / sqrt((dp/mut)*(dp/mut)');
        
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
        normdpmu_type = -1;
        
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
        normdpmu_type = -1;
   
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
        normdpmu_type = -1;
        
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
        normdpmu_type = -1;
        
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
        normdpmu_type = -1;
        
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
        normdpmu_type = -1;
        
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
        
        normdpmu_type = nan;
        
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
        normdpmu_type = norm(dp)/mu;
        
        % the function trust has a bug, therefore, in some cases
        % the problem has to be regularized
        lambda = 1e-6;
        while(norm(dp)/mu > 1.01)
%             fprintf('trust.m problem %g, regularizing with new lambda=%g\n', norm(dp)/mu, lambda);
            dp = trust(rhs,MM+lambda*eye(size(MM)),mu);
            normdpmu_type = norm(dp)/mu;
            solver_calls_tmp = solver_calls_tmp + 1;
            lambda = lambda * 10;
        end
        
        % got back to original space
        dp = Z*dp;
        dp = transpose(DM*dp);
        
    case 16 % trdog pcgr no DM 2D subspace
        g_pcgr = -0.5*g';
        pcoptions = Inf;
        pcgtol = 0.1;
        kmax = max(1,floor(length(g)/2));
        [v, dv] = definev(g_pcgr',p,lb,ub);
        dd = sqrt(abs(v));
        DM = sparse(diag(dd));
        DM = sparse(eye(size(DM)));
        DG = sparse(1:length(g),1:length(g),full(abs(g_pcgr).*dv));
%         DG = sparse(1:length(g),1:length(g),full(ones(size(g_pcgr)).*dv));
        grad = DM*g_pcgr;
        [R,permR] = feval(@aprecon,sres,pcoptions,DM,DG);
        [dp, posdef] = pcgr(DM,DG,grad,kmax,pcgtol,@atamult,sres,R,permR,'jacobprecon',pcoptions);

%         posdef = 1;
%         dp = (DM*H*DM + DG)\g';
        
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
        normdpmu_type = -1;
        
    case 17 % trdog pcgr (no DM) 3D subspace 
        g_pcgr = -0.5*g';
        pcoptions = Inf;
        pcgtol = 0.1;
        kmax = max(1,floor(length(g)/2));
        [v, dv] = definev(g_pcgr',p,lb,ub);
        dd = sqrt(abs(v));
        DM = sparse(diag(dd));
        DM = sparse(eye(size(DM)));
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
        
        % additional directions
%         tmpp = randn(1,length(v1));
%         Z(:,end+1) = tmpp/norm(tmpp);
        if(~isempty(dpmem))
            tmpp = DM*dpmem';
            Z(:,end+1) = tmpp/norm(tmpp);
        end
%         tmpp = transpose(pinv(H)*g');
%         Z(:,end+1) = tmpp/norm(tmpp);
        Z = mgrscho(Z);
        
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
        normdpmu_type = -1;
        
    case 18 % trust region solution (with DG)
        g_pcgr = -0.5*g';
        [~, dv] = definev(g_pcgr',p,lb,ub);
        DG = sparse(1:length(g),1:length(g),full(abs(g_pcgr).*dv));
        
        dp = trust(-g',H+DG,mu)';
        
        % the function trust has a bug, therefore, in some cases
        % the problem has to be regularized
        lambda = 1e-6;
        while(norm(dp)/mu > 1.01)
            %             fprintf('trust.m problem %g, regularizing with new lambda=%g\n', norm(dp)/mu, lambda);
            dp = trust(-g',H+DG+lambda*eye(size(H)),mu)';
            solver_calls_tmp = solver_calls_tmp + 1;
            lambda = lambda * 10;
        end
        normdpmu_type = -1;
        
    case 19 % trust region solution (with EV)
%         [U,S,V] = svd(H);
%         Q = S/max(S(:));
%         S2 = S;
%         S2(Q<1e-6) = 0;
%         H2 = U*S2*V';

        [V,D] = eig(H);
        q = diag(D/max(D(:))) > 0;
        
        Z = V(:,q);
        g = transpose(Z'*g');
        H = Z'*H*Z;

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
        
        dp = transpose(Z*dp');
        normdpmu_type = -1;
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



function A = mgrscho(A)
%MGRSHO Modified Gram-Schmidt orthogonalization procedure. 
% -For a basis of fundamentals on classical Gram-Schmidt process, procedure
% and its origin. Please see the text of the m-file cgrsho you can download
% from www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=12465
% The classical Gram-Scmidt algorithm is numerically unstable, mainly 
% because of all the successive subtractions in the order they appear. When
% this process is implemented on a computer, then the vectors s_n are not
% quite orthogonal because of rounding errors. This loss of orthogonality
% is particularly bad; therefore, it is said that the (naive) classical 
% Gram–Schmidt process is numerically unstable. If we write an algorithm
% based on the way we developed the Gram-Schmidt iteration (in terms of 
% projections), we get a better algorithm.
% The Gram–Schmidt process can be stabilized by a small modification. 
% Instead of computing the vector u_n as,
%
%     u_n = v_k - proj_u_1 v_n - proj_u_2 v_n -...- proj_u_n-1 v_n
%
% it is computed as,
%
%     u_n = u_n ^n-2 - proj_u_n-1 u_n ^n-2
%
% This series of computations gives the same result as the original formula
% in exact arithmetic, but it introduces smaller errors in finite-precision
% arithmetic. A stable algorithm is one which does not suffer drastically 
% from perturbations due to roundoff errors. This is called as the modified
% Gram-Schmidt orthogonalization process. 
% There are several different variations of the Gram-Schmidt process 
% including classical Gram-Schmidt (CGS), modified Gram-Schmidt (MGS) and 
% modified Gram-Schmidt with pivoting (MGSP). MGS economizes storage and is
% generally more stable than CGS.
% The Gram-Schmidt process can be used in calculating Legendre polynomials,
% Chebyshev polynomials, curve fitting of empirical data, smoothing, and
% calculating least square methods and other functional equations.
% 
% Syntax: function mgrscho(A)
%
% Input:
%    A - matrix of n linearly independent vectors of equal size. Here, them
%        must be arranged as columns.
% Output:
%    Matrix of n orthogonalized vectors.
%
% Example: Taken the problem 18, S6.3, p308, from the Mathematics 206 Solutions
% for HWK 24b. Course of Math 206 Linear Algebra by Prof. Alexia Sontag at
% Wellesley Collage, Wellesley, MA, USA. URL address:
% http://www.wellesley.edu/Math/Webpage%20Math/Old%20Math%20Site/Math206sontag/
% Homework/Pdf/hwk24b_s02_solns.pdf
%           
% We are interested to orthogonalize the vectors,
%
%    v1 = [0 2 1 0], v2 = [1 -1 0 0], v3 = [1 2 0 -1] and v4 = [1 0 0 1]
%
% by the modified Gram-Schmidt process.
%
% Vector matrix must be:
%    A = [0 1 1 1;2 -1 2 0;1 0 0 0;0 0 -1 1];
%
% Calling on Matlab the function: 
%    mgrscho(A)
%
% Answer is:
%
% ans =
%         0    0.9129    0.3162    0.2582
%    0.8944   -0.1826    0.3162    0.2582
%    0.4472    0.3651   -0.6325   -0.5164
%         0         0   -0.6325    0.7746
%
% NOTE.- Comparing the orthogonality of resulting vectors by both classical
%    Gram-Schmidt and modified Gram-Schmidt processes, using floating point
%    format with 15 digits for double and 7 digits for single. We found that
%    during the process, with the modified one there exists smaller errors. 
%
%                                  Gram-Schmidt Process
%                  --------------------------------------------------------
%        Vectors         Classical                         Modified
%       -------------------------------------------------------------------
%         A1-A2   -8.326672684688674e-017          -8.326672684688674e-017
%         A1-A3    1.665334536937735e-016           1.110223024625157e-016
%         A1-A4    1.110223024625157e-016           1.110223024625157e-016
%         A2-A3    5.551115123125783e-017          -1.110223024625157e-016
%         A2-A4   -5.551115123125783e-017          -5.551115123125783e-017
%         A3-A4              0                                 0
%       -------------------------------------------------------------------
%
%
% Created by A. Trujillo-Ortiz, R. Hernandez-Walls, A. Castro-Perez
%            and K. Barba-Rojo
%            Facultad de Ciencias Marinas
%            Universidad Autonoma de Baja California
%            Apdo. Postal 453
%            Ensenada, Baja California
%            Mexico.
%            atrujo@uabc.mx
%
% Copyright. September 30, 2006.
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A., R. Hernandez-Walls, A. Castro-Perez and K. Barba-Rojo. (2006). 
%   mgrscho:Modified Gram-Schmidt orthogonalization procedure. A MATLAB file.
%   [WWW document]. URL http://www.mathworks.com/matlabcentral/fileexchange/
%   loadFile.do?objectId=12495
%
% References:
% Gerber, H. (1990), Elementary Linear Algebra. Brooks/Cole Pub. Co. Pacific
%     Grove, CA. 
% Wong, Y.K. (1935), An Application of Orthogonalization Process to the 
%     Theory of Least Squares. Annals of Mathematical Statistics, 6:53-75.
%

if nargin ~= 1,
    error('You need to imput only one argument.');
end

[~, n]=size(A);

for j= 1:n
    R(j,j)=norm(A(:,j)); %#ok<AGROW>
    A(:,j)=A(:,j)/R(j,j);
    R(j,j+1:n)=A(:,j)'*A(:,j+1:n); %#ok<AGROW>
    A(:,j+1:n)=A(:,j+1:n)-A(:,j)*R(j,j+1:n);
end
