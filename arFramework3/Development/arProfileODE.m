% Calculation of the profile likelihood via integration method

function [xsout, chi2sout] = arProfileODE(jk, gamma, p_range, method)

global ar
global pleGlobals

if(~exist('gamma','var'))
    gamma = 1;
end

if(~exist('method','var'))
    method = 'ode23';
%     method = 'euler';
%     method = 'ode15s';
%     method = 'ode45';
%     method = 'ode113';
end

n_reopt = 1;
% quadprog_optims = optimset('Display', 'off');

dchi2 = chi2inv(0.95, 1);
dchi2_increase = 0.1;
step_factor = 1.2;
step_factor_opt = 10;

doDAE = false;

if(strcmp(method,'euler'))
    fhandel = @euler;
    doDAE = false;
elseif(strcmp(method,'ode15s'))
    fhandel = @ode15s;
elseif(strcmp(method,'ode45'))
    fhandel = @ode45;
    doDAE = false;
elseif(strcmp(method,'ode113'))
    fhandel = @ode113;
    doDAE = false;
elseif(strcmp(method,'ode23'))
    fhandel = @ode23;
    doDAE = false;
end

p = ar.p;
arCalcMerit(true);

if(~isfield(ar,'H0'))
    H0 = ar.sres(:,ar.qFit==1)'*ar.sres(:,ar.qFit==1);
else
    H0 = ar.H0;
end

% initial stepsize
% H = ar.sres(:,ar.qFit==1)' * ar.sres(:,ar.qFit==1);
% Cvar = inv(H + eye(size(H)));
% dp = Cvar(jk,jk) * 1;
dp = 1;

jkx = zeros(size(p));
jkx(jk) = 1;
jkx = jkx(ar.qFit==1);
jkx = find(jkx);

x = [p(ar.qFit==1) 0];

if(isempty(pleGlobals.chi2s{jk}))
    p_rangeu = p_range;
    p_ranged = -p_range;
else
    p_rangeu = max(pleGlobals.ps{jk}(:,jk)) - pleGlobals.p(jk);
    p_ranged = min(pleGlobals.ps{jk}(:,jk)) - pleGlobals.p(jk);
end

if(doDAE)
    dir_up = true;
    options = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'OutputFcn', @my_odeplot, ...
        'Mass', @profile_ode_dae_lhs, 'MStateDependence', 'strong', ...
        'InitialSlope', profile_ode_dae_lhs(0,x)\profile_ode_dae_rhs(0,x));
    
    dir_up = false;
    optionsd = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'OutputFcn', @my_odeplot, ...
        'Mass', @profile_ode_dae_lhs, 'MStateDependence', 'strong', ...
        'InitialSlope', profile_ode_dae_lhs(0,x)\profile_ode_dae_rhs(0,x));
    
    frhs = @profile_ode_dae_rhs;
else
%     options = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'OutputFcn', @my_odeplot);
%     optionsd = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'OutputFcn', @my_odeplot);
    options = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'OutputFcn', @my_odeplot, 'MaxStep', 1e-2);
    optionsd = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'OutputFcn', @my_odeplot, 'MaxStep', 1e-2);
    frhs = @profile_ode;
end
% options = odeset;
try
    tic;
    % up
    fprintf('going up\n');
    xs2 = ones(0,length(p));
    chi2s = [];
    lambdas = [];
    
    firstIter = true;
    dir_up = true;
    if(p_rangeu>0)
        xs = [];
        try 
            feval(fhandel, frhs, [0 p_rangeu], x, options);
        catch err_id
            disp(err_id.message);
        end
        
        if(~isempty(xs))
            xs2 = ones(size(xs,1),1) * p;
            xs2(:,ar.qFit==1) = xs(:,1:(end-1));
            lambdas = xs(:,end);
        else
            xs2 = [];
        end
        
        chi2s = zeros(1,size(xs2,1));
        for j=1:size(xs2,1)
            ar.p = xs2(j,:);
            try
                arCalcMerit(true);
                chi2s(j) = ar.chi2fit;
            catch err_id 
                chi2s(j) = nan;
            end
        end
    end
    
    ar.p = p;
    arCalcMerit(true);
    
    % down
    fprintf('going down\n');
    xs2d = ones(0,length(p));
    chi2sd = [];
    lambdasd = [];
    
    firstIter = true;
    dir_up = false;
    if(p_ranged<0)
        xs = [];
        try 
            feval(fhandel, frhs, [0 abs(p_ranged)], x, optionsd);
        catch err_id
            disp(err_id.message);
        end
        
        if(~isempty(xs))
            xs2d = ones(size(xs,1),1) * p;
            xs2d(:,ar.qFit==1) = xs(:,1:(end-1));
            lambdasd = xs(:,end);
        else
            xs2d = [];
        end
        
        chi2sd = zeros(1,size(xs2d,1));
        for j=1:size(xs2d,1)
            ar.p = xs2d(j,:);
            try
                arCalcMerit(true);
                chi2sd(j) = ar.chi2fit;
            catch err_id 
                chi2sd(j) = nan;
            end
        end
    end
    
    ar.p = p;
    arCalcMerit(true);
    
    xsout = [flipud(xs2d); xs2];
    chi2sout = [fliplr(chi2sd) chi2s];
    lambdasout = [flipud(lambdasd); lambdas];
    
    tt = toc;
    fprintf('runtime %f sec (x%f faster than conventional)\n', tt, pleGlobals.timing(jk)/tt);
    
    if(~isempty(xsout))
        figure(1)
        
        subplot(5,1,[1 2])
        plot(pleGlobals.ps{jk}(:,jk), pleGlobals.chi2s{jk}, 'or--');
        hold on
        plot(xsout(:,jk), chi2sout, 'x-b')
        legend('reoptimization', 'profile ODE');
        xlim([min(pleGlobals.ps{jk}(:,jk)) max(pleGlobals.ps{jk}(:,jk))])
        ylim([ar.chi2fit-(0.1*dchi2) ar.chi2fit+2*dchi2]);
        plot([ar.p(jk) ar.p(jk)], ylim, 'k--');
        plot(xlim, [ar.chi2fit+dchi2 ar.chi2fit+dchi2], 'r--');
        hold off
        
        subplot(5,1,[3 4])
        q_not_jk = (1:size(xsout,2))~=jk;
        plot(pleGlobals.ps{jk}(:,jk), bsxfun(@minus, pleGlobals.ps{jk}(:,q_not_jk), pleGlobals.p(q_not_jk)), 'o--');
        ylimtmp = ylim;
        hold on
        plot([ar.p(jk) ar.p(jk)], ylim, 'k--');
        plot(xsout(:,jk), bsxfun(@minus, xsout(:,q_not_jk), ar.p(q_not_jk)), 'x-')
        %     plot(xsout(:,jk), xsout(:,q_not_jk), 'x-')
        xlim([min(pleGlobals.ps{jk}(:,jk)) max(pleGlobals.ps{jk}(:,jk))])
        hold off
        legend(strrep(ar.pLabel(q_not_jk),'_','\_'))
        ylim(ylimtmp);
        
        subplot(5,1,5)
        plot(xsout(:,jk), lambdasout, 'x-')
        xlim([min(pleGlobals.ps{jk}(:,jk)) max(pleGlobals.ps{jk}(:,jk))])
        if(max(lambdasout)-min(lambdasout)>0)
            ylim([min(lambdasout) max(lambdasout)])
        end
        hold on
        plot([ar.p(jk) ar.p(jk)], ylim, 'k--');
        hold off
    end
catch err_id
    ar.p = p;
    arCalcMerit(true);
    rethrow(err_id);
end

    function euler(rhs_fun, trange, xinit, ~)
        xs = xinit;
        
        dp2 = 1;
        dp_tmp = dp;
        arCalcMerit(true);
        
        arWaitbar(0);
        
        mindp = 1e-3;
        nmax = 10;
        
        inmax = 1;
        try
            while((dir_up && xinit(jk)-x(jk) < trange(2)) || ...
                    (~dir_up && x(jk)-xinit(jk) < trange(2)))
                arWaitbar(ceil((xinit(jk)-x(jk))/trange(2)*100),100);
                
                chi2_tmp = ar.chi2fit;
                rhs_tmp = feval(rhs_fun, xinit(jkx), xinit, false);
                
                xtrial = xinit + transpose(rhs_tmp*dp_tmp);
                arCalcMerit(true, xtrial(1:(end-1)));
                chi2_trial = ar.chi2fit;
                if(chi2_trial - chi2_tmp > dchi2*dchi2_increase) % reduce step
                    while(chi2_trial - chi2_tmp > dchi2*dchi2_increase)
                        dp_tmp = dp_tmp / step_factor;
                        if(dp_tmp < mindp)
                            break;
                        end
                        xtrial = xinit + transpose(rhs_tmp*dp_tmp);
                        arCalcMerit(true, xtrial(1:(end-1)));
                        chi2_trial = ar.chi2fit;
                    end
                elseif(chi2_trial - chi2_tmp < dchi2*0.5) % increase next step
                    dp_tmp = dp_tmp * step_factor;
                end
                xinit = xtrial;
                
                % reoptimize
                [xinit, dp2] = reopt(xinit, dp2);
                
                xs = [xs;xinit]; %#ok<AGROW>
                
                inmax = inmax + 1;
                if(inmax > nmax)
                    break;
                end
            end
            arWaitbar(-1);
        catch err_id
            arWaitbar(-1);
            rethrow(err_id)
        end
    end

    function [xxstmp, dp2] = reopt(xxstmp, dp2)
        if(n_reopt>0)
            arCalcMerit(true, xxstmp(1:(end-1)));
            
            qFitReset = ar.qFit;
            ar.qFit(jk) = 0;
            qFit = ar.qFit==1;
            
            iter = 1;
            refresh = 1;
            pBest = ar.p(qFit);
%             lb = ar.lb(qFit);
%             ub = ar.ub(qFit);
            
            maxiter = 1;
            try
                while iter<=n_reopt
                    maxiter = maxiter + 1;
                    if(maxiter>100)
                        break;
                    end
                    if(refresh==1)
                        % backup old status
                        p_old = ar.p(qFit) + 0;
                        chi2_old = ar.chi2fit + 0;
                        
                        res = ar.res;
                        sres = ar.sres(:,qFit);
                    end
                    
                    % solve sub problem
                    deltap = -pinv(sres)*res';
                    
                    % set new candidate parameter values
                    ar.p(qFit) = p_old + dp2*deltap';
                    
%                     % sub problem search space
%                     lbtmp = lb-p_old;
%                     ubtmp = ub-p_old;
%                     lbtmp(lbtmp < -dp2) = -dp2;
%                     ubtmp(ubtmp > dp2) = dp2;
% 
%                     % solve sub problem
%                     deltap = lsqlin(-sres, res, [], [], [], [], ...
%                         lbtmp, ubtmp, [], quadprog_optims);
%                     
%                     % set new candidate parameter values
%                     ar.p(qFit) = p_old + deltap';
                    
                    % evaluate objective function
                    arCalcMerit(true);
                    
                    % compare old and new chi^2
                    if(ar.chi2fit < chi2_old) % if decreased, keep the step
                        iter = iter + 1;
                        
                        refresh = 1;
                        
                        % increase step size
                        dp2 = dp2 * step_factor_opt;
                        
                        pBest = ar.p(qFit);
                    else % if increased, reset & do it again with smaller step
                        refresh = 0;
                        
                        % decrease step size
                        dp2 = dp2 / step_factor_opt;
                        
                        % backup
                        ar.p(qFit) = p_old;
                    end
                end
            catch err_id
                disp(err_id.message)
            end
            ar.p(qFit) = pBest;
            ar.qFit = qFitReset;
            
            xxstmp = [ar.p(ar.qFit==1) xxstmp(end)];
        end
    end

    function dy = profile_ode(~,y,reeval)        
        if(nargin < 3)
            reeval = true;
        end
        
        if(reeval)
            arCalcMerit(true, y(1:(end-1)));
        end
        
        gdot = zeros(sum(ar.qFit==1),1);
        gdot(jkx) = 1;
        
        if(~firstIter)
            H = 2*ar.sres(:,ar.qFit==1)'*ar.sres(:,ar.qFit==1);
        else
            H = H0;
            firstIter = false;
        end
%         H = eye(size(H));
%         pReset = ar.p;
%         H = hessian;
%         ar.p = pReset;
        
        M = [H gdot; gdot' 0];
        
        g = transpose(-ar.res*ar.sres(:,ar.qFit==1));
        if(~dir_up)
            v1 = [-gamma*g; 1];
        else
            v1 = [gamma*g; 1];
        end
        v2 = [zeros(sum(ar.qFit==1),1); gamma*y(end)];
        dy = pinv(M)*v1 - v2;
        
%         v = [zeros(sum(ar.qFit==1),1); 1];
%         dy = pinv(M)*v;
        
%         [U,S,V] = svd(M);
%         s = diag(S);
%         qs = s/max(s) > 1e-3;
%         Snew = zeros(size(S));
%         for jj=1:sum(qs);
%             Snew(jj,jj) = 1/s(jj);
%         end
%         
%         Minv = V * Snew * U';
%         dy = Minv*v1 - v2;        
        
%         dy = M\v1 - v2;
%         dy = inv(M)*v1 - v2;
        
        if(~dir_up)
            dy = -dy;
        end
        
        % check constraint equation (2.1)
        % disp(norm(g + y(end)*gdot))
    end

    

    function M = profile_ode_dae_lhs(~,y)        
        arCalcMerit(true, y(1:(end-1)));
        
        gdot = zeros(sum(ar.qFit==1),1);
        gdot(jkx) = 1;
        
        H = ar.sres(:,ar.qFit==1)'*ar.sres(:,ar.qFit==1);
        
        M = [H gdot; gdot' 0];
    end


    function dy = profile_ode_dae_rhs(~,y)        
        arCalcMerit(true, y(1:(end-1)));
        
        gdot = zeros(sum(ar.qFit==1),1);
        gdot(jkx) = 1;
        
        g = transpose(-ar.res*ar.sres(:,ar.qFit==1));
        
        v1 = [-gamma*g; 1];
        v2 = [gamma*gdot*y(end); 0];
        
        dy = v1 - v2;
        if(~dir_up)
            dy = -dy;
        end
    end


    function status = my_odeplot(~,y,~,~)
        if(isempty(xs))
            xs = y';
        else
            xs = [xs;y'];
        end
        status = 0;
    end
end

