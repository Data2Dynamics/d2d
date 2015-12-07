% adaptive step reconciling chi^2 increase
% direct step method
% 
%   The new functionality enters after five steps.
%   The idea is to use previous chi2-increase to predict how many steps are
%   necessary to cross the threshold. If this exceeds the samplesize, the
%   stepsize is increase by a factor sqrt(2).
% 
%   last.dx      previous stepsize      [not used here, but provided for
%                                           alternative strategies]
%   last.dy   previous chi2-change

function [pStep, dpNew, beta, alpha] = pleInitStepDirect2(jk, pLast, dpLast, last)

global pleGlobals;

% mag factor
stepfaktor = 2;

dchi2 = pleGlobals.dchi2_point;

chi2last = feval(pleGlobals.merit_fkt);

dpNew = dpLast;
pStep = zeros(size(pLast));
pStep(jk) = dpNew;

beta = [];
alpha = [];

q_not_jk = 1:length(pLast);
q_not_jk = q_not_jk~=jk & pleGlobals.q_fit;
while(true)
    q_hit_lb = pLast+pStep <= pleGlobals.lb+pleGlobals.minstepsize;
    q_hit_ub = pLast+pStep >= pleGlobals.ub-pleGlobals.minstepsize;
    
    if(q_hit_lb(jk) || q_hit_ub(jk)) % jk hit boundaries
        dpNew = dpNew / stepfaktor;
        if(abs(dpNew)<pleGlobals.minstepsize(jk))
            if(dpLast>0)
                lbub = 'upper';
            else
                lbub = 'lower';
            end
            fprintf('PLE#%i parameter %s hit %s boundary\n', jk, pleGlobals.p_labels{jk}, lbub);
            pStep = nan(size(pLast));
            return
        end
    elseif(sum(q_hit_lb(q_not_jk) & pleGlobals.breakonlb(q_not_jk))>0 ||...
            sum(q_hit_ub(q_not_jk) & pleGlobals.breakonub(q_not_jk))>0) % other than jk hit boundaries
        fprintf('STOP: other parameters hit hard boundaries:\n')
        fprintf('\t%s\n', pleGlobals.p_labels{(q_hit_lb | q_hit_ub) & q_not_jk})
        pStep = nan(size(pLast));
        return
    elseif sum(~isnan(last.dx))>5  % if >5 last values are available
        ub = pleGlobals.ub(jk);
        lb = pleGlobals.lb(jk);
        ytarget = icdf('chi2',0.95,1)*1.2;
        minx = pleGlobals.minstepsize(jk);
        maxx = pleGlobals.maxstepsize(jk);
        ss = pleGlobals.samplesize(jk);

        dpNew = profileStepControl(last,lb,ub,ytarget,minx,maxx,ss);

        pStep(jk) = dpNew;
        return
            
    else
        feval(pleGlobals.integrate_fkt, pLast+pStep);
        chi2trial = feval(pleGlobals.merit_fkt);
        
        if(chi2trial-chi2last > dchi2*pleGlobals.relchi2stepincrease(jk))
            dpNew = dpNew / stepfaktor;
            if(abs(dpNew)<pleGlobals.minstepsize(jk))
%                 fprintf('WARNING: could not control step size (minstepsize = %e)\n', pleGlobals.minstepsize(jk))
                dpNew = dpNew * stepfaktor;
                return
            end
        else
            if(chi2trial-chi2last < dchi2*pleGlobals.relchi2stepincrease(jk) && ...
                    chi2trial-chi2last>0)
                dpNew = dpNew * stepfaktor;
                if(abs(dpNew) > pleGlobals.maxstepsize(jk))
                    dpNew = pleGlobals.maxstepsize(jk) * sign(dpNew);
                end
            end
            return
        end
    end

    pStep(jk) = dpNew;
end

