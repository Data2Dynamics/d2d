% adaptive step reconciling chi^2 increase
% linear step method

function [pStep, dpNew, beta, alpha] = pleInitStepLinear(jk, pLast, dpLast)

global pleGlobals;

% mag factor
stepfaktor = 2;

dchi2 = pleGlobals.dchi2_point;

chi2last = feval(pleGlobals.merit_fkt);

dpNew = dpLast;
if(length(dpLast) ~= length(pLast))
    dpNewtmp = zeros(size(pLast));
    dpNewtmp(jk) = dpLast;
    dpNew = dpNewtmp;
end
pStep = dpNew;

beta = [];
alpha = [];

q_not_jk = 1:length(pLast);
q_not_jk = q_not_jk~=jk & pleGlobals.q_fit;
while(true)
    q_hit_lb = pLast+dpNew <= pleGlobals.lb+pleGlobals.minstepsize;
    q_hit_ub = pLast+dpNew >= pleGlobals.ub-pleGlobals.minstepsize;
    
    if(q_hit_lb(jk) || q_hit_ub(jk)) % jk hit boundaries
        dpNew = dpNew / stepfaktor;
        if(abs(dpNew(jk))<pleGlobals.minstepsize(jk))
            if(dpLast(jk)>0)
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
    else
        feval(pleGlobals.integrate_fkt, pLast+dpNew);
        chi2trial = feval(pleGlobals.merit_fkt);
        
        if(chi2trial-chi2last > dchi2*pleGlobals.relchi2stepincrease(jk))
            dpNew = dpNew / stepfaktor;
            if(abs(dpNew(jk))<pleGlobals.minstepsize(jk))
                fprintf('WARNING: could not control step size (minstepsize = %e)\n', pleGlobals.minstepsize(jk))
                dpNew = dpNew * stepfaktor;
                return
            end
        else
            if(chi2trial-chi2last < dchi2*pleGlobals.relchi2stepincrease(jk) && ...
                    chi2trial-chi2last>0)
                dpNew = dpNew * stepfaktor;
                if(abs(dpNew(jk)) > pleGlobals.maxstepsize(jk))
                    dpNew = dpNew / abs(dpNew(jk)) * pleGlobals.maxstepsize(jk);
                end
            end
            pStep = dpNew;
            return
        end
    end
    pStep = dpNew;
end

