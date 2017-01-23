% adaptive step reconciling chi^2 increase
% linear step method

function [pStep, dpNew, beta, alpha] = pleInitStepLinear(jk, pLast, dpLast, lastAll)

global ar

% mag factor
stepfaktor = 2;

dchi2 = ar.ple.dchi2_point;

chi2last = feval(ar.ple.merit_fkt);

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
q_not_jk = q_not_jk~=jk & ar.qFit;
while(true)
    q_hit_lb = pLast+dpNew <= ar.lb+ar.ple.minstepsize;
    q_hit_ub = pLast+dpNew >= ar.ub-ar.ple.minstepsize;
    
    if(q_hit_lb(jk) || q_hit_ub(jk)) % jk hit boundaries
        dpNew = dpNew / stepfaktor;
        if(abs(dpNew(jk))<ar.ple.minstepsize(jk))
            if(dpLast(jk)>0)
                lbub = 'upper';
            else
                lbub = 'lower';
            end
            fprintf('PLE#%i parameter %s hit %s boundary\n', jk, ar.ple.p_labels{jk}, lbub);
            pStep = nan(size(pLast));
            return
        end
    elseif(sum(q_hit_lb(q_not_jk) & ar.ple.breakonlb(q_not_jk))>0 ||...
            sum(q_hit_ub(q_not_jk) & ar.ple.breakonub(q_not_jk))>0) % other than jk hit boundaries
        fprintf('STOP: other parameters hit hard boundaries:\n')
        fprintf('\t%s\n', ar.ple.p_labels{(q_hit_lb | q_hit_ub) & q_not_jk})
        pStep = nan(size(pLast));
        return
    else
        feval(ar.ple.integrate_fkt, pLast+dpNew);
        chi2trial = feval(ar.ple.merit_fkt);
        
        if(chi2trial-chi2last > dchi2*ar.ple.relchi2stepincrease(jk))
            dpNew = dpNew / stepfaktor;
            if(abs(dpNew(jk))<ar.ple.minstepsize(jk))
                fprintf('WARNING: could not control step size (minstepsize = %e)\n', ar.ple.minstepsize(jk))
                dpNew = dpNew * stepfaktor;
                return
            end
        else
            if(chi2trial-chi2last < dchi2*ar.ple.relchi2stepincrease(jk) && ...
                    chi2trial-chi2last>0)
                dpNew = dpNew * stepfaktor;
                if(abs(dpNew(jk)) > ar.ple.maxstepsize(jk))
                    dpNew = dpNew / abs(dpNew(jk)) * ar.ple.maxstepsize(jk);
                end
            end
            pStep = dpNew;
            return
        end
    end
    pStep = dpNew;
end

