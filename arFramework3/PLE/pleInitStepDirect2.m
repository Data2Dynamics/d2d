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

global ar

% mag factor
stepfaktor = 2;

dchi2 = ar.ple.dchi2_point;

chi2last = feval(ar.ple.merit_fkt);

dpNew = dpLast;
pStep = zeros(size(pLast));
pStep(jk) = dpNew;

beta = [];
alpha = [];

q_not_jk = 1:length(pLast);
q_not_jk = q_not_jk~=jk & ar.qFit;
while(true)
    q_hit_lb = pLast+pStep <= ar.lb+ar.ple.minstepsize;
    q_hit_ub = pLast+pStep >= ar.ub-ar.ple.minstepsize;
    
    if(q_hit_lb(jk) || q_hit_ub(jk)) % jk hit boundaries
        dpNew = dpNew / stepfaktor;
        if(abs(dpNew)<ar.ple.minstepsize(jk))
            if(dpLast>0)
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
    elseif sum(~isnan(last.dx))>5  % if >5 last values are available
        ub = ar.ub(jk);
        lb = ar.lb(jk);
        ytarget = icdf('chi2',0.95,1)*1.2;
        minx = ar.ple.minstepsize(jk);
        maxx = ar.ple.maxstepsize(jk);
        ss = ar.ple.samplesize(jk);

        dpNew = profileStepControl(last,lb,ub,ytarget,minx,maxx,ss);

        pStep(jk) = dpNew;
        return
            
    else
        feval(ar.ple.integrate_fkt, pLast+pStep);
        chi2trial = feval(ar.ple.merit_fkt);
        
        if(chi2trial-chi2last > dchi2*ar.ple.relchi2stepincrease(jk))
            dpNew = dpNew / stepfaktor;
            if(abs(dpNew)<ar.ple.minstepsize(jk))
%                 fprintf('WARNING: could not control step size (minstepsize = %e)\n', ar.ple.minstepsize(jk))
                dpNew = dpNew * stepfaktor;
                return
            end
        else
            if(chi2trial-chi2last < dchi2*ar.ple.relchi2stepincrease(jk) && ...
                    chi2trial-chi2last>0)
                dpNew = dpNew * stepfaktor;
                if(abs(dpNew) > ar.ple.maxstepsize(jk))
                    dpNew = ar.ple.maxstepsize(jk) * sign(dpNew);
                end
            end
            return
        end
    end

    pStep(jk) = dpNew;
end

