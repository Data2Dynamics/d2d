% Adaptive step reconciling chi^2 increase. Always keeps the increase of the
% initial chi^2 value under the boundary set by relchi2stepincrease if this
% is in accordance with minstepsize.
%
% This method uses the direction of the last step to propose a better step.

function [pStep, dpNew] = pleInitStepInertia(jk, pLast, dpLast,dpStep)

global ar

%Magic factors:
stepfaktor = ar.ple.stepfaktor(jk);
dchi2 = ar.ple.dchi2_point;

%Initialize:
chi2last = feval(ar.ple.merit_fkt);
dpNew = dpLast;
parV = 1:length(pLast);
q_not_jk = parV(:)~=jk & ar.qFit(:);

%Check whether algorithm is past his first step, where dpStep is unknown:
if sum(isnan(dpStep)) > 0
    pStep = zeros(size(pLast));
    pStep(jk) = dpNew;
else
%Rescale the last step dpStep to guess an initial step:
    c = abs(dpNew/dpStep(jk));
    pStep = c*dpStep;
end

while(true)
    c = 1; %Reset change factor
    q_hit_lb = pLast+pStep <= ar.lb+ar.ple.minstepsize;
    q_hit_ub = pLast+pStep >= ar.ub-ar.ple.minstepsize;
    %Auxiliary variables deciding which parameters hit the boundaries and 
    %which corresponding parameter steps get set to zero
    
    %Parameter jk hits boundaries:
    if(q_hit_lb(jk) || q_hit_ub(jk)) 
        c = c / stepfaktor;
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
    else
        %Other parameters hit boundaries
        if(sum(q_hit_lb(q_not_jk)) || sum(q_hit_ub(q_not_jk)))
        %Check how 'breakon' is specified
            if(sum(q_hit_lb(q_not_jk) & ar.ple.breakonlb(q_not_jk))>0 ||...
                sum(q_hit_ub(q_not_jk) & ar.ple.breakonub(q_not_jk))>0) 
                fprintf('STOP: other parameters hit hard boundaries:\n')
                fprintf('\t%s\n', ar.ple.p_labels{(q_hit_lb | q_hit_ub) & q_not_jk'})
                pStep = nan(size(pLast));
                return
            else
                %If other parameters hit hard boundaries, set their stepsize to zero
                if(sum(q_hit_lb) > 0)
                pStep(parV(q_hit_lb)) =  0;
                end
                if(sum(q_hit_ub) > 0)
                pStep(parV(q_hit_ub)) =  0;
                end
            end
        end 
        
        %Evaluate test step and obtain chi^2_trial:
        feval(ar.ple.integrate_fkt, pLast+pStep);
        chi2trial = feval(ar.ple.merit_fkt);
        
        %Decision whether stepsize needs to be decreased
        if(chi2trial-chi2last > dchi2*ar.ple.relchi2stepincrease(jk))
            c = c / stepfaktor;
            dpNew = dpNew / stepfaktor;
            if(abs(dpNew)<ar.ple.minstepsize(jk))
                fprintf('WARNING: could not control step size (minstepsize = %e)\n', ar.ple.minstepsize(jk))
                dpNew = dpNew * stepfaktor;
                % Switch to old algorithm if there is a minstepsize exception.
                % Increases stability.
                pStep = zeros(size(pLast));
                pStep(jk) = dpNew;
                return
            end
        else
            %If stepsize is small enough, the algorithm terminates
            if(chi2trial-chi2last <= dchi2*ar.ple.relchi2stepincrease(jk))
                if chi2trial-chi2last > -0.5*dchi2*ar.ple.relchi2stepincrease(jk)    
                %This condition stops increase of downward stepsizes 
                %if minimal downward change of chi^2 is already too large
                    dpNew = dpNew * stepfaktor; 
                end
                if(abs(dpNew) > ar.ple.maxstepsize(jk))
                    dpNew = ar.ple.maxstepsize(jk) * sign(dpNew);
                end
            end
            return
        end
    end
    
    %If stepsize is too large, it will be decreased and a new iteration
    %starts
    if isnan(dpStep)
        pStep(jk) = dpNew;
    else
        pStep = c*pStep;
    end
end

end
