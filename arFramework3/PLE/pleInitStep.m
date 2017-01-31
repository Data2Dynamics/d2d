% adaptive step reconciling chi^2 increase
% progressive step method

function [pStep, dpNew, beta, alpha] = pleInitStep(jk, pLast, dpLast, lastAll)

global ar

% mag factor
stepfaktor = 2;

if(~ar.ple.breakon_point)
    dchi2 = ar.ple.dchi2;
else
    dchi2 = ar.ple.dchi2_point;
end

dpNew = dpLast;

q_fit = ar.qFit;
jk_fit = find(find(q_fit)==jk);
q_step_fit = 1:sum(q_fit) ~= jk_fit;
q_step = q_fit;
q_step(jk) = false;

chi2last = feval(ar.ple.merit_fkt);

opts = optimset('Display','off');

[beta, alpha] = feval(ar.ple.diffmerit_fkt);
opterm = 2*alpha(q_step_fit,q_step_fit);
rhsterm = beta(q_step_fit) + 2*alpha(q_step_fit,jk_fit)*dpNew;
solterm = svdsolve(opterm, rhsterm);

% %%
% ar.ple.svd_threshold = 1e-6;
% 
% opterm = 2*alpha(q_step_fit,q_step_fit);
% rhsterm = beta(q_step_fit) + 2*alpha(q_step_fit,jk_fit)*dpNew;
% 
% solterm1 = svdsolve(opterm, rhsterm);
% solterm2 = pinv(opterm) * rhsterm;
% solterm3 = opterm\rhsterm;
% solterm4 = quadprog(opterm, -rhsterm, [], [], [], [], ar.lb(q_step_fit), ar.ub(q_step_fit), [], opts);
% solterm5 = quadprog(opterm, -rhsterm, [], [], [], [], ar.lb(q_step_fit)-pLast(q_step_fit), ar.ub(q_step_fit)-pLast(q_step_fit), [], opts);
% % [solterm1 solterm2 solterm3 solterm4 solterm5]
% [solterm1 solterm4 solterm5]
% 
% %%

% solterm2 = quadprog(opterm, rhsterm, [], [], [], [], ar.lb(q_step_fit)-pLast(q_step_fit), ar.ub(q_step_fit)-pLast(q_step_fit), [], opts);
pStep = zeros(size(dpNew));
pStep(jk) = dpNew;
pStep(q_step) = -transpose(solterm);

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
    else
        feval(ar.ple.integrate_fkt, pLast+pStep);
        chi2trial = feval(ar.ple.merit_fkt);
        
        if(chi2trial-chi2last > dchi2*ar.ple.relchi2stepincrease(jk))
            dpNew = dpNew / stepfaktor;
            if(abs(dpNew)<ar.ple.minstepsize(jk))
                fprintf('WARNING: could not control step size (minstepsize = %e)\n', ar.ple.minstepsize(jk))
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

    rhsterm = beta(q_step_fit) + 2*alpha(q_step_fit,jk_fit)*dpNew;
    solterm = svdsolve(opterm, rhsterm);
%     solterm = quadprog(opterm, rhsterm, [], [], [], [], ar.lb(q_step_fit), ar.ub(q_step_fit), [], opts);
    pStep = zeros(size(dpNew));
    pStep(jk) = dpNew;
    pStep(q_step) = -transpose(solterm);
end


% linear solver using SVD regularization
function solterm = svdsolve(opterm, rhsterm)
;

[U,S,V] = svd(opterm);

% repress big steps with small impact on Chi^2
s = transpose(diag(S)/max(diag(S)));
qs = s<ar.ple.svd_threshold;

% SVD solve
solterm = zeros(size(rhsterm));
for jj=1:size(U,2)
    if(~qs(jj))
        solterm = solterm + (U(:,jj)'*rhsterm/S(jj,jj))*V(:,jj);
    end
end
