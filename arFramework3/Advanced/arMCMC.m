% mcmc sampler
%
% function arMCMC(nruns, nburnin, method, append, nthinning)
%
%   nruns
%   method for proposal density:
%       1 = N(0,1)
%       2 = N(0,c) scaled
%       3 = MMALA
%       4 = Adaptive MCMC
%
% optimal acceptance rate ~23% (for multi-variate normal distributions), see in:
% Roberts, G.O.; Gelman, A.; Gilks, W.R. (1997).
% "Weak convergence and optimal scaling of random walk Metropolis algorithms".
% Ann. Appl. Probab. 7 (1): 110-120.
%       3 = MMALA
% Girolami, M. and Calderhead, B. (2011).
% "Riemann manifold Langevin and Hamiltonian Monte Carlo methods"
% JRSS B 73 (2): 123-214.
%       4 = Adaptive MCMC

function arMCMC(nruns, nburnin, method, append, nthinning)

global ar;

rng('shuffle');

if(~exist('nruns','var'))
    nruns = 1000;
end
nwindow = sum(ar.qFit == 1)*50;
if(~exist('method','var'))
    method = 1;
end
if(~exist('nburnin','var') || nburnin == 0)
    nburnin = 0;
    if(method==4)
        nburnin = nwindow * 50;
    end
end
if(~exist('append','var'))
    append = false;
end
if(~exist('nthinning','var'))
    nthinning = 1;
end

pReset = ar.p;
qFit = ar.qFit==1;

% initial values for chains
if(isfield(ar,'ps') && append)
    p_curr = ar.ps(end,qFit);
    ar.p = ar.ps(end,:);
else
    p_curr = ar.p(qFit);
end

if(isfield(ar,'ps') && append)
    jindexoffset = length(ar.chi2s);
    ar.ps = [ar.ps; nan(nruns, length(ar.p))];
    ar.ps_trial = [ar.ps_trial; nan(nruns, length(ar.p))];
    ar.chi2s = [ar.chi2s nan(1,nruns)];
    ar.chi2s_trial = [ar.chi2s_trial nan(1,nruns)];
    ar.acceptance = [ar.acceptance nan(1,nruns)];
else
    jindexoffset = 0;
    ar.ps = nan(nruns, length(ar.p));
    ar.ps_trial = nan(nruns, length(ar.p));
    ar.chi2s = nan(1,nruns);
    ar.chi2s_trial = nan(1,nruns);
    ar.acceptance = nan(1,nruns);
end

min_accept = 0.4;
max_accept = 0.7;

% bounds
lb = ar.lb(qFit);
ub = ar.ub(qFit);

% methods
if(method==1)
    fkt = @mcmc_norm_mvnrnd;
    adjust_scaling = false;
    use_sensis = false;
    burnin = false;
elseif(method==2)
    fkt = @mcmc_scaled_mvnrnd;
    adjust_scaling = true;
    use_sensis = false;
    burnin = true;
elseif(method==3)
    fkt = @mcmc_mmala;
    adjust_scaling = false;
    use_sensis = true;
    burnin = false;
elseif(method==4)
    fkt = @mcmc_adaptive;
    adjust_scaling = true;
    use_sensis = false;
    burnin = true;
    ps_hist = nan(nwindow, sum(ar.qFit == 1));
    ps_hist_index = 1;
    CAdapt = eye(sum(ar.qFit == 1));
end

if(~burnin)
    nburnin = 0;
end

svd_threshold = 0; 1e-8; % 0 ist besser!

Cmax = 1e8;
Cmin = 1e-8;
Cmod = 1.1;
Cfactor = (2.38/sqrt(sum(qFit)))^2 / sum(qFit);

if(~use_sensis)
    arCalcMerit(false);
    res_curr = [];
    Sres_curr = [];
else
    arCalcMerit(true);
    res_curr = ar.res;
    Sres_curr = ar.sres(:,qFit);
    CProposalPrior = diag((ar.ub(qFit)-ar.lb(qFit))/2);
end
L_curr = arGetMerit('chi2fit');
jrungo = -nburnin+1;
    
% additional functionality
do_chain_resets = false;
do_reflect_bounds = true;

% mcmc
arWaitbar(0);
naccepts = 30;
accepts = nan(1,naccepts);
i_accepts = 1;
fprintf('MCMC sampling...')
tic;
count_chain_reset = 0;
jthin = 1;
jcount = 1;
for jruns = 1:((nruns*nthinning)+nburnin)
    % figure(1); plot(accepts,'*-'); drawnow;
    qnonnanacc = ~isnan(accepts);
    accept_rate = sum(accepts(qnonnanacc))/length(accepts(qnonnanacc));
    if(jrungo>0)
        arWaitbar(jruns, (nruns*nthinning)+nburnin, sprintf('MCMC run (acceptance rate %4.1f%%)', ...
            accept_rate*100));
    else
        arWaitbar(jruns, (nruns*nthinning)+nburnin, sprintf('MCMC burn-in (acceptance rate %4.1f%%)', ...
            accept_rate*100));
    end
    
    [mu_curr, covar_curr] = feval(fkt, p_curr, res_curr, Sres_curr);
    p_trial = mvnrnd(mu_curr, covar_curr);
    L_trial = 0;
    
    % reflect from bejond bounds
    if(do_reflect_bounds == true)
        if(sum(p_trial<lb) + sum(p_trial>ub) > 0)
            qlb = p_trial<lb;
            p_trial(qlb) = p_trial(qlb) + 2*(lb(qlb) - p_trial(qlb));
            qub = p_trial>ub;
            p_trial(qub) = p_trial(qub) + 2*(ub(qub) - p_trial(qub));
        end
    end
    
    if(sum(p_trial<lb) + sum(p_trial>ub) == 0) % check bounds
        % fprintf('#%i bounds ok\n', jrungo);
        try
            if(~use_sensis)
                arCalcMerit(false,p_trial);
                res_trial = [];
                Sres_trial = [];
            else
                arCalcMerit(true,p_trial);
                res_trial = ar.res;
                Sres_trial = ar.sres(:,qFit);
            end
            L_trial = arGetMerit('chi2fit');
            [mu_trial, covar_trial] = feval(fkt, p_trial, res_trial, Sres_trial);
            Q_trial = mvnpdf(p_trial, mu_curr, covar_curr);
            Q_curr = mvnpdf(p_curr, mu_trial, covar_trial);
            a = exp(-0.5*(L_trial - L_curr)) * (Q_curr / Q_trial);
            randa = rand;
            qa = randa <= min([1 a]); % accept?
            % fprintf('#%i %i %e\t%e\t%e\t%e\t%e\n', jrungo, qa, a, L_trial, L_curr, Q_curr, Q_trial)
        catch err_id
            fprintf('#%i error: %s\n', jrungo, err_id.message);
            qa = false;
        end
    else
        % fprintf('#%i bound violation\n', jrungo);
        qa = false;
    end
    
    % adjust scaling during burn-in
    if(jrungo<=0 && adjust_scaling)
        if(accept_rate > max_accept && Cfactor*Cmod<Cmax)
            Cfactor = Cfactor*Cmod;
        elseif(accept_rate < min_accept && Cfactor/Cmod>Cmin)
            Cfactor = Cfactor/Cmod;
        end
        if(method==2)
            fprintf('#%i Cfactor = %g, Acceptance Rate = %4.1f%%\n', jruns, Cfactor, accept_rate*100);
        end
    end
    
    accepts(i_accepts) = qa;
    i_accepts = i_accepts + 1;
    if(i_accepts>length(accepts))
        i_accepts = 1;
    end
    
    % reset chain ?
    if(do_chain_resets == true)
        if(jtrials==100 && jruns+jindexoffset-1>0)
            p_curr = ar.ps(randi(jruns+jindexoffset-1,1),:);
            if(~use_sensis)
                arCalcMerit(false,p_curr);
                res_curr = [];
                Sres_curr = [];
            else
                arCalcMerit(true,p_curr);
                res_curr = ar.res;
                Sres_curr = ar.sres(:,qFit);
            end
            L_curr = arGetMerit('chi2fit');
            jtrials = 1;
            count_chain_reset = count_chain_reset + 1;
        end
    end
    
    % update
    if(qa)
        p_curr = p_trial;
        L_curr = L_trial;
        
        if(method==4)
            ps_hist(ps_hist_index,:) = p_curr;
            ps_hist_index = ps_hist_index + 1;
            if(ps_hist_index > nwindow)
                ps_hist_index = 1;
            end
        end
    end
    
    % save samples
    if(jrungo>0)
        jthin = jthin + 1;
        if(jthin > nthinning)
            ar.ps(jcount+jindexoffset,:) = pReset;
            ar.ps(jcount+jindexoffset,qFit) = p_curr;
            ar.ps_trial(jcount+jindexoffset,:) = pReset;
            ar.ps_trial(jcount+jindexoffset,qFit) = p_trial;
            ar.chi2s(jcount+jindexoffset) = L_curr;
            ar.chi2s_trial(jcount+jindexoffset) = L_trial;
            ar.acceptance(jcount+jindexoffset) = accept_rate;
            jthin = 1;
            jcount = jcount + 1;
        end
    end
    
    jrungo = jrungo + 1;
end
arWaitbar(-1);
fprintf('done (%s, %i chain resets) \n', secToHMS(toc), count_chain_reset);
if(isfield(ar,'mcmc_toc') && append)
    ar.mcmc_toc = toc + ar.mcmc_toc;
else
    ar.mcmc_toc = toc;
end
ar.p = pReset;


% N(0,1)
    function [mu, covar] = mcmc_norm_mvnrnd(ptmp, ~, ~)
        mu = ptmp;
        covar = eye(length(ptmp)) * 1e-8;
    end

% N(0,c) scaled
    function [mu, covar] = mcmc_scaled_mvnrnd(ptmp, ~, ~)
        mu = ptmp;
        covar = eye(length(ptmp)) * Cfactor;
    end

%     % MMALA (simplified)
%     function [mu, covar] = mcmc_mmala(ptmp, restmp, srestmp)
%         beta = - restmp*srestmp;
%         alpha = srestmp'*srestmp;
%         alpha_dash = alpha;
%         alpha_dash(eye(size(alpha_dash))==1) = ...
%             alpha(eye(size(alpha_dash))==1) * (1 + Cfactor); % (15.5.13)
%
%         % solve with SVD regulatization
%         [U,S,V] = svd(alpha_dash);
%         s = diag(S);
%         qs = (s/max(s)) > svd_threshold;
%
%         deltap = zeros(size(beta));
%         for jj=find(qs)'
%             deltap = deltap + transpose((U(:,jj)'*beta'/S(jj,jj))*V(:,jj));
%         end
%         mu = ptmp + deltap;
%
%         covar = inv(alpha_dash + inv(CProposalPrior));
%     end

% MMALA (simplified)
    function [mu, covar] = mcmc_mmala(ptmp, restmp, srestmp)
        beta = - restmp*srestmp;
        alpha = srestmp'*srestmp;
        alpha_dash = alpha + inv(CProposalPrior);
        
        % solve with SVD regulatization
        [U,S,V] = svd(alpha_dash);
        s = diag(S);
        qs = (s/max(s)) > svd_threshold;
        
        deltap = zeros(size(beta));
        for jj=find(qs)'
            deltap = deltap + transpose((U(:,jj)'*beta'/S(jj,jj))*V(:,jj));
        end
        mu = ptmp + deltap;
        
        covar = inv(alpha_dash);
    end

% Adaptive MCMC
    function [mu, covar] = mcmc_adaptive(ptmp, ~, ~)
        if(jrungo<=0)
            mu = ptmp;
            covar = eye(length(ptmp)) * Cfactor;
        else
            qnonnan = ~isnan(ps_hist(:,1));
            CAdapt = cov(ps_hist(qnonnan,:));
            mu = ptmp;
            covar = CAdapt;
        end
    end
end



