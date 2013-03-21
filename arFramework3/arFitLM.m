% Fit model parameters to data using Levenberg-Marquardt
%
% K. Levenberg: A Method for the Solution of Certain Problems in Least Squares, Quart. Appl. Math. 2, 164-168, 1944.
% D. Marquardt: An Algorithm for Least-Squares Estimation of Nonlinear Parameters, SIAM J. Appl. Math. 11, 431-441, 1963.
%
% arFitLM(silent, verbose)
%
%  VERBOSE LEVELS
%  0 = no output
%  1 = compact
%  2 = normal
%  3 = detailed with plots
%  4 = very detailed with plots
%
% Copyright (C) 2008 Andreas Raue <andreas.raue@me.com>

function arFitLM(silent, verbose)

global ar

if(nargin==0)
    silent = false;
end

if(nargin<2)
    verbose = 1;
end

if(silent)
    verbose = 0;
end

qFit = ar.qFit==1;
if(sum(qFit)==0)
    error('arFitLM: there need to be at least one free parameter');
end

sumviopars = sum(ar.p<ar.lb | ar.p>ar.ub);
if(sumviopars>0)
    error('arFitLM: %i initial guess violate bounds', sumviopars);
end

ar.fit = [];

maxiter = ar.config.optim.MaxIter;
if(isempty(maxiter))
    maxiter = 1000;
end

ar.fit.p_hist = nan(maxiter+1, length(ar.p));
ar.fit.chi2_hist = nan(1, maxiter+1);
ar.fit.lambda_hist = nan(1, maxiter+1);
ar.fit.svd_hist = nan(1, maxiter+1);

% working lambda
lambda = 0.001; % 1=conservative, 0.001=progressive

% limits for reaching convergence
lambda_max = 1e8;

% SVD regulatization threshold (NR: chapter 15.4)
svd_threshold = 0; % 1e-6; % 

% init
pBest = ar.p;
chi2_last = ar.chi2fit;

ar.fit.exitflag = 4;
ar.fit.message = 'ERROR';
exitflag_labels = {sprintf('converged (lambda > %g', lambda_max), ...
   '', ...
   sprintf('NOT-CONVERGED (after %i iterations', maxiter), ...
   'ERROR (', ...
   'ERROR (alpha_dash has NaNs', ...
   'ERROR (solving: alpha_dash * deltap = beta has NaNs'};

%% begin of computation
try
    % to make it faster, refresh only after successful step
    refresh = 1;
    
    iter = 1;
    beta = zeros(sum(qFit), 1);
    alpha = zeros(sum(qFit));
    
    % save initial to history
    ar.fit.p_hist(iter,:) = ar.p;
    ar.fit.chi2_hist(iter) = ar.chi2fit;
    ar.fit.lambda_hist(iter) = lambda;
    
    itertrials = 0;
    
    while iter<=maxiter
        itertrials = itertrials+1;
        
        if(refresh==1)
            % backup old status
            p_old = ar.p;
            chi2_old = ar.chi2fit;
            
            % calculate alpha & beta respective dChi^2/dp & d^2Chi^2/dp_i/pd_j
            arChi2(true);
            beta = - ar.res * ar.sres(:,qFit);
            alpha = ar.sres(:,qFit)' * ar.sres(:,qFit);
        end
        
        % calculate alpha'(k,l), see NR (15.5.13), Levenberg-Marquardt step
        beta_dash = beta;
        alpha_dash = alpha;
%         alpha_dash(eye(size(alpha_dash))==1) = alpha(eye(size(alpha_dash))==1) * (1+lambda);
        
        % original step proposed by Levenberg earlier.
        alpha_dash(eye(size(alpha_dash))==1) = alpha(eye(size(alpha_dash))==1) + lambda;
        
        % check for NaNs in alpha_dash
        if(sum(isnan(alpha_dash(:)))>0)
            ar.fit.exitflag = 5;
            break
        end
                
        % SVD regulatization
        [U,S,V] = svd(alpha_dash);
        s = diag(S);
        qs = (s/max(s)) > svd_threshold;
        ar.fit.svd_hist(iter) = sum(qs);
        
        % SVD solve
        deltap = zeros(size(beta_dash));
        for jj=find(qs)'
            deltap = deltap + transpose((U(:,jj)'*beta_dash'/S(jj,jj))*V(:,jj));
        end
        
        % check parameter on bounds and reduce search space accordingly
        qponbounds = (ar.p(qFit)==ar.lb(qFit) & deltap<0) | ...
            (ar.p(qFit)==ar.ub(qFit) & deltap>0);
        if(sum(~qponbounds)==0) % exit is all on bounds
            ar.fit.exitflag = 1;
            break
        end
        if(sum(qponbounds)>0)
            alpha_dash = alpha_dash(~qponbounds,~qponbounds);
            beta_dash = beta_dash(~qponbounds);
            
            % SVD regulatization
            [U,S,V] = svd(alpha_dash);
            s = diag(S);
            qs = (s/max(s)) > svd_threshold;
            ar.fit.svd_hist(iter) = sum(qs);
            
            % SVD solve
            deltap(:) = 0;
            for jj=find(qs)'
                deltap(~qponbounds) = deltap(~qponbounds) + transpose((U(:,jj)'*beta_dash'/S(jj,jj))*V(:,jj));
            end
        end
        
        % check for NaNs in proposed step direction
        if(sum(isnan(deltap))>0)
            ar.fit.exitflag = 6;
            break
        end
        
        % check bounds and shorten step if beyond
        cupper = deltap ./ (ar.ub(qFit) - p_old(qFit));
        clower = deltap ./ (ar.lb(qFit) - p_old(qFit));
        cfactor = max([cupper(~qponbounds) clower(~qponbounds)]);
        if(cfactor > 1)
            deltap = deltap / cfactor;
        end
        
        % set new candidate parameter values
        ar.p(qFit) = p_old(qFit) + deltap;
        ar.p(ar.p<ar.lb) = ar.lb(ar.p<ar.lb);
        ar.p(ar.p>ar.ub) = ar.ub(ar.p>ar.ub);
        
        sumviopars = sum(ar.p<ar.lb | ar.p>ar.ub);
        if(sumviopars>0)
            keyboard
            error('arFitLM: %i parameter violate bounds', sumviopars);
        end
        
        % evaluate objective function
        arChi2(false);
        
        if(verbose>2)
            fprintf(' %4i : lambda=%5.2e stepsize=%5.2e dchisqr=%+5.2e svd=%2i ', ...
                iter, lambda, sqrt(sum(deltap.^2)), ar.chi2fit-chi2_old, sum(qs));
            if (verbose>3)
                fprintf('%7.3e ', ar.p(qFit));
            end
        end
        
        % compare old and new chi^2
        if(ar.chi2fit>=chi2_old) % if increased, reset & do it again with bigger lambda
            if(verbose>1)
                fprintf('-');
            end
            
            lambda = lambda * 10;
            refresh = 0;
            
            % backup
            ar.p = p_old;
            ar.chi2fit = chi2_old;
            
            % stop, if lambda limit reached
            if(lambda>lambda_max)
                ar.fit.exitflag = 1;
                break
            end
        else % if decreased, keep the step
            iter = iter + 1;
            
            if(verbose>1)
                fprintf('*');
            end
            
            pBest = ar.p;
            
            % save history
            ar.fit.p_hist(iter,:) = ar.p;
            ar.fit.chi2_hist(iter) = ar.chi2fit;
            ar.fit.lambda_hist(iter) = lambda;
            
            if(ar.fit.chi2_hist(iter)>ar.fit.chi2_hist(iter-1))
                error('lala');
            end
            
            % decrease lambda
            lambda = lambda / 10;
            
            refresh = 1;
        end
        
        if(verbose>2 || (mod(itertrials,80)==0 && verbose>1))
            fprintf('\n');
        end
    end
catch err_id
    ar.fit.message = ['ERROR (' err_id.message];
    rethrow(err_id)
end

if(verbose>1)
    fprintf('\n');
end

ar.fit.maxiter = maxiter;
ar.fit.iter = iter;
ar.fit.itertrials = itertrials;

ar.p = pBest;
arChi2(false);

ar.fit.dchi2 = chi2_last - ar.chi2fit;
ar.fit.res = ar.res;
ar.fit.sres = ar.sres(:,qFit);

if(iter==maxiter+1)
    ar.fit.exitflag=3;
end

ar.fit.message = exitflag_labels{ar.fit.exitflag};
if(verbose>0 || ar.fit.exitflag>2)
    disp([ar.fit.message sprintf(', %i iterations, total chi2 improvement = %g)', ...
        ar.fit.iter, ar.fit.dchi2)]);
end

if(~silent)
    arChi2;
end

