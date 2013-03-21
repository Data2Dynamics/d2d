% Fit model parameters to data
%
% arFit2(silent, verbose)
%
%  VERBOSE LEVELS
%  0 = no output
%  1 = compact
%  2 = normal
%
% Copyright (C) 2013 Andreas Raue <andreas.raue@fdm.uni-freiburg.de>

function arFit2(silent, verbose)

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
    error('arFit: there need to be at least one free parameter');
end

sumviopars = sum(ar.p<ar.lb | ar.p>ar.ub);
if(sumviopars>0)
    error('arFit: %i initial parameter values violate bounds', sumviopars);
end

p = ar.p(qFit);
lb = ar.lb(qFit);
ub = ar.ub(qFit);
np = length(p);

ar.fit = [];

maxiter = ar.config.optim.MaxIter;
if(isempty(maxiter))
    maxiter = 1000;
end

tolfun = ar.config.optim.TolFun;
if(isempty(tolfun))
    tolfun = 1e-6;
end
tolx = ar.config.optim.TolX;
if(isempty(tolx))
    tolx = 1e-6;
end

ar.fit.p_hist = nan(maxiter+1, length(ar.p));
ar.fit.chi2_hist = nan(1, maxiter+1);
ar.fit.stepsize_hist = nan(1, maxiter+1);
ar.fit.maxstepsize_hist = nan(1, maxiter+1);
ar.fit.exitflag = -1;

% init
arChi2(true);
pBest = p;
chi2Initial = ar.chi2fit;

%% begin of computation

% step increase factor
stepsizefactor = 2;
stepsizewindowmax = 5;

% to make it faster, refresh only after successful step
refresh = 1;
iter = 1;
maxstepsize = 1;
stepsizewindow = 0;

% save initial to history
ar.fit.p_hist(iter,:) = ar.p;
ar.fit.chi2_hist(iter) = ar.chi2fit;
ar.fit.stepsize_hist(iter) = 0;
ar.fit.maxstepsize_hist(iter) = maxstepsize;

quadprog_optims = optimset('Display', 'off');
warnreset = warning;
warning('off','optim:lsqlin:LinConstraints');
warning('off','optim:quadprog:SwitchToMedScale');
warning('off','MATLAB:nearlySingularMatrix');

ar.chi2tol = 0;
exitflag = 0;
deltap = zeros(length(pBest),1);
try
    while iter<=maxiter
        if(refresh==1)
            % backup old status
            p_old = ar.p(qFit) + 0;
            chi2_old = ar.chi2fit + 0;
            
            res = ar.res;
            sres = ar.sres(:,qFit);
        end
        
        % sub problem search space
        lbtmp = lb-p_old;
        ubtmp = ub-p_old;
        lbtmp(lbtmp < -maxstepsize) = -maxstepsize;
        ubtmp(ubtmp > maxstepsize) = maxstepsize;
        
        % solve sub problem
        deltap = lsqlin(-sres, res, [], [], [], [], ...
            lbtmp, ubtmp, [], quadprog_optims);
        
        % set new candidate parameter values
        ar.p(qFit) = p_old + deltap';
        
        % evaluate objective function
        arChi2(true);
        
        if(verbose>2)
            fprintf(' %4i : chisqr=%+5.2e stepsize=%5.2e <= maxstepsize=%5.2e dchisqr=%+5.2e ', ...
                iter, ar.chi2fit, sqrt(sum(deltap.^2)), sqrt(np)*maxstepsize, ar.chi2fit-chi2_old);
            if (verbose>3)
                fprintf('%7.3e ', ar.p(qFit));
            end
        end
        
        % compare old and new chi^2
        if(ar.chi2fit < (chi2_old + ar.chi2tol)) % if decreased, keep the step
            iter = iter + 1;
            refresh = 1;
            
            % increase step size
            stepsizewindow = stepsizewindow + 1;
            %if(stepsizewindow>stepsizewindowmax && sqrt(sum(deltap.^2))/(sqrt(np)*maxstepsize)>0.9)
            if(stepsizewindow>stepsizewindowmax)
                maxstepsize = maxstepsize * stepsizefactor;
                stepsizewindow = 0;
            end
            
            if(verbose>1)
                fprintf('*');
            end
            
            pBest = ar.p(qFit);
            
            % save history
            ar.fit.p_hist(iter,:) = ar.p;
            ar.fit.chi2_hist(iter) = ar.chi2fit;
            ar.fit.maxstepsize_hist(iter) = sqrt(np)*maxstepsize;
            ar.fit.stepsize_hist(iter) = sqrt(sum(deltap.^2));
        else % if increased, reset & do it again with smaller step
            refresh = 0;
            
            % decrease step size
            maxstepsize = maxstepsize / stepsizefactor;
            stepsizewindow = 0;
            
            if(verbose>1)
                fprintf('-');
            end
            
            % backup
            ar.p(qFit) = p_old;
            
            if(sqrt(np)*maxstepsize < tolx)
                exitflag = 2;
                break;
            end
        end
        
        if(verbose>2)
            fprintf('\n');
        end
    end
    if(verbose>1)
        fprintf('\n');
    end
catch err_id
    exitflag = -1;
end

warning(warnreset);

ar.fit.maxiter = maxiter;
ar.fit.iter = iter;
ar.fit.exitflag = exitflag;

ar.p(qFit) = pBest;
arChi2(true);

% save final to history
ar.fit.p_hist(iter+1,:) = ar.p;
ar.fit.chi2_hist(iter+1) = ar.chi2fit;
ar.fit.maxstepsize_hist(iter+1) = sqrt(np)*maxstepsize;
ar.fit.stepsize_hist(iter+1) = sqrt(sum(deltap.^2));

ar.fit.dchi2 = chi2Initial - ar.chi2fit;
ar.fit.res = ar.res;
ar.fit.sres = ar.sres(:,qFit);
    
if(~silent)
    fprintf('%i iterations, total chi2 improvement = %g\n', ...
        ar.fit.iter-1, ar.fit.dchi2);
    arChi2;
end

if(exitflag == -1)
    rethrow(err_id);
end



