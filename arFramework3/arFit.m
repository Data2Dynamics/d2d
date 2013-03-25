% Fit model parameters to data using lsqnonlin
%
% arFit(silent)
%
% lsqnonlin.m exit flag description:
%       1  LSQNONLIN converged to a solution X.
%       2  Change in X smaller than the specified tolerance.
%       3  Change in the residual smaller than the specified tolerance.
%       4  Magnitude search direction smaller than the specified tolerance.
%       0  Maximum number of function evaluations or of iterations reached.
%      -1  Algorithm terminated by the output function.
%      -2  Bounds are inconsistent.
%      -4  Line search cannot sufficiently decrease the residual along the
%           current search direction.
%

function arFit(silent)

global ar

if(~isfield(ar.config, 'optimizer'))
    ar.config.optimizer = 1;
end
if(~isfield(ar.config, 'showFitting'))
    ar.config.showFitting = 0;
end

if(nargin==0)
    silent = false;
end

ar.fevals = 0;

if(ar.config.useSensis)
    ar.config.optim.Jacobian = 'on';
else
    ar.config.optim.Jacobian = 'off';
end

ar.fit.chi2_hist = nan(1,ar.config.optim.MaxIter);
ar.fit.p_hist = nan(ar.config.optim.MaxIter,length(ar.p));
ar.fit.maxstepsize_hist = nan(1,ar.config.optim.MaxIter);
ar.fit.stepsize_hist = nan(1,ar.config.optim.MaxIter);
ar.config.optim.OutputFcn = @arPlotFast;

arChi2(false);
chi2_old = ar.chi2fit;

ub = ar.ub;
lb = ar.lb;
ub(ar.type==2) = ub(ar.type==2) + 1;
lb(ar.type==2) = lb(ar.type==2) - 1;
ub = ub(ar.qFit==1);
lb = lb(ar.qFit==1);

% lsqnonlin
if(ar.config.optimizer == 1)     
    [pFit, chi2, resnorm, exitflag, output, lambda, jac] = ...
        lsqnonlin(@merit_fkt, ar.p(ar.qFit==1), lb, ub, ar.config.optim);
    
% fmincon
elseif(ar.config.optimizer == 2) 
    [pFit, chi2, exitflag, output, lambda, jac] = ...
        fmincon(@merit_fkt_fmincon, ar.p(ar.qFit==1),[],[],[],[],lb,ub, ...
        [],optimset(ar.config.optim,'DerivativeCheck','off', ...
        'GradObj','on','Hessian','on','Jacobian', ...
        [],'Display','off'));
    resnorm = merit_fkt(pFit);
    
% levenberg-marquardt
elseif(ar.config.optimizer == 3) 
    arFitLM(silent);
    return;
   
% STRSCNE
elseif(ar.config.optimizer == 4)
    warnreset = warning;
    warning('off','MATLAB:rankDeficientMatrix');
    [pFit, exitflag, output, history] = ...
        STRSCNE(ar.p(ar.qFit==1), @merit_fkt2, [-Inf,0], lb, ub, [1000,1000,1,1], @merit_dfkt2);
    warning(warnreset);
    ar.p(ar.qFit==1) = pFit;
    ar.fit.exitflag = exitflag;
    ar.fit.output = output;
    ar.fit.history = history;
    
    arChi2(false);
    fprintf('STRSCNE finished after %i iterations: code %i, total chi2 improvement = %g\n', ...
        output(1), exitflag, chi2_old - ar.chi2fit);
    
    return;

% arNLS
elseif(ar.config.optimizer == 5)     
    [pFit, chi2, resnorm, exitflag, output, lambda, jac] = ...
        arNLS(@merit_fkt, ar.p(ar.qFit==1), lb, ub, ar.config.optim);
end

if(isfield(ar, 'ms_count_snips') && ar.ms_count_snips>0)
    if(max(ar.ms_violation) > ar.ms_threshold)
        fprintf('Multiple Shooting: continuity constains violated %e > %e\n', max(ar.ms_violation), ar.ms_threshold);
    end
end

ar.p(ar.qFit==1) = pFit;
ar.fit.exitflag = exitflag;
ar.fit.output = output;
ar.fit.iter = output.iterations;
ar.fit.chi2 = chi2;
ar.fit.lambda = lambda;
ar.fit.qFit = ar.qFit;
ar.fit.res = resnorm;
ar.fit.sres = full(jac);

if(~silent || exitflag < 1)
    outputstr = '';
    switch exitflag
        case 1
            outputstr = 'Converged to a solution';
        case 2  
            outputstr = 'Change in X too small';
        case 3  
            outputstr = 'Change in RESNORM too small';
        case 4  
            outputstr = 'Computed search direction too small';
        case 0  
            outputstr = 'Too many function evaluations or iterations';
        case -1  
            outputstr = 'Stopped by output/plot function';
        case -2  
            outputstr = 'Bounds are inconsistent';
        case -3  
            outputstr = 'Regularization parameter too large (Levenberg-Marquardt)';
        case -4  
            outputstr = 'Line search failed';
        case -98
            outputstr = sprintf('Multiple Shooting: constraint strength not controlable > %e\n', 1e20);
        case -99
            outputstr = sprintf('Multiple Shooting: mean constraint violation > %e\n', ar.ms_treshold);
    end
    arChi2(false);
    
    switch ar.config.optimizer
        case 1
            fprintf('lsqnonlin ');
        case 2
            fprintf('fmincon ');
        case 6
            fprintf('arNLS ');
    end
    
    fprintf('finished after %i iterations: %s, total chi2 improvement = %g\n', ...
        ar.fit.output.iterations, outputstr, chi2_old - ar.chi2fit);
end

if(~silent)
    arChi2;
end

function [chi2, schi2,H] = merit_fkt_fmincon(pTrial)
[res, sres] = merit_fkt(pTrial);
chi2 = sum(res.^2);
sresSq = 2*(res'*ones(1,length(pTrial))) .*sres;
schi2 = sum(sresSq);

H = sresSq'*sresSq;


function [res, sres] = merit_fkt(pTrial)

global ar

arChi2(ar.config.useSensis, pTrial)

res = ar.res;
if(nargout>1 && ar.config.useSensis)
    sres = ar.sres(:, ar.qFit==1);
end


function stop = arPlotFast(~,optimValues,state)

global ar

if(strcmp(state, 'iter'))
    ar.fit.chi2_hist(optimValues.iteration+1) = ar.chi2fit;
    ar.fit.p_hist(optimValues.iteration+1,:) = ar.p;
    
    if(ar.config.optimizer == 6)
        ar.fit.maxstepsize_hist(optimValues.iteration+1) = optimValues.mu;
        ar.fit.stepsize_hist(optimValues.iteration+1) = optimValues.normdp;
    end
    
    if(ar.config.showFitting)
        arPlot(false, true, false, true, true);
        drawnow;
    end
end

stop = false;


% for STRSCNE
function res = merit_fkt2(pTrial)
global ar
arChi2(ar.config.useSensis, pTrial)
res = ar.res';

function sres = merit_dfkt2(~)
global ar
if(ar.config.useSensis)
    sres = ar.sres(:, ar.qFit==1);
end
