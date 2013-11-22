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
if(~isfield(ar.config, 'optimizerStep'))
    ar.config.optimizerStep = 0;
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

ar.fit.iter_count = 0;
ar.fit.chi2_hist = nan(1,ar.config.optim.MaxIter);
ar.fit.constr_hist = nan(1,ar.config.optim.MaxIter);
ar.fit.p_hist = nan(ar.config.optim.MaxIter,length(ar.p));
ar.fit.maxstepsize_hist = nan(1,ar.config.optim.MaxIter);
ar.fit.stepsize_hist = nan(1,ar.config.optim.MaxIter);
ar.config.optim.OutputFcn = @arPlotFast;

arChi2(true, []);
chi2_old = ar.chi2fit;

ub = ar.ub;
lb = ar.lb;
ub(ar.type==2) = ub(ar.type==2) + 1;
lb(ar.type==2) = lb(ar.type==2) - 1;
ub = ub(ar.qFit==1);
lb = lb(ar.qFit==1);

% lsqnonlin
if(ar.config.optimizer == 1)     
    [pFit, ~, resnorm, exitflag, output, lambda, jac] = ...
        lsqnonlin(@merit_fkt, ar.p(ar.qFit==1), lb, ub, ar.config.optim);
    
% fmincon
elseif(ar.config.optimizer == 2)
    options = optimset('fmincon');
    options.GradObj = 'on';
    options.GradConstr = 'on';
    options.TolFun = ar.config.optim.TolFun;
    options.TolX = ar.config.optim.TolX;
    options.Display = ar.config.optim.Display;
    options.MaxIter = ar.config.optim.MaxIter;
    options.OutputFcn = ar.config.optim.OutputFcn;
    
    options.Algorithm = 'interior-point';
    options.SubproblemAlgorithm = 'cg';
    % options.Hessian = 'fin-diff-grads';
    options.Hessian = 'user-supplied';
    options.HessFcn = @fmincon_hessianfcn;
    % options2.InitBarrierParam = 1e+6;
    % options2.InitTrustRegionRadius = 1e-1;

    [pFit, ~, exitflag, output, lambda, jac] = ...
        fmincon(@merit_fkt_fmincon, ar.p(ar.qFit==1),[],[],[],[],lb,ub, ...
        @confun,options);
    resnorm = merit_fkt(pFit);
    
% PSO
elseif(ar.config.optimizer == 3) 
    [pFit, ~, resnorm, exitflag, output, lambda, jac] = ...
        arFitPSO(lb, ub, ar.config.optim);
   
% STRSCNE
elseif(ar.config.optimizer == 4)
    warnreset = warning;
    warning('off','MATLAB:rankDeficientMatrix');
    [pFit, exitflag, output, history] = ...
        STRSCNE(ar.p(ar.qFit==1), @merit_fkt_STRSCNE, [-Inf,0], lb, ub, [1000,1000,1,1], @merit_dfkt_STRSCNE);
    warning(warnreset);
    ar.p(ar.qFit==1) = pFit;
    ar.fit.exitflag = exitflag;
    ar.fit.output = output;
    ar.fit.history = history;
    
    arChi2(false, []);
    fprintf('STRSCNE finished after %i iterations: code %i, total chi2 improvement = %g\n', ...
        output(1), exitflag, chi2_old - ar.chi2fit);
    
    return;

% arNLS
elseif(ar.config.optimizer == 5)     
    [pFit, ~, resnorm, exitflag, output, lambda, jac] = ...
        arNLS(@merit_fkt, ar.p(ar.qFit==1), lb, ub, ar.config.optim, ar.config.optimizerStep);
    
% fmincon as least squares fit
elseif(ar.config.optimizer == 6)
    options = optimset('fmincon');
    options.GradObj = 'on';
    options.TolFun = ar.config.optim.TolFun;
    options.TolX = ar.config.optim.TolX;
    options.Display = ar.config.optim.Display;
    options.MaxIter = ar.config.optim.MaxIter;
    options.OutputFcn = ar.config.optim.OutputFcn;
    
    [pFit, ~, exitflag, output, lambda, jac] = ...
        fmincon(@merit_fkt_fmincon_lsq, ar.p(ar.qFit==1),[],[],[],[],lb,ub, ...
        [],options);
    resnorm = merit_fkt(pFit);
    
else
    error('ar.config.optimizer invalid');
end

if(isfield(ar, 'ms_count_snips') && ar.ms_count_snips>0)
    if(max(ar.ms_violation) > ar.ms_threshold)
        fprintf('Multiple Shooting: continuity constains violated %e > %e\n', max(ar.ms_violation), ar.ms_threshold);
    end
end

ar.p(ar.qFit==1) = pFit;
arChi2(true, []);

ar.fit.exitflag = exitflag;
ar.fit.output = output;
ar.fit.iter = output.iterations;
ar.fit.chi2 = ar.chi2fit;
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
        case 50
            outputstr = 'MaxIter Reached';
        case 51
            outputstr = 'FitCount Reache';
    end
    
    fprintf('%s finished after %i iterations: %s, total improvement = %g\n', ...
        ar.config.optimizers{ar.config.optimizer}, ...
        ar.fit.output.iterations, outputstr, chi2_old - ar.chi2fit);
end

if(~silent)
    arChi2(true);
end


% lsqnonlin and arNLS
function [res, sres] = merit_fkt(pTrial)
global ar
arChi2(ar.config.useSensis, pTrial)
res = [ar.res ar.constr];
if(nargout>1 && ar.config.useSensis)
    sres = [];
    if(~isempty(ar.sres))
        sres = ar.sres(:, ar.qFit==1);
    end
    if(~isempty(ar.sconstr))
        sres = [sres; ar.sconstr(:, ar.qFit==1)];
    end
end

% fmincon
function [l, g, H] = merit_fkt_fmincon(pTrial)
global ar
arChi2(ar.config.useSensis, pTrial)
l = sum(ar.res.^2);
if(nargout>1)
    g = ar.res*ar.sres(:, ar.qFit==1);
end
if(nargout>2)
    H = ar.sres(:, ar.qFit==1)'*ar.sres(:, ar.qFit==1);
end

% fmincon as lsq
function [l, g] = merit_fkt_fmincon_lsq(pTrial)
global ar
arChi2(ar.config.useSensis, pTrial)
res = [ar.res ar.constr];
if(nargout>1 && ar.config.useSensis)
    sres = [];
    if(~isempty(ar.sres))
        sres = ar.sres(:, ar.qFit==1);
    end
    if(~isempty(ar.sconstr))
        sres = [sres; ar.sconstr(:, ar.qFit==1)];
    end
end
l = sum(res.^2);
if(nargout>1)
    g = res*sres;
end

function [c, ceq, gc, gceq] = confun(pTrial)
global ar
global fmincon_pTrial
if(sum(pTrial(:)~=fmincon_pTrial)>0)
    arChi2(ar.config.useSensis, pTrial)
    fmincon_pTrial = pTrial(:);
end
% Nonlinear inequality constraints
c = [];
% Nonlinear equality constraints
ceq = ar.constr;
if(nargout>2)
    gc = [];
    gceq = ar.sconstr(:, ar.qFit==1)';
end

function hessian = fmincon_hessianfcn(pTrial, lambda)
global ar
global fmincon_pTrial
if(sum(pTrial(:)~=fmincon_pTrial)>0)
    arChi2(ar.config.useSensis, pTrial)
    fmincon_pTrial = pTrial(:);
end
H = ar.sres(:, ar.qFit==1)'*ar.sres(:, ar.qFit==1);
Hconstr = zeros(size(H));
for jc = 1:length(ar.constr)
    Hconstr = Hconstr + lambda.eqnonlin(jc)*(ar.sconstr(jc, ar.qFit==1)'*ar.sconstr(jc, ar.qFit==1));
end
hessian = H + Hconstr;


% STRSCNE
function res = merit_fkt_STRSCNE(pTrial)
global ar
arChi2(ar.config.useSensis, pTrial)
res = [ar.res ar.constr]';

% derivatives for STRSCNE
function sres = merit_dfkt_STRSCNE(~)
global ar
if(ar.config.useSensis)
    sres = ar.sres(:, ar.qFit==1);
    if(~isempty(ar.sconstr))
        sres = [sres; ar.sconstr(:, ar.qFit==1)];
    end
end




function stop = arPlotFast(~,optimValues,state)

global ar

stop = false;

if(strcmp(state, 'iter'))
    if(ar.fit.iter_count>0)
        if((ar.chi2fit+ar.chi2constr) > (ar.fit.chi2_hist(ar.fit.iter_count) + ...
                ar.fit.constr_hist(ar.fit.iter_count)))
            return;
        end
    end
        
    ar.fit.chi2_hist(ar.fit.iter_count+1) = ar.chi2fit;
    ar.fit.constr_hist(ar.fit.iter_count+1) = ar.chi2constr;
    ar.fit.p_hist(ar.fit.iter_count+1,:) = ar.p;
    ar.fit.opti_hist(ar.fit.iter_count+1,:) = ar.firstorderopt;
    
    if(ar.config.optimizer == 5)
        if(isscalar(optimValues.mu))
            ar.fit.maxstepsize_hist(ar.fit.iter_count+1) = optimValues.mu;
        else
            ar.fit.maxstepsize_hist(ar.fit.iter_count+1) = det(optimValues.mu);
        end
        ar.fit.stepsize_hist(ar.fit.iter_count+1) = optimValues.normdp;
    else
        ar.fit.maxstepsize_hist(ar.fit.iter_count+1) = nan;
        if(ar.fit.iter_count>0)
            ar.fit.stepsize_hist(ar.fit.iter_count+1) = norm(ar.fit.p_hist(ar.fit.iter_count,:) - ar.p);
        end
    end
    
    ar.fit.iter_count = ar.fit.iter_count + 1;
    
    if(ar.config.showFitting)
        arPlot(false, true, false, true, true);
        drawnow;
    end
end

