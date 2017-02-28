% Fit model parameters to data using lsqnonlin (by default)
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
% Select different optimizers by editing ar.config.optimizer
%       1 - lsqnonlin (default)
%       2 - fmincon
%       4 - STRSCNE (Bellavia et al, A Scaled Trust Region Solver for Constrained Nonlinear Equations)
%       5 - arNLS (additional method choices under ar.config.optimizerStep; see help arNLS)
%       6 - fmincon_as_lsq
%       7 - arNLS with SR1 updates
%       8 - NL2SOL (Denis et al, Algorithm 573:  NL2SOL--An Adaptive Nonlinear Least-Squares)
%       9 - TRESNEI (B.Morini, M.Porcelli "TRESNEI, a Matlab trust-region solver for systems 
%       of nonlinear equalities and inequalities")
%      10 - Ceres (Sameer Agarwal and Keir Mierle and Others, Google Solver)
%      11 - lsqnonlin_repeated - repeated runs of lsqnonlin until convergence
%      12 - fminsearchbnd
%      13 - patternsearch
%      14 - patternsearch combined with fminsearchbnd
%      15 - particleswarm
%      16 - simulated annealing
%      17 - ga (geneticalgorithm)
%
%  Convergence of the optimization algorithm can be stored by setting
%  ar.config.logFitting = 1. Then interesing variables for each iteration
%  is stored in ar.fit.optimLog

function varargout = arFit(varargin)

global ar
global fit

if(nargin==0)
    qglobalar = true;
    silent = false;
else
    if(isstruct(varargin{1}))
        qglobalar = false;
        ar = varargin{1};
        if(nargin>1)
            varargin = varargin(2:end);
        else
            varargin = {};
        end
    else
        qglobalar = true;
    end
    
    if(~isempty(varargin))
        silent = varargin{1};
    else
        silent = false;
    end
end

if(~isfield(ar.config, 'optimizer'))
    ar.config.optimizer = 1;
end
if(~isfield(ar.config, 'optimizerStep'))
    ar.config.optimizerStep = 0;
end
if(~isfield(ar.config, 'showFitting'))
    ar.config.showFitting = 0;
end
if(~isfield(ar.config, 'logFitting'))
    ar.config.logFitting = 0;
end
ar.config.optim.OutputFcn = cell(0);
if(ar.config.logFitting)
    ar.config.optim.OutputFcn = [ar.config.optim.OutputFcn, {@arLogFitDetailed}];
elseif(isfield(ar,'fit'));
    if(isfield(ar.fit,'optimLog'))
        ar.fit = rmfield(ar.fit,'optimLog');
    end
end
if(ar.config.showFitting)
    ar.config.optim.OutputFcn = [ar.config.optim.OutputFcn, {@arPlotFast}];
end

if(ar.config.useSensis)
    ar.config.optim.Jacobian = 'on';
else
    ar.config.optim.Jacobian = 'off';
end

removeL1path = false;
if any(ar.type==3)
    if(~isfield(ar.config, 'l1trdog'))
        ar.config.l1trdog = 1;
        l1trdog(); % Modify trdog for L1 regularization
    end
    
    removeL1path = true;
    ar_path = fileparts(which('arInit.m'));
    addpath([ar_path '/l1/trdog'])
end

fit = struct([]);
fit(1).iter_count = 0;
fit.fevals = 0;

if(ar.config.logFitting)
    fit.chi2_hist = nan(1,ar.config.optim.MaxIter);
    fit.constr_hist = nan(1,ar.config.optim.MaxIter);
    fit.opti_hist = nan(1,ar.config.optim.MaxIter);
    fit.p_hist = nan(ar.config.optim.MaxIter,length(ar.p));
    fit.maxstepsize_hist = nan(1,ar.config.optim.MaxIter);
    fit.stepsize_hist = nan(1,ar.config.optim.MaxIter);
end

ar = arCalcMerit(ar, true, ar.p(ar.qFit==1));
chi2_old = arGetMerit('chi2fit');

if(sum(ar.qFit==1)<=0)
    error('No parameters are allowed to be fitted. Check ar.qFit.')
end

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
    % options.Algorithm = 'trust-region-reflective';
    options.SubproblemAlgorithm = 'cg';
    % options.Hessian = 'fin-diff-grads';
    options.Hessian = 'user-supplied';
    options.HessFcn = @fmincon_hessianfcn;
    % options2.InitBarrierParam = 1e+6;
    % options2.InitTrustRegionRadius = 1e-1;
    
    switch options.Algorithm
        case 'interior-point'
            myconfun = @confun;
        case 'trust-region-reflective'
            myconfun = [];
    end

    [pFit, ~, exitflag, output, lambda, grad] = ...
        fmincon(@merit_fkt_fmincon, ar.p(ar.qFit==1),[],[],[],[],lb,ub, ...
        myconfun,options);
    resnorm = merit_fkt(pFit);
    jac = [];
    fit.grad = grad;
    
% PSO
elseif(ar.config.optimizer == 3) 
    [pFit, ~, resnorm, exitflag, output, lambda, jac] = ...
        arFitPSO(lb, ub);
   
% STRSCNE
elseif(ar.config.optimizer == 4)
    
    if(~isempty(ar.config.optim.Display))
        silent = strcmp(ar.config.optim.Display,'iter')==0;
    else
        silent = 1;
    end
    if(~isempty(ar.config.optim.MaxIter))
        maxiter = ar.config.optim.MaxIter;
    else
        maxiter = 1e3;
    end
    if(~isempty(ar.config.optim.MaxFunEvals))
        maxfneval = ar.config.optim.MaxFunEvals;
    else
        maxfneval = 100*length(ar.p(ar.qFit==1));
    end
    delta = 1;
    warnreset = warning;
    warning('off','MATLAB:rankDeficientMatrix');
    [pFit, exitflag, output] = ...
        STRSCNE(ar.p(ar.qFit==1), @merit_fkt_STRSCNE, [-Inf,0], ...
        lb, ub, [maxiter,maxfneval,delta,~silent], @merit_dfkt_STRSCNE);
    warning(warnreset);
    resnorm = merit_fkt(pFit);
    jac = merit_dfkt_STRSCNE(pFit);
    lambda = [];
    % convert to lsqnonlin exitflag
    if(exitflag==1 || exitflag == 2)
        exitflag = 0;
    else
        exitflag = NaN;
    end

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
        [],options);% x l, u are transposed because TRESNEI.m uses column notation, opposed to row notation in D2D

    resnorm = merit_fkt(pFit);
    
% arNLS boosted by SR1 updates
elseif(ar.config.optimizer == 7)
    [pFit, ~, resnorm, exitflag, output, lambda, jac] = ...
        arNLS(@merit_fkt_sr1, ar.p(ar.qFit==1), lb, ub, ar.config.optim, ar.config.optimizerStep);

% NL2SOL
elseif(ar.config.optimizer == 8)
    if ~exist('mexnl2sol', 'file')
        compileNL2SOL;
    end
    [pFit, ~, resnorm, exitflag, output.iterations, lambda, jac] = ...
        mexnl2sol(@merit_fkt, ar.p(ar.qFit==1), lb, ub, ar.config.optim, 1);

% TRESNEI
elseif(ar.config.optimizer == 9)
    [pFit, exitflag, output, lambda, jac] = ...
        arTRESNEI(@merit_fkt, ar.p(ar.qFit==1), lb, ub, ar.config.optim);
    resnorm = merit_fkt(pFit);  
 
% Ceres
elseif(ar.config.optimizer == 10)
    if ~exist('ceresd2d', 'file')
         compileCeres;
    end
    [pFit, ~, ~, exitflag, output.iterations, jac, ceresexitmessage] = ...
        ceresd2d(@merit_fkt, ar.p(ar.qFit==1), lb, ub, ar.config.optimceres);
    resnorm = merit_fkt(pFit);
    lambda = [];
    fit.ceresexitmessage = ceresexitmessage;
    
% Repeated optimization with lsqnonlin
elseif(ar.config.optimizer == 11)
    [pFit, resnorm, res, exitflag, output, lambda, jac] = ...
        lsqnonlin(@merit_fkt, ar.p(ar.qFit==1), lb, ub, ar.config.optim);
    
    % Keep going until we really converged
    resnormOld = inf;
    while (((resnormOld-resnorm)/resnorm) > 1e-5)
        resnormOld = resnorm;
        [pFit, resnorm, res, exitflag, output, lambda, jac] = ...
            lsqnonlin(@merit_fkt, pFit, lb, ub, ar.config.optim);
    end
    
% fminsearchbnd
elseif(ar.config.optimizer == 12)
    options = optimset('fminsearch');
    if(~isempty(ar.config.optim.Display))
        options.Display = ar.config.optim.Display;
    end
    if(~isempty(ar.config.optim.TolX))
        options.TolX = ar.config.optim.TolX;
    end
    if(~isempty(ar.config.optim.TolFun))
        options.TolFun = ar.config.optim.TolFun;
    end    
    [pFit, ~, exitflag, output] = ...
        fminsearchbnd(@merit_fkt_chi2, ar.p(ar.qFit==1), lb, ub, options);
    resnorm = merit_fkt(pFit);
    lambda = [];
    jac = [];

% patternsearch    
elseif(ar.config.optimizer == 13)
    options = psoptimset;
    if(~isempty(ar.config.optim.Display))
        options.Display = ar.config.optim.Display;
    end
    [pFit, ~, exitflag, output] = ...
        patternsearch(@merit_fkt_chi2, ar.p(ar.qFit==1), [], [], [], [], lb, ub, [], options);
    resnorm = merit_fkt(pFit);
    lambda = [];
    jac = [];

% patternsearch hybrid
elseif(ar.config.optimizer == 14)
    options = psoptimset;
    if(~isempty(ar.config.optim.Display))
        options.Display = ar.config.optim.Display;
    end
    pFit = ...
        patternsearch(@merit_fkt_chi2, ar.p(ar.qFit==1), [], [], [], [], lb, ub, [], options);

    options = optimset('fminsearch');
    if(~isempty(ar.config.optim.Display))
        options.Display = ar.config.optim.Display;
    end
    if(~isempty(ar.config.optim.TolX))
        options.TolX = ar.config.optim.TolX;
    end
    if(~isempty(ar.config.optim.TolFun))
        options.TolFun = ar.config.optim.TolFun;
    end 

    [pFit, ~, exitflag, output] = ...
        fminsearchbnd(@merit_fkt_chi2, pFit, lb, ub, options);
    resnorm = merit_fkt(pFit);
    lambda = [];
    jac = [];   

% particleswarm
elseif(ar.config.optimizer == 15)
    options = optimoptions('particleswarm');
    if(~isempty(ar.config.optim.Display))
        options.Display = ar.config.optim.Display;
    end
    [pFit, ~, exitflag, output] = ...
        particleswarm(@merit_fkt_chi2, length(ar.p(ar.qFit==1)), lb, ub, options);
    resnorm = merit_fkt(pFit);
    lambda = [];
    jac = [];
    
% simulated annealing
elseif(ar.config.optimizer == 16)
    options = saoptimset('simulannealbnd');
    if(~isempty(ar.config.optim.Display))
        options.Display = ar.config.optim.Display;
    end
    [pFit, ~, exitflag, output] = ...
        simulannealbnd(@merit_fkt_chi2, ar.p(ar.qFit==1), lb, ub, options);
    resnorm = merit_fkt(pFit);
    lambda = [];
    jac = [];
      
% ga
elseif(ar.config.optimizer == 17)
    options = gaoptimset;
    if(~isempty(ar.config.optim.Display))
        options.Display = ar.config.optim.Display;
    end
    [pFit, ~, exitflag, output] = ...
        ga(@merit_fkt_chi2, sum(ar.qFit==1), [], [], [], [], lb, ub, [], options);
    resnorm = merit_fkt(pFit);
    lambda = [];
    jac = [];
    output.iterations = output.generations;

else
    error('ar.config.optimizer invalid');    
end

if(isfield(ar, 'ms_count_snips') && ar.ms_count_snips>0)
    if(max(ar.ms_violation) > ar.ms_threshold)
        arFprintf(1, 'Multiple Shooting: continuity constains violated %e > %e\n', max(ar.ms_violation), ar.ms_threshold);
    end
end

if removeL1path
    rmpath([ar_path '/l1/trdog'])
end

ar.p(ar.qFit==1) = pFit;
ar = arCalcMerit(ar, true, ar.p(ar.qFit==1));

fit.exitflag = exitflag;
fit.output = output;
fit.iter = output.iterations;
fit.chi2 = arGetMerit('chi2fit');
fit.lambda = lambda;
fit.qFit = ar.qFit;
fit.res = resnorm;
fit.sres = full(jac);
fit.improve = chi2_old - fit.chi2;
if(isfield(ar,'fit') && isfield(ar.fit,'optimLog'))
    fit.optimLog = ar.fit.optimLog;
end

ar.fit = fit;

if(~silent || exitflag < 1)
    arFitPrint;
end

if(~silent)
    ar = arCalcMerit(ar, true);
end

if(nargout>0 && ~qglobalar)
    varargout{1} = ar;
else
    varargout = {};
end


% lsqnonlin and arNLS
function [res, sres] = merit_fkt(pTrial)
global ar

% Only compute sensis when requested
if ( isfield( ar.config, 'sensiSkip' ) )
    sensiskip = ar.config.sensiSkip;
else
    sensiskip = false;
end
sensi = ar.config.useSensis && (~sensiskip || (nargout > 1));

arCalcMerit(sensi, pTrial);
arLogFit(ar);
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

np = sum(ar.qFit==1);
if ( numel(res) < np )
    tres = zeros(1,np);
    tres(1:length(res)) = res;
    tsres = zeros(np);
    tsres(1:length(res), 1:np) = sres;
    
    res = tres;
    sres = tsres;
end

% arNLS boosted by SR1 updates
function [res, sres, H, ssres] = merit_fkt_sr1(p, pc, ~, sresc, ssresc)
global ar
arCalcMerit(ar.config.useSensis, p);
arLogFit(ar);
res = [ar.res ar.constr];
if(nargout>1 && ar.config.useSensis)
    sres = [];
    if(~isempty(ar.sres))
        sres = ar.sres(:, ar.qFit==1);
    end
    if(~isempty(ar.sconstr))
        sres = [sres; ar.sconstr(:, ar.qFit==1)];
    end
    
    if(nargout>2)
        if(nargin<2)
            H = 2*(sres'*sres);
            ssres = zeros([length(p) length(p) length(res)]);
            return;
        end
        
        % SR1 update on residuals
        ssres = nan(size(ssresc));
        for j=1:length(res)
            ssres(:,:,j) = sr1_update(ssresc(:,:,j), pc, sresc(j,:), p, sres(j,:));
            % ssres(j,:,:) = sr1_update(ssresc(j,:,:), pc, 2*sresc(j,:), p, 2*sres(j,:));
            % ssres(:,:,j) = sr1_update(ssresc(:,:,j), pc, 0.5*sresc(j,:), p, 0.5*sres(j,:));
            % ssres(j,:,:) = sr1_update(ssresc(:,:,j), pc, -2*sresc(j,:), p, -2*sres(j,:));
            % ssres(j,:,:) = sr1_update(ssresc(:,:,j), pc, -sresc(j,:), p, -sres(j,:));
        end
        
%         figure(1);
%         subplot(1,2,1)
%         imagesc(2*(sres'*sres));
%         subplot(1,2,2)
%         imagesc(squeeze(sum(bsxfun(@times, res', ssres),1)));
        
        H = 2*(sres'*sres) + sum(bsxfun(@times, shiftdim(res,-1), ssres),3);
    end
end

% SR1 updates
function H = sr1_update(H, pC, gC, pO, gO)
sk = transpose(pC - pO);
yk = transpose(gC - gO);
rk = yk - H*sk;
c1 = 0.5;
if(abs(rk'*sk) > c1*norm(rk)*norm(sk))
    H = H + (rk*rk') / (rk'*sk);
end

% fmincon
function [l, g, H] = merit_fkt_fmincon(pTrial)
global ar
arCalcMerit(ar.config.useSensis, pTrial);
arLogFit(ar);
l = sum(ar.res.^2);
if(nargout>1)
    g = ar.res*ar.sres(:, ar.qFit==1);
end
if(nargout>2)
    type3_ind = ar.type == 3;
    type3_ind = type3_ind(ar.qFit==1);
    
    H = ar.sres(:, ar.qFit==1)'*ar.sres(:, ar.qFit==1);
    H(type3_ind,type3_ind) = H(type3_ind,type3_ind) .* ~eye(sum(type3_ind));
end

% fmincon as lsq
function [l, g] = merit_fkt_fmincon_lsq(pTrial)
global ar
arCalcMerit(ar.config.useSensis, pTrial);
arLogFit(ar);
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
arCalcMerit(ar.config.useSensis, pTrial);
arLogFit(ar);
% Nonlinear inequality constraints
c = [];
% Nonlinear equality constraints
ceq = ar.constr;
if(nargout>2)
    gc = [];
    if(~isempty(ar.sconstr))
        gceq = ar.sconstr(:, ar.qFit==1)';
    else
        gceq = [];
    end
end

function hessian = fmincon_hessianfcn(pTrial, lambda)
global ar
arCalcMerit(ar.config.useSensis, pTrial);
arLogFit(ar);
H = ar.sres(:, ar.qFit==1)'*ar.sres(:, ar.qFit==1);
Hconstr = zeros(size(H));
for jc = 1:length(ar.constr)
    Hconstr = Hconstr + lambda.eqnonlin(jc)*(ar.sconstr(jc, ar.qFit==1)'*ar.sconstr(jc, ar.qFit==1));
end
hessian = H + Hconstr;

% STRSCNE
function res = merit_fkt_STRSCNE(pTrial)
global ar
arCalcMerit(ar.config.useSensis, pTrial);
arLogFit(ar);
res = [ar.res ar.constr]';

% derivatives for STRSCNE
function sres = merit_dfkt_STRSCNE(pTrial)
global ar
if(ar.config.useSensis)
    arCalcMerit(ar.config.useSensis, pTrial);
    sres = ar.sres(:, ar.qFit==1);
    if(~isempty(ar.sconstr))
        sres = [sres; ar.sconstr(:, ar.qFit==1)];
    end
end

% fminsearch, particleswarm, simulannealbnd, patternsearch
function chi2 = merit_fkt_chi2(pTrial)
global ar

try
    arCalcMerit(false, pTrial);
    arLogFit(ar);
    chi2 = sum([ar.res ar.constr].^2);
catch
    % workaround for particleswarm
    chi2 = rand(1)*1e23;
end

% plot fitting
function stop = arPlotFast(~,~,state)
global ar
stop = false;
if(strcmp(state, 'iter'))    
    if(ar.config.showFitting)
        % note that the ~ar.config.useNewPlots is a temporary fix. arPlot2,
        % which is called when ar.config.useNewPlots is set to 1 does not 
        % properly redraw -2 log(L) labels when forced to do a fastplot
        % TO DO: The proper solution would be to make sure that fastplot 
        % also updates the -2 log(L) text label
        arPlot(false, ~ar.config.useNewPlots, true, true, true);
        drawnow;
    end
end

% log fitting
function arLogFit(ar)

global fit

fit.fevals = fit.fevals + 1;

if(fit.iter_count>0)
    if(arGetMerit('chi2all') > (fit.chi2_hist(fit.iter_count) + ...
            fit.constr_hist(fit.iter_count)))
        return;
    end
end

fit.chi2_hist(fit.iter_count+1)   = arGetMerit('chi2fit');
fit.constr_hist(fit.iter_count+1) = arGetMerit('chi2constr');
fit.p_hist(fit.iter_count+1,:) = ar.p;
fit.opti_hist(fit.iter_count+1,:) = ar.firstorderopt;
fit.maxstepsize_hist(fit.iter_count+1) = nan;
if(fit.iter_count>0)
    q = fit.p_hist(fit.iter_count,:);
    fit.stepsize_hist(fit.iter_count+1) = norm(q(:) - ar.p(:));
end
fit.iter_count = fit.iter_count + 1;


function stop = arLogFitDetailed(~,optimValues,state)
stop = false;
global ar

if(ar.config.optimizer ==1)
    fn = {'iteration','funccount','stepsize','firstorderopt','cgiterations','positivedefinite','ratio','degenerate','trustregionradius','resnorm','gradient_norm'};
end

optimValues.gradient_norm = norm(optimValues.gradient);
% optimValues.x = x;
% optimValues = rmfield(optimValues,'gradient');
% optimValues = rmfield(optimValues,'residual');
switch state
    case 'init'
        ar.fit.optimLog.values = NaN(ar.config.optim.MaxIter,length(fn));
        ar.fit.optimLog.labels = fn;
    case 'iter'
        for i=1:length(fn)
            ar.fit.optimLog.values(optimValues.funccount,i) = optimValues.(fn{i});
        end
    case 'done'
        ar.fit.optimLog.values = ar.fit.optimLog.values(sum(~isnan(ar.fit.optimLog.values),2)>0,:);
    otherwise
end

