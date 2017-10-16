% Initialize Profile Likelihood Exploit
%
% arPLEInit(force, breakon_point, mode)
%   force:              [true] = exising PLEs are deleted
%   breakon_point:      [false] = calc simultaneous CIs, true = calc pointwise CIs
%   mode:               [1] = direct step, 2 = progressive step
% 
%   The profile likelihood calculation by the functions ple* was intended
%   as running independent of D2D, i.e. it was intendent to be also used by
%   other tools. 
% 
%   This function extracts info from the global variable ar and calls
%   pleInit.

function arPLEInit(force, breakon_point, mode)

global ar
if(~exist('force','var') || isempty(force))
    force = true; % existing PLEs are removed
end

if(isempty(ar) )
    error('please initialize by arInit')
end

if(~exist('breakon_point', 'var') || isempty(breakon_point))
    breakon_point = true;
end
if(~exist('mode', 'var') || isempty(mode))
    mode = 1;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

alpha = 0.95;

%     pleInit(ar.p, ar.qFit==1, ar.lb, ar.ub, ar.qLog10, @arPLEIntegrate, @arPLEMerit, ...
%     @arPLEDiffMerit, @arPLEFit, @arPLESetOptim, ar.pLabel, 1-ar.ppl.alpha_level, force);

%% Beginning of old function pleInit.m

if ~isfield(ar,'ple') || ~isstruct(ar.ple) || force
    ar.ple = struct('p', ar.p, ...
        'integrate_fkt', @arPLEIntegrate, ...
        'merit_fkt', @arPLEMerit, ...
        'diffmerit_fkt', @arPLEDiffMerit, ...
        'fit_fkt', @arPLEFit, ...
        'setoptim_fkt', @arPLESetOptim, ...
        'initstep_fkt', @pleInitStepDirect, ...  % is overwritten in arPLEInit, if mode is provided
        'dof', sum(ar.qFit==1), ...
        'dof_point', 1,... 
        'attempts', 4',...
        'showCalculation', true,...
        'ylabel', '\chi^2_{PL}',...
        'breakon_point', true,... %false',...
        'plot_point', true,...
        'plot_simu', false,... %true',...
        'dist_thres', 0.01,...
        'grad_thres', 1,...
        'closetobound', 0.001,...
        'allowbetteroptimum', false);
end


% step sizes
ar.ple.samplesize = 50 * ones(size(ar.ple.p));
ar.ple.relchi2stepincrease = 0.1 * ones(size(ar.ple.p));
ar.ple.maxstepsize = (ar.ub-ar.lb)./ar.ple.samplesize;
ar.ple.minstepsize = ones(size(ar.ple.maxstepsize))*1e-3;
ar.ple.breakonlb = false(size(ar.ple.p));
ar.ple.breakonub = false(size(ar.ple.p));


ar.ple.alpha_level = 1-alpha;
ar.ple.dchi2 = chi2inv(1-ar.ple.alpha_level, ar.ple.dof);
ar.ple.dchi2_point = chi2inv(1-ar.ple.alpha_level, ar.ple.dof_point);

% magic factors
ar.ple.chi2_strID_ratio = 1e-1;
ar.ple.svd_threshold = 1e-6; % SVD regulatization threshold (NR: chapter 15.4)
ar.ple.optimset_tol = 1e-1;

% % labels
% if(~exist('p_labels', 'var'))
%     p_labels = {};
%     for j=1:length(ar.ple.p)
%         p_labels{j} = sprintf('p%02i', j); %#ok<AGROW>
%     end
% end
% ar.ple.p_labels = p_labels;
ar.ple.p_labels = ar.pLabel;

ar.ple.conf_lb = nan(1,length(ar.ple.p));
ar.ple.conf_ub = nan(1,length(ar.ple.p));
ar.ple.conf_lb_point = nan(1,length(ar.ple.p));
ar.ple.conf_ub_point = nan(1,length(ar.ple.p));
ar.ple.conf_rel = nan(1,length(ar.ple.p));
ar.ple.conf_rel_point = nan(1,length(ar.ple.p));

ar.ple.IDstatus = nan(1,length(ar.ple.p));
ar.ple.IDstatus_point = nan(1,length(ar.ple.p));
ar.ple.IDlabel = {'', 'pra.nID', 'str.nID', 'single str.nID'};
ar.ple.savePath = ['PLE-' datestr(now, 30)];

if(~isfield(ar.ple,'ps') || force)
    ar.ple.ps = {};
end
if(~isfield(ar.ple,'psinit') || force)
    ar.ple.psinit = {};
end
if(~isfield(ar.ple,'psinitstep') || force)
    ar.ple.psinitstep = {};
end
if(~isfield(ar.ple,'chi2s') || force)
    ar.ple.chi2s = {};
end
if(~isfield(ar.ple,'chi2sinit') || force)
    ar.ple.chi2sinit = {};
end

ar.ple.finished = 0;

%% End of old function pleInit


ar.ple.violations = @arPLEMeritViolations;
ar.ple.priors = @arPLEPrior;
ar.ple.priorsAll = @arPLEPriorAll;

if(breakon_point)
    ar.ple.breakon_point = true;
    ar.ple.plot_simu = false;
else
    ar.ple.breakon_point = false;
    ar.ple.plot_simu = true;
end

if( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors == 1) || ...
        (ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0) )
    ar.ple.ylabel = '-2 log(PL)';
end

if(isfield(ar, 'pTrue'))
    ar.ple.p_true = ar.pTrue;
end
if(mode==1)
    ar.ple.initstep_fkt = @pleInitStepDirect;
    ar.ple.mode = 1;
elseif(mode==2)
    ar.ple.initstep_fkt = @pleInitStep;
    ar.ple.mode = 2;
elseif(mode==3)
    ar.ple.initstep_fkt = @pleInitStepLinear;
    ar.ple.mode = 3;
elseif(mode==4)
    ar.ple.initstep_fkt = @pleInitStepComposite;
    ar.ple.mode = 4;
elseif(mode==5)
    ar.ple.initstep_fkt = @pleInitStepDirect2;
    ar.ple.mode = 5;
end

ar.ple.savePath = [arSave '/PLE'];
ar.ple.errors = [];



function arPLEIntegrate(p)
global ar
try
    arCalcMerit(false, p(ar.qFit==1));
catch exception
    if ( ~isfield( ar, 'ple_errors' ) )
        ar.ple.errors = ar.p;
    end
    ar.ple.errors(end+1,:) = ar.p;
    fprintf('ERROR INTEGRATOR (#%i): %s\n', size(ar.ple.errors,1), exception.message);
    rethrow(exception)  
end

function chi2 = arPLEMerit
global ar
if ar.config.fiterrors == 1 || (ar.config.fiterrors == 0 && sum(ar.qFit(ar.qError==1)<2)>0)
    chi2 = 2*ar.ndata*log(sqrt(2*pi)) + arGetMerit('chi2fit');
% elseif(ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0)
%     chi2 = 2*ar.ndata_err*log(sqrt(2*pi)) + arGetMerit('chi2fit');
else
    chi2 = arGetMerit('chi2');
end

function chi2 = arPLEMeritViolations
global ar
chi2 = ar.chi2constr + ar.chi2prior;
if(chi2==0)
    chi2 = 0;
end

function chi2 = arPLEPrior(jk)
global ar
if(ar.type(jk)==1)
    chi2 = ((ar.mean(jk)-ar.p(jk))./ar.std(jk))^2;
elseif(ar.type(jk)==3)
    chi2 = abs((ar.mean(jk)-ar.p(jk))./ar.std(jk));
else
    chi2 = 0;
end

function chi2all = arPLEPriorAll
global ar
type1 = find(ar.type==1);
chi21 = ((ar.mean(type1)-ar.p(type1))./ar.std(type1)).^2;

type3 = find(ar.type==3);
chi23 = abs((ar.mean(type3)-ar.p(type3))./ar.std(type3));

chi2all = sum(chi21) + sum(chi23);
    

function [beta, alpha] = arPLEDiffMerit
global ar
beta = transpose(-ar.res*ar.sres(:,ar.qFit==1));
alpha = ar.sres(:,ar.qFit==1)'*ar.sres(:,ar.qFit==1);


function [p, gradient] = arPLEFit(jk)
global ar
if(nargin==1)
    qFitReset = ar.qFit;
    ar.qFit(jk) = 0;
end
try
    arFit(true);
    
    gradient = nan(1, length(ar.p));
    if ~isempty(ar.fit.sres)
        gradient(ar.qFit==1) = 2 * (ar.fit.res * ar.fit.sres);
    else
        gradient(ar.qFit==1) = ar.fit.grad;
    end
    
    if(nargin==1)
        ar.qFit = qFitReset;
    end
    p = ar.p;
    arPLEIntegrate(ar.p);  % Recalculate objective function, with FitErrorCorrection calculated with qFitReset
catch exception
    disp(['ERROR FIT: ' exception.message]);
    if(nargin==1)
        ar.qFit = qFitReset;
    end
    rethrow(exception)  
end

function arPLESetOptim(p)
global ar

ar.p = p + 0;
arCalcMerit(false);
ar.ple.p = p+0;
ar.ple.merit = arPLEMerit+0;
