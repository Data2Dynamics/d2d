% Initialize Profile Likelihood Exploit
%
% arPLEInit(force, breakon_point, mode)
%   force:              [true] = exising PLEs are deleted
%   breakon_point:      [false] = calc simultaneous CIs, true = calc pointwise CIs
%   mode:               [1] = direct step, 2 = progressive step

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

pleInit(ar.p, ar.qFit==1, ar.lb, ar.ub, ar.qLog10, @arPLEIntegrate, @arPLEMerit, ...
    @arPLEDiffMerit, @arPLEFit, @arPLESetOptim, ar.pLabel, 1-ar.ppl.alpha_level, force);

global pleGlobals;

pleGlobals.violations = @arPLEMeritViolations;
pleGlobals.priors = @arPLEPrior;
pleGlobals.priorsAll = @arPLEPriorAll;

if(breakon_point)
    pleGlobals.breakon_point = true;
    pleGlobals.plot_simu = false;
else
    pleGlobals.breakon_point = false;
    pleGlobals.plot_simu = true;
end

if(ar.config.fiterrors == 1)
    pleGlobals.ylabel = '-2 log(PL)';
end

if(isfield(ar, 'pTrue'))
    pleGlobals.p_true = ar.pTrue;
end
if(mode==1)
    pleGlobals.initstep_fkt = @pleInitStepDirect;
    pleGlobals.mode = 1;
elseif(mode==2)
    pleGlobals.initstep_fkt = @pleInitStep;
    pleGlobals.mode = 2;
elseif(mode==3)
    pleGlobals.initstep_fkt = @pleInitStepLinear;
    pleGlobals.mode = 3;
elseif(mode==4)
    pleGlobals.initstep_fkt = @pleInitStepComposite;
    pleGlobals.mode = 4;
elseif(mode==5)
    pleGlobals.initstep_fkt = @pleInitStepDirect2;
    pleGlobals.mode = 5;
end

pleGlobals.savePath = [arSave '/PLE'];
ar.ple_errors = [];

function arPLEIntegrate(p)
global ar
try
    arChi2(false, p(ar.qFit==1));
catch exception
    ar.ple_errors(end+1,:) = ar.p;
    fprintf('ERROR INTEGRATOR (#%i): %s\n', size(ar.ple_errors,1), exception.message);
    rethrow(exception)  
end

function chi2 = arPLEMerit
global ar
if(ar.config.fiterrors == 1)
    chi2 = 2*ar.ndata*log(sqrt(2*pi)) + ar.chi2fit;
else
    chi2 = ar.chi2;
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
global pleGlobals
ar.p = p + 0;
arChi2(false);
pleGlobals.p = p+0;
pleGlobals.chi2 = arPLEMerit+0;



