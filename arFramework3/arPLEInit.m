% Initialize Profile Likelihood Exploit
%
% arPLEInit(breakon_point, mode)
%   breakon_point:      [false] = calc simultaneous CIs, true = calc pointwise CIs
%   mode:               [1] = direct step, 2 = progressive step

function arPLEInit(breakon_point, mode)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('breakon_point', 'var'))
    breakon_point = true;
end
if(~exist('mode', 'var'))
    mode = 1;
end

pleInit(ar.p, ar.qFit==1, ar.lb, ar.ub, ar.qLog10, @arPLEIntegrate, @arPLEMerit, ...
    @arPLEDiffMerit, @arPLEFit, @arPLESetOptim, ar.pLabel, 1-ar.ppl.alpha_level);

global pleGlobals;

pleGlobals.violations = @arPLEMeritViolations;
pleGlobals.priors = @arPLEPrior;

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
chi2 = ar.chi2ss + ar.chi2prior;
if(chi2==0)
    chi2 = 0;
end

function chi2 = arPLEPrior(jk)
global ar
if(ar.type(jk)==1)
    chi2 = ((ar.mean(jk)-ar.p(jk))./ar.std(jk))^2;
else
    chi2 = 0;
end

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
    gradient(ar.qFit==1) = 2 * (ar.fit.res * ar.fit.sres);
    
    if(nargin==1)
        ar.qFit = qFitReset;
    end
    p = ar.p;
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



