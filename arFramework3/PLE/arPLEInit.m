% arPLEInit([force], [breakon_point], [mode],[autosave])
%
% Initialize Profile Likelihood Exploit
%   force           existing PLEs are deleted  [true]
%   breakon_point   calc simultaneous (false) or pointwise (true) CIs  [true]
%   mode            Specifies step choice algorithm.
%                       (1): Changing only profile parameters
%                       (2): Use information of direction of previous step
%   autosave        Indicate whether Profiles should be saved automatically
%                       [true]
% 
% The profile likelihood calculation by the functions ple* was intended
% as running independent of D2D, i.e. it was intended to be also used by
% other tools. 
%
% Mode 2 proposes steps more efficiently and thus performs better in most 
% scenarios. In cases with many local optima or poorly specified integration
% tolerances, mode 1 may sample profiles more accurately, although not
% sparse enough to reach the threshold quickly.
%
% See also: ple

function arPLEInit(force, breakon_point, mode,autosave)

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
    mode = 2;  % Before March, 20th 2020, mode=1 was default.
end
if(~exist('autosave', 'var') || isempty(autosave))
    autosave= true;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

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
        'breakon_point', true,... 
        'plot_point', true,...
        'plot_simu', false,... 
        'dist_thres', 0.01,... % threshold for values of other parameters to assess whether they are at bounds (evaluated in plePlot* => plotting with 'ko')
        'grad_thres', 0.01,... % threshold for gradient of log-likelihood when assessing whether bounds of other parameters constrain fitting (evaluated in plePlot* => plotting with 'ko'). Before it was 1, unknown why, 1 seems to be too large
        'closetobound', 0.001,... 
        'continuousSave', false,...
        'allowbetteroptimum', false,...
        'autosave',autosave,...
        'usesensis',false);
end


%Configurations for profile steps
ar.ple.samplesize = 100 * ones(size(ar.ple.p));
ar.ple.relchi2stepincrease = 0.1 * ones(size(ar.ple.p));
ar.ple.maxstepsize = 2*(ar.ub-ar.lb)./ar.ple.samplesize;
ar.ple.minstepsize = ones(size(ar.ple.maxstepsize))*5*1e-4;
ar.ple.breakonlb = false(size(ar.ple.p));
ar.ple.breakonub = false(size(ar.ple.p));
ar.ple.stepfaktor = 1.5*ones(size(ar.ple.p)); 

%Chosen confidence level. alpha and chi2 need to both be set separately if they
%are changed.
ar.ple.alpha_level = 0.05;
ar.ple.dchi2 = arChi2inv(1-ar.ple.alpha_level, ar.ple.dof);
ar.ple.dchi2_point = arChi2inv(1-ar.ple.alpha_level, ar.ple.dof_point);

ar.ple.chi2_strID_ratio = 1e-1; 
%relative distance to threshold which is accepted as structurally NI
ar.ple.optimset_tol = 1e-1;
%absolute difference between initial point and new minimum required
%to accept the new minimum

ar.ple.p_labels = ar.pLabel;
ar.ple.lb = ar.lb;
ar.ple.ub = ar.ub;
ar.ple.qLog10 = ar.qLog10;
ar.ple.conf_lb = nan(1,length(ar.ple.p));
ar.ple.conf_ub = nan(1,length(ar.ple.p));
ar.ple.conf_lb_point = nan(1,length(ar.ple.p));
ar.ple.conf_ub_point = nan(1,length(ar.ple.p));
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

if ar.config.fiterrors==1 || (ar.config.fiterrors==0 && sum(ar.qFit==1 & ar.qError==1)>0)
    ar.ple.ylabel = '-2 log(PL)';
end

if(isfield(ar, 'pTrue'))
    ar.ple.p_true = ar.pTrue;
end

%Specifies Profile-Step-Choice-Algorithm
if(mode==1)
    ar.ple.initstep_fkt = @pleInitStepDirect;
    ar.ple.mode = 1;
elseif(mode==2)
    ar.ple.initstep_fkt = @pleInitStepInertia;
    ar.ple.mode = 2;
end

if (autosave)
    ar.ple.savePath = [arSave '/PLE'];
end
ar.ple.errors = [];

%Next Section defines functions with handles stored in the ple-struct.
function arPLEIntegrate(p)
global ar
try
    arCalcMerit(ar.ple.usesensis, p(ar.qFit==1));
    %Using sensitivitites in this function is a bit more time-intensive and 
    %makes no difference in most cases, but not using them may lead to wrongly
    %sampled profiles in some models.
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
    try
        switch ar.L1subtype(jk)
            case 1
                chi2 = abs((ar.mean(jk)-ar.p(jk))./ar.std(jk));
            case 2
                chi2 = abs(ar.lnuweights(jk).*(ar.mean(jk)-ar.p(jk))./ar.std(jk));
            case 3
                chi2 = abs((ar.mean(jk)-ar.p(jk)).^ar.expo(jk)./ar.std(jk));
            case 4
                chi2 = abs((ar.mean(jk)-ar.p(jk))./ar.std(jk).*(1-ar.alpha(jk))) ...
                    + ar.alpha(jk)./ar.std(jk) .* (ar.mean(jk) - ar.p(jk)).^ 2;
        end
    catch ME
        if strcmpi(ME.identifier,'MATLAB:UndefinedFunction')
            fprintf('Please provide a subtype and the corresponding parameters')
        end
        chi2 = 0;
    end
elseif(ar.type(jk)==5)
    bk = (ar.grplas.grouping == ar.grplas.grouping(jk) & ar.type == 5);
    chi2 = sqrt(((ar.mean(bk)-ar.p(bk))./ar.std(bk)) ...
        * ar.grplas.A(bk,bk) ...
        * ((ar.mean(bk)-ar.p(bk))./ar.std(bk))');
else
    chi2 = 0;
end

function chi2all = arPLEPriorAll
global ar
type1 = find(ar.type==1);
chi21 = ((ar.mean(type1)-ar.p(type1))./ar.std(type1)).^2;

chi23 = 0;
t3st1 = (ar.type == 3 & ar.L1subtype == 1);
t3st2 = (ar.type == 3 & ar.L1subtype == 2);
t3st3 = (ar.type == 3 & ar.L1subtype == 3);
t3st4 = (ar.type == 3 & ar.L1subtype == 4);
if any(t3st1)
    chi23(t3st1) = abs((ar.mean(t3st1)-ar.p(t3st1))./ar.std(t3st1));
end
if any(t3st2)
    chi23(t3st2) = abs(ar.lnuweights(t3st2).*(ar.mean(t3st2)-ar.p(t3st2))./ar.std(t3st2));
end
if any(t3st3)
    chi23(t3st3) = abs((ar.mean(t3st3)-ar.p(t3st3)).^ar.expo(t3st3)./ar.std(t3st3));
end

if any(t3st4)
    chi23(t3st4) = ar.alpha(t3st4) ./ ar.std(t3st4) ...
        .* abs(ar.mean(t3st4) - ar.p(t3st4)) + (1-ar.alpha(t3st4)) ./ ar.std(t3st4) ...
        .* (ar.mean(t3st4) - ar.p(t3st4)) .^ 2;
end

if any(ar.type == 5)
    chi25 = zeros(size(ar.grplas.groups));
    
    for i = 1:length(ar.grplas.groups)
        bk = (ar.grplas.grouping == ar.grplas.groups(i) & ar.type == 5);
        chi25(i) = sqrt(((ar.mean(bk)-ar.p(bk))./ar.std(bk)) ...
            * ar.grplas.A(bk,bk) ...
            * ((ar.mean(bk)-ar.p(bk))./ar.std(bk))');
    end
else
    chi25 = 0;
end


chi2all = sum(chi21) + sum(chi23) + sum(chi25);
    

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
    arPLEIntegrate(ar.p);  
    % Recalculate objective function, with FitErrorCorrection calculated with qFitReset
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
