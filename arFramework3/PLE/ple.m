% ple([i], [samplesize], [relchi2stepincrease], [maxstepsize], [minstepsize], [breakonbounds], [doLeftRightBranch])
% 
% Profile Likelihood Exploit
%
%   i                     parameter index or vector of par indices [all if not specified]
%                         alternatively the name of a parameter or a cell of
%                         names can be provided. If one parameter name is
%                         provided which is not in ar.ple.p_labels, then the string
%                         is interpreted as a regular expression.
%   samplesize            number of sampling steps in each direction     [100]
%   relchi2stepincrease   percentage chi^2 increase of a step wrt bound  [0.1]
%   maxstepsize           maximum size of a step              [(ar.ub-ar.lb)/50]
%   minstepsize           minumum size of a step              [5*e-4]
%   breakonlb             stop if hit lb                      [false]
%   breakonub             stop if hit ub                      [false]
%                         The "breakon" parameters are only relevant
%                         for the non-profile parameters. The profile will
%                         break regardless for the profile parameter if it
%                         reaches the boundary.
%   doLeftRightBranch     calculate left/right branch         [true true]
% 
% The profile likelihood calculation by the functions ple* was intended
% as running independent of D2D, i.e. it was intended to be also used by
% other tools. Therefore, the function does not use info in global ar and stores the
% results an a separate variable.
%
% Before usage, workspace has to be initialized.
%
% See also: arPLEInit, pleExtend, pleSmooth

function ple(jk, samplesize, relchi2stepincrease, ...
    maxstepsize, minstepsize, breakonlb, breakonub, doLeftRightBranch)

global ar

if(~isfield(ar,'ple') || isempty(ar.ple))
    error('PLE ERROR: please initialize')
end 

% Assign the vector of relevant parameter indices to jk:
if(~exist('jk','var') || isempty(jk))
    jk = find(ar.qFit==1);
elseif(ischar(jk))
    tref = strmatch(jk,ar.ple.p_labels,'exact');
    if(isempty(tref))
        jk = find(~cellfun(@isempty,regexp(ar.ple.p_labels,jk)));
    else 
        jk = tref;
    end
    
    if isempty(jk)
        disp(['Pattern ''',jk,''' not found in ar.ple.p_labels']);
        return;
    end
elseif(iscell(jk)) % cell of pLabels
    [~,jk] = intersect(ar.ple.p_labels,jk);
elseif(isnumeric(jk))
else
    error('Argument has to be a string or an array of indices.')
end

ar.ple.finished = 0;
%Used to check whether current profile has run to its conclusion.

if(exist('samplesize', 'var') && ~isempty(samplesize))
    ar.ple.samplesize(jk) = samplesize;
end
if(exist('relchi2stepincrease', 'var') && ~isempty(relchi2stepincrease))
    ar.ple.relchi2stepincrease(jk) = relchi2stepincrease;
end
if(exist('maxstepsize', 'var') && ~isempty(maxstepsize))
    ar.ple.maxstepsize(jk) = maxstepsize;
end
if(exist('minstepsize', 'var') && ~isempty(minstepsize))
    ar.ple.minstepsize(jk) = minstepsize;
end
if(exist('breakonlb', 'var') && ~isempty(breakonlb))
    ar.ple.breakonlb(jk) = breakonlb;
end
if(exist('breakonub', 'var') && ~isempty(breakonub))
    ar.ple.breakonub(jk) = breakonub;
end
if(~isfield(ar.ple, 'showCalculation'))
    ar.ple.showCalculation = true;
    % Field showCalculation can be used to turn off real-time plotting. Doing
    % so reduces computation time. 
end
if(~exist('doLeftRightBranch','var') || isempty(doLeftRightBranch))
    doLeftRightBranch = [true true];
end

% Specifications on what to do if there are more than 1 parameter. The
% function ple will evaluate itself for every parameter separately.
if(nargin<1)
    fprintf('PLE for %i parameters ...\n', sum(ar.qFit==1))
    jindex = find(ar.qFit);
    do_plotting = ar.ple.showCalculation;
    ar.ple.showCalculation = false;
    for j=1:length(jindex)
        ple(jindex(j));
        if(do_plotting)
            plePlotMulti;
        end
    end
    ar.ple.showCalculation = do_plotting;
    return
elseif(length(jk)>1)
    fprintf('PLE for %i parameters ...\n', length(jk))
    do_plotting = ar.ple.showCalculation;
    ar.ple.showCalculation = false;
    for j=1:length(jk)
        ple(jk(j));
        if(do_plotting)
            plePlotMulti;
        end
    end
    ar.ple.showCalculation = do_plotting;  
    return
end
if(ar.qFit(jk)~=1)
    fprintf('\nPLE#%i SKIPPED: parameter %s is fixed\n', jk, ar.ple.p_labels{jk});
    return
else
    fprintf('\nPLE#%i for parameter %s\n', jk, ar.ple.p_labels{jk});
end

% respress warning "Matrix is close to singular or badly scaled."
warn_reset = warning;
warning('off', 'MATLAB:nearlySingularMatrix');

% Setup containers
p = ar.ple.p; %ar.ple.p = ar.p
ar.ple.ps{jk} = nan(2*ar.ple.samplesize(jk)+1, length(p));
ar.ple.psinit{jk} = nan(2*ar.ple.samplesize(jk)+1, length(p));
ar.ple.psinitstep{jk} = nan(2*ar.ple.samplesize(jk)+1, length(p));
ar.ple.dpStep{jk} = nan(2*ar.ple.samplesize(jk)+1, length(p));
ar.ple.dpLast{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.chi2s{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.chi2sviolations{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.chi2spriors{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.chi2spriorsAll{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.chi2sinit{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.gradient{jk} = nan(2*ar.ple.samplesize(jk)+1, length(p));
jindex = ar.ple.samplesize(jk)+1;
% ar.ple.psinit     = parameter after step, before optimization
% ar.ple.chi2sinit  = corresponding merit function value
% ar.ple.psinitstep = parameter step
% ar.ple.dpStep     = actual step after fitting
% ar.ple.dpLast     = step choice algorithm output, test algorithm intern
%                       properties
% ar.ple.gradient   = constrained gradient
% ar.ple.gradient is somehow used in plePlotMulti to check for some
% boundary condition (???)
% Proposal: Transfer all those commented fields into a struct ar.ple.test
%   to accentuate their role as fields which are mostly used for testing.

% Initial Integration:
feval(ar.ple.integrate_fkt, ar.ple.p);
ar.ple.chi2sinit{jk}(jindex) = feval(ar.ple.merit_fkt);
ar.ple.psinit{jk}(jindex,:) = ar.ple.p;
% Fields ending with _fkt contain function handles specified in arPLEInit.
% Here: Integration of the ODEs and updating of corresponding fields.

% Initial fit:
try
    [p, gradient_start] = feval(ar.ple.fit_fkt, jk);
    ar.ple.ps{jk}(jindex,:) = p;
    ar.ple.gradient{jk}(jindex,:) = gradient_start;
catch exception
    fprintf('ERROR PLE: at initial fit (%s)\n', exception.message);
    return
end
% ar.ple.fit_fkt fits all but the profile parameter.

% ps_start if defined after the local fit with the profile parameter already 
% constrained. Use this to initialize left and right branch correctly.
ps_start = p;
ar.ple.psinitstep{jk}(jindex,:) = zeros(size(p));
ar.ple.merit = feval(ar.ple.merit_fkt);
ar.ple.chi2s{jk}(jindex) = ar.ple.merit;
ar.ple.dpLast{jk}(jindex) = 0;
ar.ple.dpStep{jk}(jindex,:) = zeros(size(p)); 

if(isfield(ar.ple,'violations'))
    ar.ple.chi2sviolations{jk}(jindex) = feval(ar.ple.violations);
end
if(isfield(ar.ple,'priors'))
    ar.ple.chi2spriors{jk}(jindex) = feval(ar.ple.priors, jk);
end
if(isfield(ar.ple,'priorsAll'))
    ar.ple.chi2spriorsAll{jk}(jindex) = feval(ar.ple.priorsAll);
end

%% Algorithm

if(~ar.ple.breakon_point)
    dchi2 = ar.ple.dchi2;
else
    dchi2 = ar.ple.dchi2_point;
end

% Initialize process towards upper bound:
pLast = ps_start;
feval(ar.ple.integrate_fkt, pLast);
dpLast = ar.ple.maxstepsize(jk)/10; %Initially proposed step.
% Note that pLast is a vector, but dpLast is a scalar.
dpStep = NaN; %There is no previous step

arWaitbar(0);

try
    if doLeftRightBranch(2)
    for j=1:ar.ple.samplesize(jk)
        jindex = (ar.ple.samplesize(jk)+1) + j;
        arWaitbar(j, ar.ple.samplesize(jk), sprintf('PLE#%i estimating upper confidence bound for %s', ...
            jk, strrep(ar.ple.p_labels{jk},'_', '\_')));
        
        % Choose a parameter step:
        [pStep, dpLast] = feval(ar.ple.initstep_fkt, jk, pLast, dpLast,dpStep); 
        % ar.ple.initstep_fkt is the algorithm for choosing a step.
        % dpLast is used as the starting stepsize for the next step choice
        % algorithm.
        
        if(sum(isnan(pStep))>0)
            break;
        end
        % This exit criterion is realized in the step choice algorithm.
        
        % Make the step and integrate the ODE:
        pTrial = pLast + pStep;
        
        feval(ar.ple.integrate_fkt, pTrial);
        ar.ple.dpLast{jk}(jindex) = dpLast;
        ar.ple.chi2sinit{jk}(jindex) = feval(ar.ple.merit_fkt);
        ar.ple.psinit{jk}(jindex,:) = pTrial;
        ar.ple.psinitstep{jk}(jindex,:) = pStep;
        
        % Fit to new constrained optimum:
        [p, gradient] = feval(ar.ple.fit_fkt, jk);
        
        pLast = p; %Make next step from this parameter     
        ar.ple.ps{jk}(jindex,:) = pLast;
        dpStep = ar.ple.ps{jk}(jindex,:)-ar.ple.ps{jk}(jindex-1,:);
        ar.ple.dpStep{jk}(jindex,:) = dpStep;
        ar.ple.gradient{jk}(jindex,:) = gradient;
        ar.ple.chi2s{jk}(jindex) = feval(ar.ple.merit_fkt);

        if(isfield(ar.ple,'violations'))
            ar.ple.chi2sviolations{jk}(jindex) = feval(ar.ple.violations);
        end
        if(isfield(ar.ple,'priors'))
            ar.ple.chi2spriors{jk}(jindex) = feval(ar.ple.priors, jk);
        end
        if(isfield(ar.ple,'priorsAll'))
            ar.ple.chi2spriorsAll{jk}(jindex) = feval(ar.ple.priorsAll);
        end
        
        if(ar.ple.showCalculation)
            try
                plePlot(jk);
            end
        end
        if(isfield(ar.ple, 'continuousSave') && ar.ple.continuousSave)
            pleSave(ar.ple);
        end
            
        if(feval(ar.ple.merit_fkt) > ar.ple.merit+dchi2*1.2)
            break
        end
        %Determine how far the threshold should be exceeded before
        %termination.
    end
    end
catch exception
    fprintf('ERROR PLE: going to upper bound (%s)\n', exception.message);
end

arWaitbar(-1);

% Initialize process towards lower bound:
% Reset starting parameters.
pLast = ps_start; 
feval(ar.ple.integrate_fkt, pLast);
dpLast = -ar.ple.maxstepsize(jk)/10; 
dpStep = NaN;

arWaitbar(0);
try
    if doLeftRightBranch(1)
    for j=1:ar.ple.samplesize(jk)
        jindex = (ar.ple.samplesize(jk)+1) - j;
        arWaitbar(j, ar.ple.samplesize(jk), sprintf('PLE#%i estimating lower confidence bound for %s', ...
            jk, strrep(ar.ple.p_labels{jk},'_', '\_')));
        
        % estimate intial step
        [pStep, dpLast] = feval(ar.ple.initstep_fkt, jk, pLast, dpLast, dpStep);
        if(sum(isnan(pStep))>0)
            break;
        end
        pTrial = pLast + pStep;
        
        feval(ar.ple.integrate_fkt, pTrial);
        ar.ple.dpLast{jk}(jindex) = dpLast;
        ar.ple.chi2sinit{jk}(jindex) = feval(ar.ple.merit_fkt);
        ar.ple.psinit{jk}(jindex,:) = pTrial;
        ar.ple.psinitstep{jk}(jindex,:) = pStep;
        
        % Fit
        [p, gradient] = feval(ar.ple.fit_fkt, jk); 
        
        pLast = p;        
        ar.ple.ps{jk}(jindex,:) = pLast;
        dpStep = ar.ple.ps{jk}(jindex,:)-ar.ple.ps{jk}(jindex+1,:);
        ar.ple.dpStep{jk}(jindex,:) = dpStep;
        ar.ple.gradient{jk}(jindex,:) = gradient;
        ar.ple.chi2s{jk}(jindex) = feval(ar.ple.merit_fkt);
        if(isfield(ar.ple,'violations'))
            ar.ple.chi2sviolations{jk}(jindex) = feval(ar.ple.violations);
        end
        if(isfield(ar.ple,'priors'))
            ar.ple.chi2spriors{jk}(jindex) = feval(ar.ple.priors, jk);
        end
        if(isfield(ar.ple,'priorsAll'))
            ar.ple.chi2spriorsAll{jk}(jindex) = feval(ar.ple.priorsAll);
        end
        
        if(ar.ple.showCalculation)
            try 
                plePlot(jk);
            end
        end
        if(isfield(ar.ple, 'continuousSave') && ar.ple.continuousSave)
            pleSave(ar.ple);
        end
            
        if(feval(ar.ple.merit_fkt) > ar.ple.merit+dchi2*1.2)
            break
        end
    end
    end
catch exception
    fprintf('ERROR PLE: going to lower bound (%s)\n', exception.message);
end

arWaitbar(-1);

% reset warnings
warning(warn_reset);

% reset parameters
feval(ar.ple.integrate_fkt, ar.ple.p);
% Caution: This resets to the parameters which were in ar.p at the time of 
%            intializing. However, the profile minimum parameters are
%            obtained after an additional fit starting from this point and
%            thus may differ.

%% Finalizing

chi2 = ar.ple.chi2s{jk};
ps = ar.ple.ps{jk};
%These variables now contain parameters and chi2s for the whole profile

q_chi2good = chi2<=min(chi2)+ar.ple.dchi2 & ~isnan(chi2);
q_chi2good_point = chi2<=min(chi2)+ar.ple.dchi2_point & ~isnan(chi2);
% Identifiers which points are under the threshold

% Define new optimum
if(ar.ple.merit-min(chi2) > ar.ple.optimset_tol)
    [minchi2, iminchi2] = min(chi2);
    fprintf('PLE#%i found better optimum with chi^2 decrease of %e\n', jk, ...
        ar.ple.merit - minchi2);
    
    if(ar.ple.allowbetteroptimum)
        ar.ple.merit = minchi2;
        ar.ple.p = ps(iminchi2,:);
        feval(ar.ple.setoptim_fkt, ps(iminchi2,:));
    end
end

% IDstatus is used to define what kind of identifiability is given.
% Check if there is a calculation error
if(sum(~isnan(ar.ple.chi2s{jk})) < 3)
    ar.ple.IDstatus(jk) = 4;
    ar.ple.IDstatus_point(jk) = 4;
else

    % Calculate CI simultaneous if requested.
    
    if(~ar.ple.breakon_point)
        if(min(ps(q_chi2good,jk))==min(ps(~isnan(chi2),jk)))
            ar.ple.conf_lb(jk) = -Inf;
            % bounds are infinite if chi2 does not exceed boundary 
        else
            kind = find(ps(:,jk)==min(ps(q_chi2good,jk)));
            try
                ar.ple.conf_lb(jk) = interp1(chi2([kind kind-1]), ps([kind kind-1], jk), min(chi2)+ar.ple.dchi2);
            % interpolate to increase accuracy of confidence interval
            catch
                ar.ple.conf_lb(jk) = NaN;  % e.g. isinf(chi2(kind-1))
            end
        end
        if(max(ps(q_chi2good,jk))==max(ps(~isnan(chi2),jk)))
            ar.ple.conf_ub(jk) = Inf;
        else
            kind = find(ps(:,jk)==max(ps(q_chi2good,jk)));
            try
                ar.ple.conf_ub(jk) = interp1(chi2([kind kind+1]), ps([kind kind+1], jk), min(chi2)+ar.ple.dchi2);
            catch
                ar.ple.conf_ub(jk) = NaN;
            end

        end
    end

    % calculate CI point-wise in every case

    if(min(ps(q_chi2good_point,jk))==min(ps(~isnan(chi2),jk)))
        ar.ple.conf_lb_point(jk) = -Inf; 
        % bounds are infinite if chi2 does not exceed boundary 
    else
        kind = find(ps(:,jk)==min(ps(q_chi2good_point,jk)));
        try
            ar.ple.conf_lb_point(jk) = interp1(chi2([kind kind-1]), ps([kind kind-1], jk), min(chi2)+ar.ple.dchi2_point);
            % interpolate to increase accuracy of confidence interval
        catch
            ar.ple.conf_lb_point(jk) = NaN;
        end
    end
    if(max(ps(q_chi2good_point,jk))==max(ps(~isnan(chi2),jk)))
        ar.ple.conf_ub_point(jk) = Inf;
    else
        kind = find(ps(:,jk)==max(ps(q_chi2good_point,jk)));
        try
            ar.ple.conf_ub_point(jk) = interp1(chi2([kind kind+1]), ps([kind kind+1], jk), min(chi2)+ar.ple.dchi2_point);
        catch
            ar.ple.conf_ub_point(jk) = NaN;
        end
    end

    % check ID point-wise
    if(max(chi2(~isnan(chi2)))<(ar.ple.dchi2_point*ar.ple.chi2_strID_ratio)+min(chi2(~isnan(chi2))))
        ar.ple.IDstatus_point(jk) = 3; % structurally NI
    else
        if((ar.ple.conf_lb_point(jk)==-Inf || ar.ple.conf_ub_point(jk)==Inf))% && ...
                %(isnan(ar.ple.IDstatus_point(jk)) || (ar.ple.IDstatus_point(jk)<3)))
            ar.ple.IDstatus_point(jk) = 2; % practically NI
        else
            ar.ple.IDstatus_point(jk) = 1; % identifiable
        end
    end

    % check ID simultaneous
    if(~ar.ple.breakon_point)
        % structural
        if(max(chi2(~isnan(chi2)))<(ar.ple.dchi2*ar.ple.chi2_strID_ratio)+min(chi2(~isnan(chi2))))
            ar.ple.IDstatus(jk) = 3;
        else
            % practical
            if((ar.ple.conf_lb(jk)==-Inf || ar.ple.conf_ub(jk)==Inf))% && ...
                %(isnan(ar.ple.IDstatus(jk)) || (ar.ple.IDstatus(jk)<3)))
                ar.ple.IDstatus(jk) = 2;
            else
                ar.ple.IDstatus(jk) = 1;
            end
        end
    end

    % calulate coverage
    if(isfield(ar.ple, 'p_true'))
        ar.ple.cover = ar.ple.p_true <= ar.ple.conf_ub & ...
            ar.ple.p_true >= ar.ple.conf_lb;
        ar.ple.cover_point = ar.ple.p_true <= ar.ple.conf_ub_point & ...
            ar.ple.p_true >= ar.ple.conf_lb_point;
    end
end

%Since ple is called repeatedly, this is executed for every parameter.
ar.ple.finished = 1;   
if (ar.ple.autosave)
    pleSave(ar)
end

