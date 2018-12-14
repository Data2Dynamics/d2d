% ple([i], [samplesize], [relchi2stepincrease], [maxstepsize], [minstepsize], [breakonbounds])
% 
% Profile Likelihood Exploit
%
%   i                     parameter index or vector of par indices [all if not specified]
%                         alternatively the name of a parameter or a cell of
%                         names can be provided. If one parameter name is
%                         provided which is not in ar.ple.p_labels, then the string
%                         is interpreted as a regular expression.
%   samplesize            number of sampling steps            [100]
%   relchi2stepincrease   percentage chi^2 increase of a step [0.1]
%   maxstepsize           maximum size of a step              [0.2 * p_jk]
%   minstepsize           minumum size of a step              [1e-6]
%   breakonlb             stop if hit lb                      [false]
%   breakonub             stop if hit ub                      [true]
% 
% The profile likelihood calculation by the functions ple* was intended
% as running independent of D2D, i.e. it was intended to be also used by
% other tools. Therefore, the function does not use info in global ar and stores the
% results an a separate variable.
%
% Before usage, workspace has to be initialized.
%
% See also: arPLEInit, pleExtend

function ple(jk, samplesize, relchi2stepincrease, ...
    maxstepsize, minstepsize, breakonlb, breakonub)

global ar

if(~isfield(ar,'ple') || isempty(ar.ple))
    error('PLE ERROR: please initialize')
end 
if(~isfield(ar.ple, 'showCalculation'))
    ar.ple.showCalculation = true;
end

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
    ar.ple.tmean = mean(ar.ple.timing(ar.qFit==1));
    ar.ple.tstd = std(ar.ple.timing(ar.qFit==1));    
    fprintf('\nPLE mean elapsed time for %i parameter: %s +/- %s\n', ...
        sum(ar.qFit==1), secToHMS(ar.ple.tmean), secToHMS(ar.ple.tstd));
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
    ar.ple.tmean = mean(ar.ple.timing(ar.qFit(jk)));
    ar.ple.tstd = std(ar.ple.timing(ar.qFit(jk)));    
    fprintf('\nPLE mean elapsed time for %i parameter: %s +/- %s\n', ...
        length(jk), secToHMS(ar.ple.tmean), secToHMS(ar.ple.tstd));
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

ar.ple.conf_labels = {'Hessian', 'PLE'};

% setup containers
p = ar.ple.p;
ar.ple.ps{jk} = nan(2*ar.ple.samplesize(jk)+1, length(p));
ar.ple.psinit{jk} = nan(2*ar.ple.samplesize(jk)+1, length(p));
ar.ple.psinitstep{jk} = nan(2*ar.ple.samplesize(jk)+1, length(p));
ar.ple.chi2s{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.chi2sviolations{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.chi2spriors{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.chi2spriorsAll{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.chi2sinit{jk} = nan(1,2*ar.ple.samplesize(jk)+1);
ar.ple.gradient{jk} = nan(2*ar.ple.samplesize(jk)+1, length(p));
jindex = ar.ple.samplesize(jk)+1;

% initial fit
feval(ar.ple.integrate_fkt, ar.ple.p);
ar.ple.chi2sinit{jk}(jindex) = feval(ar.ple.merit_fkt);
ar.ple.psinit{jk}(jindex,:) = ar.ple.p;
try
    [p, gradient_start] = feval(ar.ple.fit_fkt, jk);
    ar.ple.ps{jk}(jindex,:) = p;
    ar.ple.gradient{jk}(jindex,:) = gradient_start;
catch exception
    fprintf('ERROR PLE: at initial fit (%s)\n', exception.message);
    return
end
ps_start = p;
ar.ple.psinitstep{jk}(jindex,:) = zeros(size(p));
ar.ple.merit = feval(ar.ple.merit_fkt);
last.y = ar.ple.chi2sinit{jk}(jindex)-ar.ple.merit;
ar.ple.chi2s{jk}(jindex) = ar.ple.merit;
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

pLast = ps_start;
feval(ar.ple.integrate_fkt, pLast);
dpLast = ar.ple.maxstepsize(jk)/2.1; % exploiting upper bound, base step
last.dx = NaN(1,ar.ple.samplesize(jk));
last.dy = NaN(1,ar.ple.samplesize(jk));
last.dy(1) = dpLast;

estimatetime = 0;
fittime = 0;
arWaitbar(0);

try
    for j=1:ar.ple.samplesize(jk)
        jindex = (ar.ple.samplesize(jk)+1) + j;
        arWaitbar(j, ar.ple.samplesize(jk), sprintf('PLE#%i estimating upper confidence bound for %s', ...
            jk, strrep(ar.ple.p_labels{jk},'_', '\_')));
        
        tic;
        % estimate (intial) step
        last.x = pLast(jk);
        [pStep, dpLast] = feval(ar.ple.initstep_fkt, jk, pLast, dpLast, last);
        last.dx(j) = dpLast;   
        if(sum(isnan(pStep))>0)
            break;
        end
        pTrial = pLast + pStep;
        
        feval(ar.ple.integrate_fkt, pTrial);
        ar.ple.chi2sinit{jk}(jindex) = feval(ar.ple.merit_fkt);
        last.y = ar.ple.chi2sinit{jk}(jindex)-ar.ple.merit;
        ar.ple.psinit{jk}(jindex,:) = pTrial;
        ar.ple.psinitstep{jk}(jindex,:) = pStep;
        
        estimatetime = estimatetime + toc;
        
        tic;
        % Fit
        [p, gradient] = feval(ar.ple.fit_fkt, jk);
        if(length(dpLast)>1)
            dpLast = p - pLast;
        end
        pLast = p;
        
        ar.ple.ps{jk}(jindex,:) = pLast;
        ar.ple.gradient{jk}(jindex,:) = gradient;
        ar.ple.chi2s{jk}(jindex) = feval(ar.ple.merit_fkt);
        last.y = ar.ple.chi2s{jk}(jindex)-ar.ple.merit;
        
        last.dy(j) = ar.ple.chi2s{jk}(jindex) - ar.ple.chi2s{jk}(jindex-1);
        if(isfield(ar.ple,'violations'))
            ar.ple.chi2sviolations{jk}(jindex) = feval(ar.ple.violations);
        end
        if(isfield(ar.ple,'priors'))
            ar.ple.chi2spriors{jk}(jindex) = feval(ar.ple.priors, jk);
        end
        if(isfield(ar.ple,'priorsAll'))
            ar.ple.chi2spriorsAll{jk}(jindex) = feval(ar.ple.priorsAll);
        end
        fittime = fittime + toc;
        
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
catch exception
    fprintf('ERROR PLE: going to upper bound (%s)\n', exception.message);
end

arWaitbar(-1);

pLast = ps_start;
feval(ar.ple.integrate_fkt, pLast);
dpLast = -ar.ple.maxstepsize(jk)/2.1; % exploiting lower bound, base step
last.dx = NaN(1,ar.ple.samplesize(jk));
last.dy = NaN(1,ar.ple.samplesize(jk));

arWaitbar(0);
try
    for j=1:ar.ple.samplesize(jk)
        jindex = (ar.ple.samplesize(jk)+1) - j;
        arWaitbar(j, ar.ple.samplesize(jk), sprintf('PLE#%i estimating lower confidence bound for %s', ...
            jk, strrep(ar.ple.p_labels{jk},'_', '\_')));
        
        tic;
        % estimate intial step
        last.x = pLast(jk);
        [pStep, dpLast] = feval(ar.ple.initstep_fkt, jk, pLast, dpLast, last);
        last.dx(j) = dpLast;
        if(sum(isnan(pStep))>0)
            break;
        end
        pTrial = pLast + pStep;
        
        feval(ar.ple.integrate_fkt, pTrial);
        ar.ple.chi2sinit{jk}(jindex) = feval(ar.ple.merit_fkt);
        last.y = ar.ple.chi2sinit{jk}(jindex)-ar.ple.merit;
        ar.ple.psinit{jk}(jindex,:) = pTrial;
        ar.ple.psinitstep{jk}(jindex,:) = pStep;
        
        estimatetime = estimatetime + toc;
        
        tic;
        % Fit
        [p, gradient] = feval(ar.ple.fit_fkt, jk); 
        if(length(dpLast)>1)
            dpLast = p - pLast;
        end
        pLast = p;
        
        ar.ple.ps{jk}(jindex,:) = pLast;
        ar.ple.gradient{jk}(jindex,:) = gradient;
        ar.ple.chi2s{jk}(jindex) = feval(ar.ple.merit_fkt);
        last.y = ar.ple.chi2s{jk}(jindex)-ar.ple.merit;
        last.dy(j) = ar.ple.chi2s{jk}(jindex) - ar.ple.chi2s{jk}(jindex+1);
        if(isfield(ar.ple,'violations'))
            ar.ple.chi2sviolations{jk}(jindex) = feval(ar.ple.violations);
        end
        if(isfield(ar.ple,'priors'))
            ar.ple.chi2spriors{jk}(jindex) = feval(ar.ple.priors, jk);
        end
        if(isfield(ar.ple,'priorsAll'))
            ar.ple.chi2spriorsAll{jk}(jindex) = feval(ar.ple.priorsAll);
        end
        fittime = fittime + toc;
        
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
catch exception
    fprintf('ERROR PLE: going to lower bound (%s)\n', exception.message);
end

arWaitbar(-1);

% reset warnings
warning(warn_reset);

% reset parameters
feval(ar.ple.integrate_fkt, ar.ple.p);

%% Finalizing

chi2 = ar.ple.chi2s{jk};
ps = ar.ple.ps{jk};

ar.ple.estimatetime(jk) = estimatetime;
ar.ple.fittime(jk) = fittime;
ar.ple.timing(jk) = fittime+estimatetime;

q_chi2good = chi2<=min(chi2)+ar.ple.dchi2 & ~isnan(chi2);
q_chi2good_point = chi2<=min(chi2)+ar.ple.dchi2_point & ~isnan(chi2);

% define new optimum
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

% check calculation error
if(sum(~isnan(ar.ple.chi2s{jk})) < 3)
    ar.ple.IDstatus(jk) = 4;
    ar.ple.IDstatus_point(jk) = 4;
else

    % calculate CI simultaneous
    if(~ar.ple.breakon_point)
        ar.ple.conf_lb(jk) = min(ps(q_chi2good,jk));
        ar.ple.conf_ub(jk) = max(ps(q_chi2good,jk));
        
        if(min(ps(q_chi2good,jk))==min(ps(~isnan(chi2),jk)))
            ar.ple.conf_lb(jk) = -Inf;
        else
            kind = find(ps(:,jk)==min(ps(q_chi2good,jk)));
            ar.ple.conf_lb(jk) = interp1(chi2([kind kind-1]), ps([kind kind-1], jk), min(chi2)+ar.ple.dchi2);
        end
        if(max(ps(q_chi2good,jk))==max(ps(~isnan(chi2),jk)))
            ar.ple.conf_ub(jk) = Inf;
        else
            kind = find(ps(:,jk)==max(ps(q_chi2good,jk)));
            ar.ple.conf_ub(jk) = interp1(chi2([kind kind+1]), ps([kind kind+1], jk), min(chi2)+ar.ple.dchi2);
        end
    end

    % calculate CI point-wise
    ar.ple.conf_lb_point(jk) = min(ps(q_chi2good_point,jk));
    ar.ple.conf_ub_point(jk) = max(ps(q_chi2good_point,jk));

    if(min(ps(q_chi2good_point,jk))==min(ps(~isnan(chi2),jk)))
        ar.ple.conf_lb_point(jk) = -Inf;
    else
        kind = find(ps(:,jk)==min(ps(q_chi2good_point,jk)));
        try
            ar.ple.conf_lb_point(jk) = interp1(chi2([kind kind-1]), ps([kind kind-1], jk), min(chi2)+ar.ple.dchi2_point);
        catch
            ar.ple.conf_lb_point(jk) = NaN;  % e.g. isinf(chi2(kind-1))
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
    % structural
    if(max(chi2(~isnan(chi2)))<(ar.ple.dchi2_point*ar.ple.chi2_strID_ratio)+min(chi2(~isnan(chi2))))
        ar.ple.IDstatus_point(jk) = 3;
    else
        % practical
        if((ar.ple.conf_lb_point(jk)==-Inf || ar.ple.conf_ub_point(jk)==Inf))% && ...
                %(isnan(ar.ple.IDstatus_point(jk)) || (ar.ple.IDstatus_point(jk)<3)))
            ar.ple.IDstatus_point(jk) = 2;
        else
            ar.ple.IDstatus_point(jk) = 1;
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

    % calulate relative CIs
    ar.ple.conf_rel = abs((ar.ple.conf_ub(:) - ar.ple.conf_lb(:))/2*100 ./ ar.ple.p(:));
    ar.ple.conf_rel_point = abs((ar.ple.conf_ub_point(:) - ar.ple.conf_lb_point(:))/2*100 ./ ar.ple.p(:));

    % calulate coverage
    if(isfield(ar.ple, 'p_true'))
        ar.ple.cover = ar.ple.p_true <= ar.ple.conf_ub & ...
            ar.ple.p_true >= ar.ple.conf_lb;
        ar.ple.cover_point = ar.ple.p_true <= ar.ple.conf_ub_point & ...
            ar.ple.p_true >= ar.ple.conf_lb_point;
    end
    
    % calculate elapsed times
    rel_estimate = 100 * estimatetime / (estimatetime+fittime);
    rel_fit = 100 * fittime / (estimatetime+fittime);
    fprintf('PLE#%i elapsed time %s (step: %2.0f%%, fit: %2.0f%%)\n', jk, ...
        secToHMS(ar.ple.timing(jk)), rel_estimate, rel_fit);
end


%% Output
% newIDlables = {'identifiable','practically non-identifiable','structurally non-identifiable','calculation error'};
% if(~ar.ple.breakon_point)
%     fprintf('PLE#%i suggesting: %s (simultaneous)\n', jk, newIDlables{ar.ple.IDstatus(jk)});
% else
%     fprintf('PLE#%i suggesting: %s (point-wise)\n', jk, newIDlables{ar.ple.IDstatus_point(jk)});
% end

ar.ple.finished = 1;    
pleSave(ar)

