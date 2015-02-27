% Profile Likelihood Exploit
%
% ple([i, samplesize, relchi2stepincrease, maxstepsize, minstepsize, breakonbounds])
%
% i:                    i'th parameter, see pwInfo
%                       [if omitted, all free parameters are considered]
%                       alternatively the name of a parameter or a cell of
%                       names can be provided. If one parameter name is
%                       provided which is not in ar.pLabel, then the string
%                       is interpreted as a regular expression.
% samplesize:           number of sampling steps            [100]
% relchi2stepincrease:  percentage chi^2 increase of a step [0.1]
% maxstepsize:          maximum size of a step              [0.2 * p_jk]
% minstepsize:          minumum size of a step              [1e-6]
% breakonlb:            stop if hit lb                      [false]
% breakonub:            stop if hit ub                      [true]

function ple(jk, samplesize, relchi2stepincrease, ...
    maxstepsize, minstepsize, breakonlb, breakonub)

global pleGlobals;

if(isempty(pleGlobals))
    error('PLE ERROR: please initialize')
end 
if(~isfield(pleGlobals, 'showCalculation'))
    pleGlobals.showCalculation = true;
end

if(ischar(jk))
    global ar
    tref = strmatch(jk,ar.pLabel,'exact');
    if(isempty(tref))
        jk = find(~cellfun(@isempty,regexp(ar.pLabel,jk)));
    else 
        jk = tref;
    end
    
    if isempty(jk)
        disp('Pattern ''',jk,''' not found in ar.pLabel');
        return;
    end
elseif(iscell(jk)) % cell of pLabels
    [~,jk] = intersect(ar.pLabel,jk);
elseif(isnumeric(jk))
else
    error('Argument has to be a string or an array of indices.')
end

if(nargin<1)
    fprintf('PLE for %i parameters ...\n', sum(pleGlobals.q_fit))
    jindex = find(pleGlobals.q_fit);
    do_plotting = pleGlobals.showCalculation;
    pleGlobals.showCalculation = false;
    for j=1:length(jindex)
        ple(jindex(j))
        if(do_plotting)
            plePlotMulti;
        end
    end
    pleGlobals.showCalculation = do_plotting;
    pleGlobals.tmean = mean(pleGlobals.timing(pleGlobals.q_fit));
    pleGlobals.tstd = std(pleGlobals.timing(pleGlobals.q_fit));    
    fprintf('\nPLE mean elapsed time for %i parameter: %s +/- %s\n', ...
        sum(pleGlobals.q_fit), secToHMS(pleGlobals.tmean), secToHMS(pleGlobals.tstd));
    return
elseif(length(jk)>1)
    fprintf('PLE for %i parameters ...\n', length(jk))
    do_plotting = pleGlobals.showCalculation;
    pleGlobals.showCalculation = false;
    for j=1:length(jk)
        ple(jk(j))
        if(do_plotting)
            plePlotMulti;
        end
    end
    pleGlobals.showCalculation = do_plotting;
    pleGlobals.tmean = mean(pleGlobals.timing(pleGlobals.q_fit(jk)));
    pleGlobals.tstd = std(pleGlobals.timing(pleGlobals.q_fit(jk)));    
    fprintf('\nPLE mean elapsed time for %i parameter: %s +/- %s\n', ...
        length(jk), secToHMS(pleGlobals.tmean), secToHMS(pleGlobals.tstd));
    return
end
if(~pleGlobals.q_fit(jk))
    fprintf('\nPLE#%i SKIPPED: parameter %s is fixed\n', jk, pleGlobals.p_labels{jk});
    return
else
    fprintf('\nPLE#%i for parameter %s\n', jk, pleGlobals.p_labels{jk});
end

if(exist('samplesize', 'var'))
    pleGlobals.samplesize(jk) = samplesize;
end
if(exist('relchi2stepincrease', 'var'))
    pleGlobals.relchi2stepincrease(jk) = relchi2stepincrease;
end
if(exist('maxstepsize', 'var'))
    pleGlobals.maxstepsize(jk) = maxstepsize;
end
if(exist('minstepsize', 'var'))
    pleGlobals.minstepsize(jk) = minstepsize;
end
if(exist('breakonlb', 'var'))
    pleGlobals.breakonlb(jk) = breakonlb;
end
if(exist('breakonub', 'var'))
    pleGlobals.breakonub(jk) = breakonub;
end

% respress warning "Matrix is close to singular or badly scaled."
warn_reset = warning;
warning('off', 'MATLAB:nearlySingularMatrix');

pleGlobals.conf_labels = {'Hessian', 'PLE'};

% setup containers
p = pleGlobals.p;
pleGlobals.ps{jk} = nan(2*pleGlobals.samplesize(jk)+1, length(p));
pleGlobals.psinit{jk} = nan(2*pleGlobals.samplesize(jk)+1, length(p));
pleGlobals.psinitstep{jk} = nan(2*pleGlobals.samplesize(jk)+1, length(p));
pleGlobals.chi2s{jk} = nan(1,2*pleGlobals.samplesize(jk)+1);
pleGlobals.chi2sviolations{jk} = nan(1,2*pleGlobals.samplesize(jk)+1);
pleGlobals.chi2spriors{jk} = nan(1,2*pleGlobals.samplesize(jk)+1);
pleGlobals.chi2sinit{jk} = nan(1,2*pleGlobals.samplesize(jk)+1);
pleGlobals.gradient{jk} = nan(2*pleGlobals.samplesize(jk)+1, length(p));
jindex = pleGlobals.samplesize(jk)+1;

% initial fit
feval(pleGlobals.integrate_fkt, pleGlobals.p);
pleGlobals.chi2sinit{jk}(jindex) = feval(pleGlobals.merit_fkt);
pleGlobals.psinit{jk}(jindex,:) = pleGlobals.p;
try
    [p, gradient_start] = feval(pleGlobals.fit_fkt, jk);
    pleGlobals.ps{jk}(jindex,:) = p;
    pleGlobals.gradient{jk}(jindex,:) = gradient_start;
catch exception
    fprintf('ERROR PLE: at initial fit (%s)\n', exception.message);
    return
end
ps_start = p;
pleGlobals.psinitstep{jk}(jindex,:) = zeros(size(p));
pleGlobals.chi2 = feval(pleGlobals.merit_fkt);
pleGlobals.chi2s{jk}(jindex) = pleGlobals.chi2;
if(isfield(pleGlobals,'violations'))
    pleGlobals.chi2sviolations{jk}(jindex) = feval(pleGlobals.violations);
end
if(isfield(pleGlobals,'priors'))
    pleGlobals.chi2spriors{jk}(jindex) = feval(pleGlobals.priors, jk);
end

%% Algorithm

if(~pleGlobals.breakon_point)
    dchi2 = pleGlobals.dchi2;
else
    dchi2 = pleGlobals.dchi2_point;
end

pLast = ps_start;
feval(pleGlobals.integrate_fkt, pLast);
dpLast = pleGlobals.maxstepsize(jk)/2.1; % exploiting upper bound, base step

estimatetime = 0;
fittime = 0;
arWaitbar(0);

try
    for j=1:pleGlobals.samplesize(jk)
        jindex = (pleGlobals.samplesize(jk)+1) + j;
        arWaitbar(j, pleGlobals.samplesize(jk), sprintf('PLE#%i estimating upper confidence bound for %s', ...
            jk, strrep(pleGlobals.p_labels{jk},'_', '\_')));
        
        tic;
        % estimate intial step
        [pStep, dpLast] = feval(pleGlobals.initstep_fkt, jk, pLast, dpLast);
        if(sum(isnan(pStep))>0)
            break;
        end
        pTrial = pLast + pStep;
        
        feval(pleGlobals.integrate_fkt, pTrial);
        pleGlobals.chi2sinit{jk}(jindex) = feval(pleGlobals.merit_fkt);
        pleGlobals.psinit{jk}(jindex,:) = pTrial;
        pleGlobals.psinitstep{jk}(jindex,:) = pStep;
        
        estimatetime = estimatetime + toc;
        
        tic;
        % Fit
        [p, gradient] = feval(pleGlobals.fit_fkt, jk);
        if(length(dpLast)>1)
            dpLast = p - pLast;
        end
        pLast = p;
        
        pleGlobals.ps{jk}(jindex,:) = pLast;
        pleGlobals.gradient{jk}(jindex,:) = gradient;
        pleGlobals.chi2s{jk}(jindex) = feval(pleGlobals.merit_fkt);
        if(isfield(pleGlobals,'violations'))
            pleGlobals.chi2sviolations{jk}(jindex) = feval(pleGlobals.violations);
        end
        if(isfield(pleGlobals,'priors'))
            pleGlobals.chi2spriors{jk}(jindex) = feval(pleGlobals.priors, jk);
        end
        fittime = fittime + toc;
        
        if(pleGlobals.showCalculation)
            try
                plePlot(jk);
            end
        end
            
        if(feval(pleGlobals.merit_fkt) > pleGlobals.chi2+dchi2*1.2)
            break
        end
    end
catch exception
    fprintf('ERROR PLE: going to upper bound (%s)\n', exception.message);
end

arWaitbar(-1);

pLast = ps_start;
feval(pleGlobals.integrate_fkt, pLast);
dpLast = -pleGlobals.maxstepsize(jk)/2.1; % exploiting lower bound, base step

arWaitbar(0);
try
    for j=1:pleGlobals.samplesize(jk)
        jindex = (pleGlobals.samplesize(jk)+1) - j;
        arWaitbar(j, pleGlobals.samplesize(jk), sprintf('PLE#%i estimating lower confidence bound for %s', ...
            jk, strrep(pleGlobals.p_labels{jk},'_', '\_')));
        
        tic;
        % estimate intial step
        [pStep, dpLast] = feval(pleGlobals.initstep_fkt, jk, pLast, dpLast);
        if(sum(isnan(pStep))>0)
            break;
        end
        pTrial = pLast + pStep;
        
        feval(pleGlobals.integrate_fkt, pTrial);
        pleGlobals.chi2sinit{jk}(jindex) = feval(pleGlobals.merit_fkt);
        pleGlobals.psinit{jk}(jindex,:) = pTrial;
        pleGlobals.psinitstep{jk}(jindex,:) = pStep;
        
        estimatetime = estimatetime + toc;
        
        tic;
        % Fit
        [p, gradient] = feval(pleGlobals.fit_fkt, jk); 
        if(length(dpLast)>1)
            dpLast = p - pLast;
        end
        pLast = p;
        
        pleGlobals.ps{jk}(jindex,:) = pLast;
        pleGlobals.gradient{jk}(jindex,:) = gradient;
        pleGlobals.chi2s{jk}(jindex) = feval(pleGlobals.merit_fkt);
        if(isfield(pleGlobals,'violations'))
            pleGlobals.chi2sviolations{jk}(jindex) = feval(pleGlobals.violations);
        end
        if(isfield(pleGlobals,'priors'))
            pleGlobals.chi2spriors{jk}(jindex) = feval(pleGlobals.priors, jk);
        end
        fittime = fittime + toc;
        
        if(pleGlobals.showCalculation)
            try 
                plePlot(jk);
            end
        end
            
        if(feval(pleGlobals.merit_fkt) > pleGlobals.chi2+dchi2*1.2)
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
feval(pleGlobals.integrate_fkt, pleGlobals.p);

%% Finalizing

chi2 = pleGlobals.chi2s{jk};
ps = pleGlobals.ps{jk};

pleGlobals.estimatetime(jk) = estimatetime;
pleGlobals.fittime(jk) = fittime;
pleGlobals.timing(jk) = fittime+estimatetime;

q_chi2good = chi2<=min(chi2)+pleGlobals.dchi2 & ~isnan(chi2);
q_chi2good_point = chi2<=min(chi2)+pleGlobals.dchi2_point & ~isnan(chi2);

% define new optimum
if(pleGlobals.chi2-min(chi2) > pleGlobals.optimset_tol)
    [minchi2, iminchi2] = min(chi2);
    fprintf('PLE#%i found better optimum with chi^2 decrease of %e\n', jk, ...
        pleGlobals.chi2 - minchi2);
    
    if(pleGlobals.allowbetteroptimum)
        pleGlobals.chi2 = minchi2;
        pleGlobals.p = ps(iminchi2,:);
        feval(pleGlobals.setoptim_fkt, ps(iminchi2,:));
    end
end

% check calculation error
if(sum(~isnan(pleGlobals.chi2s{jk})) < 3)
    pleGlobals.IDstatus(jk) = 4;
    pleGlobals.IDstatus_point(jk) = 4;
else

    % calculate CI simultaneous
    if(~pleGlobals.breakon_point)
        pleGlobals.conf_lb(jk) = min(ps(q_chi2good,jk));
        pleGlobals.conf_ub(jk) = max(ps(q_chi2good,jk));
        
        if(min(ps(q_chi2good,jk))==min(ps(~isnan(chi2),jk)))
            pleGlobals.conf_lb(jk) = -Inf;
        else
            kind = find(ps(:,jk)==min(ps(q_chi2good,jk)));
            pleGlobals.conf_lb(jk) = interp1(chi2([kind kind-1]), ps([kind kind-1], jk), min(chi2)+pleGlobals.dchi2);
        end
        if(max(ps(q_chi2good,jk))==max(ps(~isnan(chi2),jk)))
            pleGlobals.conf_ub(jk) = Inf;
        else
            kind = find(ps(:,jk)==max(ps(q_chi2good,jk)));
            pleGlobals.conf_ub(jk) = interp1(chi2([kind kind+1]), ps([kind kind+1], jk), min(chi2)+pleGlobals.dchi2);
        end
    end

    % calculate CI point-wise
    pleGlobals.conf_lb_point(jk) = min(ps(q_chi2good_point,jk));
    pleGlobals.conf_ub_point(jk) = max(ps(q_chi2good_point,jk));

    if(min(ps(q_chi2good_point,jk))==min(ps(~isnan(chi2),jk)))
        pleGlobals.conf_lb_point(jk) = -Inf;
    else
        kind = find(ps(:,jk)==min(ps(q_chi2good_point,jk)));
        pleGlobals.conf_lb_point(jk) = interp1(chi2([kind kind-1]), ps([kind kind-1], jk), min(chi2)+pleGlobals.dchi2_point);
    end
    if(max(ps(q_chi2good_point,jk))==max(ps(~isnan(chi2),jk)))
        pleGlobals.conf_ub_point(jk) = Inf;
    else
        kind = find(ps(:,jk)==max(ps(q_chi2good_point,jk)));
        pleGlobals.conf_ub_point(jk) = interp1(chi2([kind kind+1]), ps([kind kind+1], jk), min(chi2)+pleGlobals.dchi2_point);
    end

    % check ID point-wise
    % structural
    if(max(chi2(~isnan(chi2)))<(pleGlobals.dchi2_point*pleGlobals.chi2_strID_ratio)+min(chi2(~isnan(chi2))))
        pleGlobals.IDstatus_point(jk) = 3;
    else
        % practical
        if((pleGlobals.conf_lb_point(jk)==-Inf || pleGlobals.conf_ub_point(jk)==Inf))% && ...
                %(isnan(pleGlobals.IDstatus_point(jk)) || (pleGlobals.IDstatus_point(jk)<3)))
            pleGlobals.IDstatus_point(jk) = 2;
        else
            pleGlobals.IDstatus_point(jk) = 1;
        end
    end

    % check ID simultaneous
    if(~pleGlobals.breakon_point)
        % structural
        if(max(chi2(~isnan(chi2)))<(pleGlobals.dchi2*pleGlobals.chi2_strID_ratio)+min(chi2(~isnan(chi2))))
            pleGlobals.IDstatus(jk) = 3;
        else
            % practical
            if((pleGlobals.conf_lb(jk)==-Inf || pleGlobals.conf_ub(jk)==Inf))% && ...
                %(isnan(pleGlobals.IDstatus(jk)) || (pleGlobals.IDstatus(jk)<3)))
                pleGlobals.IDstatus(jk) = 2;
            else
                pleGlobals.IDstatus(jk) = 1;
            end
        end
    end

    % calulate relative CIs
    pleGlobals.conf_rel = abs((pleGlobals.conf_ub - pleGlobals.conf_lb)/2*100 ./ pleGlobals.p);
    pleGlobals.conf_rel_point = abs((pleGlobals.conf_ub_point - pleGlobals.conf_lb_point)/2*100 ./ pleGlobals.p);

    % calulate coverage
    if(isfield(pleGlobals, 'p_true'))
        pleGlobals.cover = pleGlobals.p_true <= pleGlobals.conf_ub & ...
            pleGlobals.p_true >= pleGlobals.conf_lb;
        pleGlobals.cover_point = pleGlobals.p_true <= pleGlobals.conf_ub_point & ...
            pleGlobals.p_true >= pleGlobals.conf_lb_point;
    end
    
    % calculate elapse times
    rel_estimate = 100 * estimatetime / (estimatetime+fittime);
    rel_fit = 100 * fittime / (estimatetime+fittime);
    fprintf('PLE#%i elapsed time %s (step: %2.0f%%, fit: %2.0f%%)\n', jk, ...
        secToHMS(pleGlobals.timing(jk)), rel_estimate, rel_fit);
end


%% Output
% newIDlables = {'identifiable','practically non-identifiable','structurally non-identifiable','calculation error'};
% if(~pleGlobals.breakon_point)
%     fprintf('PLE#%i suggesting: %s (simultaneous)\n', jk, newIDlables{pleGlobals.IDstatus(jk)});
% else
%     fprintf('PLE#%i suggesting: %s (point-wise)\n', jk, newIDlables{pleGlobals.IDstatus_point(jk)});
% end

if(~exist([cd '/' pleGlobals.savePath], 'dir'))
    mkdir([cd '/' pleGlobals.savePath])
end
save([pleGlobals.savePath '/results.mat'], 'pleGlobals');

