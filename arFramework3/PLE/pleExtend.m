% pleExtend([i], [samplesize])
% 
% Profile Likelihood Exploit. 
% Extend existing profile in one direction of parameter axis.
%
%   i                     parameter index or vector of par indices [all if not specified]
%   samplesize            number of sampling steps [100]
%                         (use negative numbers to extend profiles to the left)
%
% See also: arPLEInit, ple

function pleExtend(jk, samplesize)

global ar

if(isempty(ar.ple))
    error('PLE ERROR: please initialize')
end 
if(~isfield(ar.ple, 'showCalculation'))
    ar.ple.showCalculation = true;
end

if(isempty(ar.ple.chi2s))
    error('pleExtend.m: ar.ple does not contain the result of ple. Perform ple first.')
end

if(nargin == 0 || isempty(jk))
    if(~exist('samplesize','var'))
        samplesize = 100;
    end
    fprintf('PLE extending for %i parameters ...\n', sum(ar.qFit))
    jindex = find(ar.qFit(1:length(ar.ple.chi2s)) & ~cellfun(@isempty,ar.ple.chi2s));
    do_plotting = ar.ple.showCalculation;
    ar.ple.showCalculation = false;
    for j=1:length(jindex)
        pleExtend(jindex(j), abs(samplesize))
        if(do_plotting)
            plePlotMulti;
        end
        pleExtend(jindex(j), -abs(samplesize))
        if(do_plotting)
            plePlotMulti;
        end
    end
    ar.ple.showCalculation = do_plotting;
    return
end

if(nargin == 1)
    pleExtend(jk, +100);
    pleExtend(jk, -100);
    return
end

pleClearNaNs(jk)

if(samplesize>0)
    updowntag = 'upper';
    lastchi2 = ar.ple.chi2s{jk}(end);
else
    updowntag = 'lower';
    lastchi2 = ar.ple.chi2s{jk}(1);
end

if(~ar.qFit(jk))
    fprintf('\nPLE#%i extending to %s bound SKIPPED: parameter %s is fixed\n', ...
        jk, updowntag, ar.ple.p_labels{jk});
    return
else
    if(isnan(lastchi2))
        fprintf('\nPLE#%i extending for parameter %s SKIPPED: %s bound already reached?\n', ...
            jk, ar.ple.p_labels{jk}, updowntag);
        return
    else
        fprintf('\nPLE#%i extending for parameter %s to %s bound\n', ...
            jk, ar.ple.p_labels{jk}, updowntag);
    end
end

% respress warning "Matrix is close to singular or badly scaled."
warn_reset = warning;
warning('off', 'MATLAB:nearlySingularMatrix');

%% Algorithm
ar.ple.finished = 0;
if(~ar.ple.breakon_point)
    dchi2 = ar.ple.dchi2;
else
    dchi2 = ar.ple.dchi2_point;
end

arWaitbar(0);

p = ar.ple.p;

if(samplesize>0)
    
    pLast = ar.ple.ps{jk}(end,:);
    feval(ar.ple.integrate_fkt, pLast);
    dpLast = ar.ple.dpLast{jk}(end); 
    dpStep = ar.ple.dpStep{jk}(end,:);
    if dpLast == 0
        dpLast = ar.ple.minstepsize(jk);
        dpStep = NaN;
    end
    %The exception is necessary if the optimum is at the parameter bound,
    %because no step has been made yet.
    
    old_length = length(ar.ple.chi2sinit{jk});
    ar.ple.chi2sinit{jk} = [ar.ple.chi2sinit{jk} nan(1,samplesize)];
    ar.ple.chi2s{jk} = [ar.ple.chi2s{jk} nan(1,samplesize)];
    ar.ple.chi2sviolations{jk} = [ar.ple.chi2sviolations{jk} nan(1,samplesize)];
    ar.ple.chi2spriors{jk} = [ar.ple.chi2spriors{jk} nan(1,samplesize)];
    ar.ple.chi2spriorsAll{jk} = [ar.ple.chi2spriorsAll{jk} nan(1,samplesize)];
    ar.ple.psinit{jk} = [ar.ple.psinit{jk}; nan(samplesize,length(p))];
    ar.ple.psinitstep{jk} = [ar.ple.psinitstep{jk}; nan(samplesize,length(p))];
    ar.ple.ps{jk} = [ar.ple.ps{jk}; nan(samplesize,length(p))];
    ar.ple.gradient{jk} = [ar.ple.gradient{jk}; nan(samplesize,length(p))];
    ar.ple.dpStep{jk} = [ar.ple.dpStep{jk};nan(samplesize,length(p))];
    ar.ple.dpLast{jk} = [ar.ple.dpLast{jk} nan(1,samplesize)];
   
    try
        for j=1:samplesize
            jindex = old_length + j;
            arWaitbar(j, samplesize, sprintf('PLE#%i estimating upper confidence bound for %s', ...
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
    catch exception
        fprintf('ERROR PLE: going to upper bound (%s)\n', exception.message);
    end
else
    samplesize = abs(samplesize);
    
    pLast = ar.ple.ps{jk}(1,:);
    feval(ar.ple.integrate_fkt, pLast);
    dpLast = ar.ple.dpLast{jk}(1); 
    dpStep = ar.ple.dpStep{jk}(1,:);
    if dpLast == 0
        dpLast = -ar.ple.minstepsize(jk);
        dpStep = NaN;
    end
    
    ar.ple.chi2sinit{jk} = [nan(1,samplesize) ar.ple.chi2sinit{jk}];
    ar.ple.chi2s{jk} = [nan(1,samplesize) ar.ple.chi2s{jk}];
    ar.ple.chi2sviolations{jk} = [nan(1,samplesize) ar.ple.chi2sviolations{jk}];
    ar.ple.chi2spriors{jk} = [nan(1,samplesize) ar.ple.chi2spriors{jk}];
    ar.ple.chi2spriorsAll{jk} = [nan(1,samplesize) ar.ple.chi2spriorsAll{jk}];
    ar.ple.psinit{jk} = [nan(samplesize,length(p)); ar.ple.psinit{jk}];
    ar.ple.psinitstep{jk} = [nan(samplesize,length(p)); ar.ple.psinitstep{jk}];
    ar.ple.ps{jk} = [nan(samplesize,length(p)); ar.ple.ps{jk}];
    ar.ple.gradient{jk} = [nan(samplesize,length(p)); ar.ple.gradient{jk}];
    ar.ple.dpStep{jk} = [nan(samplesize,length(p));ar.ple.dpStep{jk}];
    ar.ple.dpLast{jk} = [nan(1,samplesize) ar.ple.dpLast{jk}]; 
       
    try
        for j=1:samplesize
            jindex = samplesize + 1 - j;
            arWaitbar(j, samplesize, sprintf('PLE#%i estimating lower confidence bound for %s', ...
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
            %Determine how far the threshold should be exceeded before
            %termination.
        end
    catch exception
        fprintf('ERROR PLE: going to lower bound (%s)\n', exception.message);
    end
    
end

arWaitbar(-1);

% reset warnings
warning(warn_reset);

% reset parameters
feval(ar.ple.integrate_fkt, ar.ple.p);

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