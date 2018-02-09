% Profile Likelihood Exploit
%
% pleExtend(i, samplesize)
%
% i:                    i'th parameter, see pwInfo
%                       [if omitted, all free parameters are considered]
% samplesize:           number of sampling steps            [100]

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

ar.ple.conf_labels = {'Hessian', 'PLE'};

%% Algorithm
ar.ple.finished = 0;
if(~ar.ple.breakon_point)
    dchi2 = ar.ple.dchi2;
else
    dchi2 = ar.ple.dchi2_point;
end

estimatetime = 0;
fittime = 0;
arWaitbar(0);

p = ar.ple.p;

if(samplesize>0)
    
    pLast = ar.ple.ps{jk}(end,:);
    feval(ar.ple.integrate_fkt, pLast);
    dpLast = ar.ple.maxstepsize(jk)/2.1; % exploiting upper bound, base step
    
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
    
    attempts = 0;
    try
        for j=1:samplesize
            jindex = old_length + j;
            arWaitbar(j, samplesize, sprintf('PLE#%i estimating upper confidence bound for %s', ...
                jk, strrep(ar.ple.p_labels{jk},'_', '\_')));
            
            tic;
            % estimate intial step
            [pStep, dpLast] = feval(ar.ple.initstep_fkt, jk, pLast, dpLast);
            if(sum(isnan(pStep))>0)
                break;
            end
            pLast = pLast + pStep;
            
            feval(ar.ple.integrate_fkt, pLast);
            ar.ple.chi2sinit{jk}(jindex) = feval(ar.ple.merit_fkt);
            ar.ple.psinit{jk}(jindex,:) = pLast;
            ar.ple.psinitstep{jk}(jindex,:) = pStep;
            
            estimatetime = estimatetime + toc;
            
            tic;
            % Fit
            try
                [p, gradient] = feval(ar.ple.fit_fkt, jk);
                pLast = p;
                
                ar.ple.ps{jk}(jindex,:) = pLast;
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
                fittime = fittime + toc;

                if(ar.ple.showCalculation)
                    plePlot(jk);
                end
                if(isfield(ar.ple, 'continuousSave') && ar.ple.continuousSave)
                    pleSave(ar.ple);
                end

                if(feval(ar.ple.merit_fkt) > ar.ple.merit+dchi2*1.2)
                    break
                end
                
                % Reset attempts
                attempts = 0;
                
            catch exception
                % Do a few more attempts
                attempts = attempts + 1;
                if ( attempts > ar.ple.attempts )
                    rethrow(exception);
                end
            end
        end
    catch exception
        fprintf('ERROR PLE: going to upper bound (%s)\n', exception.message);
    end
else
    samplesize = abs(samplesize);
    
    pLast = ar.ple.ps{jk}(1,:);
    feval(ar.ple.integrate_fkt, pLast);
    dpLast = -ar.ple.maxstepsize(jk)/2.1; % exploiting lower bound, base step
    
    ar.ple.chi2sinit{jk} = [nan(1,samplesize) ar.ple.chi2sinit{jk}];
    ar.ple.chi2s{jk} = [nan(1,samplesize) ar.ple.chi2s{jk}];
    ar.ple.chi2sviolations{jk} = [nan(1,samplesize) ar.ple.chi2sviolations{jk}];
    ar.ple.chi2spriors{jk} = [nan(1,samplesize) ar.ple.chi2spriors{jk}];
    ar.ple.chi2spriorsAll{jk} = [nan(1,samplesize) ar.ple.chi2spriorsAll{jk}];
    ar.ple.psinit{jk} = [nan(samplesize,length(p)); ar.ple.psinit{jk}];
    ar.ple.psinitstep{jk} = [nan(samplesize,length(p)); ar.ple.psinitstep{jk}];
    ar.ple.ps{jk} = [nan(samplesize,length(p)); ar.ple.ps{jk}];
    ar.ple.gradient{jk} = [nan(samplesize,length(p)); ar.ple.gradient{jk}];
    
    attempts = 0;
    try
        for j=1:samplesize
            jindex = samplesize + 1 - j;
            arWaitbar(j, samplesize, sprintf('PLE#%i estimating lower confidence bound for %s', ...
                jk, strrep(ar.ple.p_labels{jk},'_', '\_')));
            
            tic;
            % estimate intial step
            [pStep, dpLast] = feval(ar.ple.initstep_fkt, jk, pLast, dpLast);
            if(sum(isnan(pStep))>0)
                break;
            end
            pLast = pLast + pStep;
            
            feval(ar.ple.integrate_fkt, pLast);
            ar.ple.chi2sinit{jk}(jindex) = feval(ar.ple.merit_fkt);
            ar.ple.psinit{jk}(jindex,:) = pLast;
            ar.ple.psinitstep{jk}(jindex,:) = pStep;
            
            estimatetime = estimatetime + toc;
            
            tic;
            % Fit           
            try
                [p, gradient] = feval(ar.ple.fit_fkt, jk);
                pLast = p;
            
                ar.ple.ps{jk}(jindex,:) = pLast;
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
                fittime = fittime + toc;

                if(ar.ple.showCalculation)
                    plePlot(jk);
                end
                if(isfield(ar.ple, 'continuousSave') && ar.ple.continuousSave)
                    pleSave(ar.ple);
                end

                if(feval(ar.ple.merit_fkt) > ar.ple.merit+dchi2*1.2)
                    break
                end
                
                % Reset attempts
                attempts = 0;
                
            catch exception
                % Do a few more attempts
                attempts = attempts + 1;
                if ( attempts > ar.ple.attempts )
                    rethrow(exception);
                end
            end
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

ar.ple.estimatetime(jk) = ar.ple.estimatetime(jk) + estimatetime;
ar.ple.fittime(jk) = ar.ple.fittime(jk) + fittime;
ar.ple.timing(jk) = ar.ple.timing(jk) + fittime+estimatetime;

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
        ar.ple.conf_lb_point(jk) = interp1(chi2([kind kind-1]), ps([kind kind-1], jk), min(chi2)+ar.ple.dchi2_point);
    end
    if(max(ps(q_chi2good_point,jk))==max(ps(~isnan(chi2),jk)))
        ar.ple.conf_ub_point(jk) = Inf;
    else
        kind = find(ps(:,jk)==max(ps(q_chi2good_point,jk)));
        ar.ple.conf_ub_point(jk) = interp1(chi2([kind kind+1]), ps([kind kind+1], jk), min(chi2)+ar.ple.dchi2_point);
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
    ar.ple.conf_rel = abs((ar.ple.conf_ub - ar.ple.conf_lb)/2*100 ./ ar.ple.p);
    ar.ple.conf_rel_point = abs((ar.ple.conf_ub_point - ar.ple.conf_lb_point)/2*100 ./ ar.ple.p);

    % calulate coverage
    if(isfield(ar.ple, 'p_true'))
        ar.ple.cover = ar.ple.p_true <= ar.ple.conf_ub & ...
            ar.ple.p_true >= ar.ple.conf_lb;
        ar.ple.cover_point = ar.ple.p_true <= ar.ple.conf_ub_point & ...
            ar.ple.p_true >= ar.ple.conf_lb_point;
    end
    
    % calculate elapse times
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
pleSave(ar.ple)

