% Profile Likelihood Exploit
%
% pleExtend(i, samplesize)
%
% i:                    i'th parameter, see pwInfo
%                       [if omitted, all free parameters are considered]
% samplesize:           number of sampling steps            [100]

function pleExtend(jk, samplesize)

global pleGlobals;

if(isempty(pleGlobals))
    error('PLE ERROR: please initialize')
end 
if(~isfield(pleGlobals, 'showCalculation'))
    pleGlobals.showCalculation = true;
end

if(nargin == 0 || isempty(jk))
    if(~exist('samplesize','var'))
        samplesize = 100;
    end
    fprintf('PLE extending for %i parameters ...\n', sum(pleGlobals.q_fit))
    jindex = find(pleGlobals.q_fit(1:length(pleGlobals.chi2s)) & ~cellfun(@isempty,pleGlobals.chi2s));
    do_plotting = pleGlobals.showCalculation;
    pleGlobals.showCalculation = false;
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
    pleGlobals.showCalculation = do_plotting;
    return
end

if(nargin == 1)
    pleExtend(jk, +100);
    pleExtend(jk, -100);
    return
end

if(samplesize>0)
    updowntag = 'upper';
    lastchi2 = pleGlobals.chi2s{jk}(end);
else
    updowntag = 'lower';
    lastchi2 = pleGlobals.chi2s{jk}(1);
end

if(~pleGlobals.q_fit(jk))
    fprintf('\nPLE#%i extending to %s bound SKIPPED: parameter %s is fixed\n', ...
        jk, updowntag, pleGlobals.p_labels{jk});
    return
else
    if(isnan(lastchi2))
        fprintf('\nPLE#%i extending for parameter %s SKIPPED: %s bound already reached?\n', ...
            jk, pleGlobals.p_labels{jk}, updowntag);
        return
    else
        fprintf('\nPLE#%i extending for parameter %s to %s bound\n', ...
            jk, pleGlobals.p_labels{jk}, updowntag);
    end
end

% respress warning "Matrix is close to singular or badly scaled."
warn_reset = warning;
warning('off', 'MATLAB:nearlySingularMatrix');

pleGlobals.conf_labels = {'Hessian', 'PLE'};

%% Algorithm

if(~pleGlobals.breakon_point)
    dchi2 = pleGlobals.dchi2;
else
    dchi2 = pleGlobals.dchi2_point;
end

estimatetime = 0;
fittime = 0;
arWaitbar(0);

p = pleGlobals.p;

if(samplesize>0)
    
    pLast = pleGlobals.ps{jk}(end,:);
    feval(pleGlobals.integrate_fkt, pLast);
    dpLast = pleGlobals.maxstepsize(jk)/2.1; % exploiting upper bound, base step
    
    old_length = length(pleGlobals.chi2sinit{jk});
    pleGlobals.chi2sinit{jk} = [pleGlobals.chi2sinit{jk} nan(1,samplesize)];
    pleGlobals.chi2s{jk} = [pleGlobals.chi2s{jk} nan(1,samplesize)];
    pleGlobals.chi2sviolations{jk} = [pleGlobals.chi2sviolations{jk} nan(1,samplesize)];
    pleGlobals.chi2spriors{jk} = [pleGlobals.chi2spriors{jk} nan(1,samplesize)];
    pleGlobals.psinit{jk} = [pleGlobals.psinit{jk}; nan(samplesize,length(p))];
    pleGlobals.psinitstep{jk} = [pleGlobals.psinitstep{jk}; nan(samplesize,length(p))];
    pleGlobals.ps{jk} = [pleGlobals.ps{jk}; nan(samplesize,length(p))];
    pleGlobals.gradient{jk} = [pleGlobals.gradient{jk}; nan(samplesize,length(p))];
    
    try
        for j=1:samplesize
            jindex = old_length + j;
            arWaitbar(j, samplesize, sprintf('PLE#%i estimating upper confidence bound for %s', ...
                jk, strrep(pleGlobals.p_labels{jk},'_', '\_')));
            
            tic;
            % estimate intial step
            [pStep, dpLast] = feval(pleGlobals.initstep_fkt, jk, pLast, dpLast);
            if(sum(isnan(pStep))>0)
                break;
            end
            pLast = pLast + pStep;
            
            feval(pleGlobals.integrate_fkt, pLast);
            pleGlobals.chi2sinit{jk}(jindex) = feval(pleGlobals.merit_fkt);
            pleGlobals.psinit{jk}(jindex,:) = pLast;
            pleGlobals.psinitstep{jk}(jindex,:) = pStep;
            
            estimatetime = estimatetime + toc;
            
            tic;
            % Fit
            [p, gradient] = feval(pleGlobals.fit_fkt, jk);
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
                plePlot(jk);
            end
            
            if(feval(pleGlobals.merit_fkt) > pleGlobals.chi2+dchi2*1.2)
                break
            end
        end
    catch exception
        fprintf('ERROR PLE: going to upper bound (%s)\n', exception.message);
    end
    
else
    samplesize = abs(samplesize);
    
    pLast = pleGlobals.ps{jk}(1,:);
    feval(pleGlobals.integrate_fkt, pLast);
    dpLast = -pleGlobals.maxstepsize(jk)/2.1; % exploiting lower bound, base step
    
    pleGlobals.chi2sinit{jk} = [nan(1,samplesize) pleGlobals.chi2sinit{jk}];
    pleGlobals.chi2s{jk} = [nan(1,samplesize) pleGlobals.chi2s{jk}];
    pleGlobals.chi2sviolations{jk} = [nan(1,samplesize) pleGlobals.chi2sviolations{jk}];
    pleGlobals.chi2spriors{jk} = [nan(1,samplesize) pleGlobals.chi2spriors{jk}];
    pleGlobals.psinit{jk} = [nan(samplesize,length(p)); pleGlobals.psinit{jk}];
    pleGlobals.psinitstep{jk} = [nan(samplesize,length(p)); pleGlobals.psinitstep{jk}];
    pleGlobals.ps{jk} = [nan(samplesize,length(p)); pleGlobals.ps{jk}];
    pleGlobals.gradient{jk} = [nan(samplesize,length(p)); pleGlobals.gradient{jk}];
    
    try
        for j=1:samplesize
            jindex = samplesize + 1 - j;
            arWaitbar(j, samplesize, sprintf('PLE#%i estimating lower confidence bound for %s', ...
                jk, strrep(pleGlobals.p_labels{jk},'_', '\_')));
            
            tic;
            % estimate intial step
            [pStep, dpLast] = feval(pleGlobals.initstep_fkt, jk, pLast, dpLast);
            if(sum(isnan(pStep))>0)
                break;
            end
            pLast = pLast + pStep;
            
            feval(pleGlobals.integrate_fkt, pLast);
            pleGlobals.chi2sinit{jk}(jindex) = feval(pleGlobals.merit_fkt);
            pleGlobals.psinit{jk}(jindex,:) = pLast;
            pleGlobals.psinitstep{jk}(jindex,:) = pStep;
            
            estimatetime = estimatetime + toc;
            
            tic;
            % Fit
            [p, gradient] = feval(pleGlobals.fit_fkt, jk);
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
                plePlot(jk);
            end
            
            if(feval(pleGlobals.merit_fkt) > pleGlobals.chi2+dchi2*1.2)
                break
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
feval(pleGlobals.integrate_fkt, pleGlobals.p);

%% Finalizing

chi2 = pleGlobals.chi2s{jk};
ps = pleGlobals.ps{jk};

pleGlobals.estimatetime(jk) = pleGlobals.estimatetime(jk) + estimatetime;
pleGlobals.fittime(jk) = pleGlobals.fittime(jk) + fittime;
pleGlobals.timing(jk) = pleGlobals.timing(jk) + fittime+estimatetime;

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

