% fit sequence
%
% arFits(ps, log_fit_history, backup_save)
%
% ps:                           parameter values      
% log_fit_history               [false]
% backup_save                   [false]
% 
% if ps contains rows with only NaN, then this fit is not performed and the
% old fit result is maintained, if existing. This enables overwriting fits
% where integration was not feasible (e.g. in arFitLHS).

function arFits(ps, log_fit_history, backup_save)

global ar

if(~exist('log_fit_history','var'))
    log_fit_history = false;
end
if(~exist('backup_save','var'))
    backup_save = false;
end

dop = find(sum(~isnan(ps),2)>0);

n = length(dop);

if(~isfield(ar,'ps_start'))
    ar.ps_start = ps;
else
    ar.ps_start(dop,:) = ps(dop,:);
end

if(~isfield(ar,'ps'))
    ar.ps = nan(size(ps));
else
    ar.ps(dop,:) = nan(size(ps(dop,:)));
end

if(~isfield(ar,'ps_errors'))
    ar.ps_errors = nan(size(ps));
else
    ar.ps_errors(dop,:) = nan(size(ps(dop,:)));
end

if(~isfield(ar,'chi2s_start'))
    ar.chi2s_start = nan(1,size(ps,1));
else
    ar.chi2s_start(dop) = nan(1,size(ps(dop,:),1));
end

if(~isfield(ar,'chi2sconstr_start'))
    ar.chi2sconstr_start = nan(1,size(ps,1));
else
    ar.chi2sconstr_start(dop) = nan(1,size(ps(dop,:),1));
end
if(~isfield(ar,'chi2s'))
    ar.chi2s = nan(1,size(ps,1));
else
    ar.chi2s(dop) = nan(1,size(ps(dop,:),1));
end
if(~isfield(ar,'chi2sconstr'))
    ar.chi2sconstr = nan(1,size(ps,1));
else
    ar.chi2sconstr(dop) = nan(1,size(ps(dop,:),1));
end
if(~isfield(ar,'exitflag'))
    ar.exitflag = nan(1,size(ps,1));
else
    ar.exitflag(dop) = nan(1,size(ps(dop,:),1));
end
if(~isfield(ar,'timing'))
    ar.timing = nan(1,size(ps,1));
else
    ar.timing(dop) = nan(1,size(ps(dop,:),1));
end
if(~isfield(ar,'fun_evals'))
    ar.fun_evals = nan(1,size(ps,1));
else
    ar.fun_evals(dop) = nan(1,size(ps(dop,:),1));
end
if(~isfield(ar,'optim_crit'))
    ar.optim_crit = nan(1,size(ps,1));
else
    ar.optim_crit(dop) = nan(1,size(ps(dop,:),1));
end

arChi2(true,ar.p(ar.qFit==1));
pReset = ar.p;
chi2Reset = ar.chi2fit + ar.chi2constr;

if(log_fit_history)
    ar.fit_hist = [];
end

arWaitbar(0);
for j=1:n
    arWaitbar(j, n);
    ar.p = ps(dop(j),:);
    if(isfield(ar.config,'useDouble') && ar.config.useDouble==1)
        ar.p(ar.iref) = ar.p(ar.iprimary);
    end
    
    tic;
    try
        arChi2(true, []);
        ar.chi2s_start(dop(j)) = ar.chi2fit;
        ar.chi2sconstr_start(dop(j)) = ar.chi2constr;
        arFit(true);
        ar.ps(dop(j),:) = ar.p;
        ar.chi2s(dop(j)) = ar.chi2fit;
        ar.chi2sconstr(dop(j)) = ar.chi2constr;
        ar.exitflag(dop(j)) = ar.fit.exitflag;
        ar.fun_evals(dop(j)) = ar.fit.fevals;
        ar.optim_crit(dop(j)) = ar.firstorderopt;
    catch exception
        ar.ps_errors(dop(j),:) = ar.p;
        fprintf('fit #%i: %s\n', dop(j), exception.message);
    end
    ar.timing(dop(j)) = toc;
    if(log_fit_history)
        name = ar.config.optimizers{ar.config.optimizer};
        if(ar.config.optimizer==5)
            tmpnames = arNLS;
            name = [name '_' tmpnames{ar.config.optimizerStep+1}]; %#ok<AGROW>
        end
        
        ar.fit_hist(dop(j)).hist = ar.fit;
        ar.fit_hist(dop(j)).optimizer = ar.config.optimizer;
        if(ar.config.optimizer==5)
            ar.fit_hist(dop(j)).optimizerStep = ar.config.optimizerStep;
        else
            ar.fit_hist(dop(j)).optimizerStep = nan;
        end
        ar.fit_hist(dop(j)).config = ar.config.optim;
        ar.fit_hist(dop(j)).name = [name '_' sprintf('run%i', dop(j))];
        
        [~,imin] = min(ar.fit.chi2_hist + ar.fit.constr_hist);
        ar.fit_hist(dop(j)).p = ar.fit.p_hist(imin,:);
    end
    if(backup_save)
        save('arFits_backup.mat','ar');
    end    
end

fprintf('total fitting time: %fsec\n', sum(ar.timing(~isnan(ar.timing))));
fprintf('mean fitting time: %fsec\n', 10^mean(log10(ar.timing(~isnan(ar.timing)))));
arWaitbar(-1);

if(chi2Reset>min(ar.chi2s + ar.chi2sconstr))
    [chi2min,imin] = min(ar.chi2s + ar.chi2sconstr);
    ar.p = ar.ps(imin,:);
    if(ar.config.fiterrors == 1)
        fprintf('selected best fit #%i with %f (old = %f)\n', ...
            imin, 2*ar.ndata*log(sqrt(2*pi)) + chi2min, 2*ar.ndata*log(sqrt(2*pi)) + chi2Reset);
    else
        fprintf('selected best fit #%i with %f (old = %f)\n', ...
            imin, chi2min, chi2Reset);
    end
else
    fprintf('did not find better fit\n');
    ar.p = pReset;
end
arChi2(true,[]);

