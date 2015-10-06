% fit sequence on MATLAB cluster
%
% arFitsCluster(ps, log_fit_history, backup_save)
%
% ps:                           parameter values      
% log_fit_history               [false]
% backup_save                   [false]

function arFitsCluster(ps, log_fit_history, backup_save)

global ar

if(~exist('log_fit_history','var'))
    log_fit_history = false;
end
if(~exist('backup_save','var'))
    backup_save = false;
end

n = size(ps,1);
ar.ps_start = ps;
ps_errors = nan(size(ps));
chi2s_start = nan(1,n);
chi2sconstr_start = nan(1,n);
chi2s = nan(1,n);
chi2sconstr = nan(1,n);
exitflag = nan(1,n);
timing = nan(1,n);
fun_evals = nan(1,n);
optim_crit = nan(1,n);

arChi2(true,ar.p(ar.qFit==1));
pReset = ar.p;
chi2Reset = ar.chi2fit + ar.chi2constr;

if(log_fit_history)
    ar.fit_hist = [];
end

if(backup_save && exist('arFitCluster_Backup','dir')~=7)
    mkdir('arFitCluster_Backup');
end

ar1 = ar;
pct = parfor_progress(n); %#ok<NASGU>
parfor j=1:n
    thisworker = getCurrentWorker; % Worker object
    
    % load from backup if exist
    if(backup_save && exist(sprintf('./arFitCluster_Backup/fit_%i.mat', j),'file'))
        S = load(sprintf('./arFitCluster_Backup/fit_%i.mat', j));
        chi2s_start(j) = S.x.chi2fit;
        chi2sconstr_start(j) = S.x.chi2constr;
        ps(j,:) = S.x.p;
        chi2s(j) = S.x.chi2fit;
        chi2sconstr(j) = S.x.chi2constr;
        exitflag(j) = S.x.exitflag;
        fun_evals(j) = S.x.fevals;
        optim_crit(j) = S.x.firstorderopt;
        timing(j) = S.x.timing;
        pct = parfor_progress/100;
        fprintf('%i/%i fit #%i (%s, %s): loaded from backup\n', round(n*pct), n, j, ...
            thisworker.Name, datestr(now, 0));
        continue
    end
    
    ar2 = ar1;
    ar2.p = ps(j,:);
    tic;
    try
        ar2 = arChi2(ar2, true, ar2.p(ar2.qFit==1));
        chi2s_start(j) = ar2.chi2fit;
        chi2sconstr_start(j) = ar2.chi2constr;
        ar2 = arFit(ar2, true);
        ps(j,:) = ar2.p;
        chi2s(j) = ar2.chi2fit;
        chi2sconstr(j) = ar2.chi2constr;
        exitflag(j) = ar2.fit.exitflag;
        fun_evals(j) = ar2.fit.fevals;
        optim_crit(j) = ar2.firstorderopt;
        pct = parfor_progress/100;
        if(ar2.config.fiterrors == 1)
            fprintf('%i/%i fit #%i (%s, %s): objective function %g\n', round(n*pct), n, j, ...
                thisworker.Name, datestr(now, 0), ...
                2*ar2.ndata*log(sqrt(2*pi)) + ar2.chi2fit + ar2.chi2constr);
        else
            fprintf('%i/%i fit #%i (%s, %s): objective function %g\n', round(n*pct), n, j, ...
                thisworker.Name, datestr(now, 0), ...
                ar2.chi2fit + ar2.chi2constr);
        end
        
        if(log_fit_history)
            name = ar2.config.optimizers{ar2.config.optimizer};
            if(ar2.config.optimizer==5)
                tmpnames = arNLS;
                name = [name '_' tmpnames{ar2.config.optimizerStep+1}];
            end
            
            fit_hist(j).hist = ar2.fit;
            fit_hist(j).optimizer = ar2.config.optimizer;
            if(ar2.config.optimizer==5)
                fit_hist(j).optimizerStep = ar2.config.optimizerStep;
            else
                fit_hist(j).optimizerStep = nan;
            end
            fit_hist(j).config = ar2.config.optim;
            fit_hist(j).name = [name '_' sprintf('run%i', j)];
            
            [~,imin] = min(ar2.fit.chi2_hist + ar2.fit.constr_hist);
            fit_hist(j).p = ar2.fit.p_hist(imin,:);
        end
    catch exception
        ps_errors(j,:) = ar2.p;
        pct = parfor_progress/100;
        fprintf('%i/%i fit #%i (%s, %s): %s\n', round(n*pct), n, j, ...
            thisworker.Name, datestr(now, 0), ...
            exception.message);
    end
    timing(j) = toc;
    
    % save to backup if not exist
    if(backup_save && ~exist(sprintf('./arFitCluster_Backup/fit_%i.mat', j),'file'))
        ar3 = struct([]);
        ar3(1).chi2fit = chi2s_start(j);
        ar3.chi2constr = chi2sconstr_start(j);
        ar3.p = ps(j,:);
        ar3.chi2fit = chi2s(j);
        ar3.chi2constr = chi2sconstr(j);
        ar3.exitflag = exitflag(j);
        ar3.fevals = fun_evals(j);
        ar3.firstorderopt = optim_crit(j);
        ar3.timing = timing(j);
        parsave(sprintf('./arFitCluster_Backup/fit_%i.mat', j), ar3)
    end
end
pct = parfor_progress(0); %#ok<NASGU>

ar.chi2s_start = chi2s_start;
ar.chi2sconstr_start = chi2sconstr_start;
ar.ps = ps;
ar.chi2s = chi2s;
ar.chi2sconstr = chi2sconstr;
ar.exitflag = exitflag;
ar.fun_evals = fun_evals;
ar.optim_crit = optim_crit;
ar.ps_errors = ps_errors;
ar.timing = timing;
if(log_fit_history)
    ar.fit_hist = fit_hist;
end

fprintf('total fitting time: %fsec\n', sum(ar.timing(~isnan(ar.timing))));
fprintf('median fitting time: %fsec\n', median(ar.timing(~isnan(ar.timing))));

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
arChi2(true,ar.p(ar.qFit==1));

%% backup save on cluster before transfer
save('arFitsCluster_backup.mat', 'ar');

