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
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
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

startTime = clock;
arShowProgressParFor(n);

ar1 = ar;
parfor j=1:n
    
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
        
        arShowProgressParFor(j, n, startTime, 'loaded from backup')
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
        
        if(ar2.config.useFitErrorMatrix==0 && ar2.config.fiterrors == 1)
            objective_functino_val = 2*ar2.ndata*log(sqrt(2*pi)) + ar2.chi2fit + ar2.chi2constr;
        elseif(ar2.config.useFitErrorMatrix==1 && sum(sum(ar2.config.fiterrors_matrix==1))>0)
            objective_functino_val = 2*ar2.ndata_err*log(sqrt(2*pi)) + ar2.chi2fit + ar2.chi2constr;
        else
            objective_functino_val = ar2.chi2fit + ar2.chi2constr;
        end
        return_message = sprintf('objective function %g', objective_functino_val);
        
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
        return_message = exception.message;
    end
    
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
    
    arShowProgressParFor(j, n, startTime, return_message)
    timing(j) = toc;
end
arShowProgressParFor(0);

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

fprintf('total fitting time: %s\n', secToHMS(sum(ar.timing(~isnan(ar.timing)))));
fprintf('median fitting time: %s\n', secToHMS(median(ar.timing(~isnan(ar.timing)))));

if(chi2Reset>min(ar.chi2s + ar.chi2sconstr))
    [chi2min,imin] = min(ar.chi2s + ar.chi2sconstr);
    ar.p = ar.ps(imin,:);
    if(ar.config.useFitErrorMatrix==0 && ar.config.fiterrors == 1)
        fprintf('selected best fit #%i with %f (old = %f)\n', ...
            imin, 2*ar.ndata*log(sqrt(2*pi)) + chi2min, 2*ar.ndata*log(sqrt(2*pi)) + chi2Reset);
    elseif(ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0)
        fprintf('selected best fit #%i with %f (old = %f)\n', ...
            imin, 2*ar.ndata_err*log(sqrt(2*pi)) + chi2min, 2*ar.ndata_err*log(sqrt(2*pi)) + chi2Reset);
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
% Sometime the automatic data transfer from the cluster does not work.
% So it is better to save a local copy of the final results that 
% can be retrieved by file transfer.

save('arFitsCluster_backup.mat', 'ar');


