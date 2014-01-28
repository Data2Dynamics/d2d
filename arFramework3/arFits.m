% fit sequence
%
% arFits(ps, log_fit_history, backup_save)
%
% ps:                           parameter values      
% log_fit_history               [false]
% backup_save                   [false]

function arFits(ps, log_fit_history, backup_save)

global ar

if(~exist('log_fit_history','var'))
    log_fit_history = false;
end
if(~exist('backup_save','var'))
    backup_save = false;
end

n = size(ps,1);
ar.ps_start = ps;
ar.ps = nan(size(ps));
ar.ps_errors = nan(size(ps));
ar.chi2s_start = nan(1,n);
ar.chi2sconstr_start = nan(1,n);
ar.chi2s = nan(1,n);
ar.chi2sconstr = nan(1,n);
ar.exitflag = nan(1,n);
ar.timing = nan(1,n);
ar.fun_evals = nan(1,n);
ar.optim_crit = nan(1,n);

arChi2(true,[]);
pReset = ar.p;
chi2Reset = ar.chi2fit + ar.chi2constr;

if(log_fit_history)
    ar.fit_hist = [];
end

arWaitbar(0);
for j=1:n
    arWaitbar(j, n);
    ar.p = ps(j,:);
    tic;
    try
        arChi2(true, []);
        ar.chi2s_start(j) = ar.chi2fit;
        ar.chi2sconstr_start(j) = ar.chi2constr;
        arFit(true);
        ar.ps(j,:) = ar.p;
        ar.chi2s(j) = ar.chi2fit;
        ar.chi2sconstr(j) = ar.chi2constr;
        ar.exitflag(j) = ar.fit.exitflag;
        ar.fun_evals(j) = ar.fit.fevals;
        ar.optim_crit(j) = ar.firstorderopt;
    catch exception
        ar.ps_errors(j,:) = ar.p;
        fprintf('fit #%i: %s\n', j, exception.message);
    end
    ar.timing(j) = toc;
    if(log_fit_history)
        name = ar.config.optimizers{ar.config.optimizer};
        if(ar.config.optimizer==5)
            tmpnames = arNLS;
            name = [name '_' tmpnames{ar.config.optimizerStep+1}]; %#ok<AGROW>
        end
        
        ar.fit_hist(j).hist = ar.fit;
        ar.fit_hist(j).optimizer = ar.config.optimizer;
        if(ar.config.optimizer==5)
            ar.fit_hist(j).optimizerStep = ar.config.optimizerStep;
        else
            ar.fit_hist(j).optimizerStep = nan;
        end
        ar.fit_hist(j).config = ar.config.optim;
        ar.fit_hist(j).name = [name '_' sprintf('run%i', j)];
        
        [~,imin] = min(ar.fit.chi2_hist + ar.fit.constr_hist);
        ar.fit_hist(j).p = ar.fit.p_hist(imin,:);
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

