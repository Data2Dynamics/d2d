% arLoadFitBackup([n])
% loads existing backups of arFits function and saves them to one workspace (arFitsCluster_backup.mat)
%
% n:    number of fits   [size(ps,1)]

function arLoadFitBackup(n)

global ar

if isfield(ar,'ps_start')
    ps = ar.ps_start;
else
    ps = arRandomPars(n, []);
    ar.ps_start = ps;
end
if(~exist('n','var') || isempty(n))
    n = size(ps,1);
end

ps_errors = nan(size(ps));
chi2s_start = nan(1,n);
chi2sconstr_start = nan(1,n);
chi2s = nan(1,n);
chi2sconstr = nan(1,n);
exitflag = nan(1,n);
timing = nan(1,n);
fun_evals = nan(1,n);
optim_crit = nan(1,n);


startTime = clock;
arShowProgressParFor(n);

for j=1:n
    % load from backup if exist
    if (exist(sprintf('./arFitCluster_Backup/fit_%i.mat', j),'file'))
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
end

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

save('arFitsCluster_backup.mat', 'ar');
