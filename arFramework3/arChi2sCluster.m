% chi2 sequence on MATLAB cluster
%
% arChi2sCluster(ps, sensis)
%
% ps:        parameter values      

function arChi2sCluster(ps, sensis, silent)

global ar

if(~exist('sensis','var'))
    sensis = false;
end

if(~exist('silent','var'))
    silent = false;
end

pReset = ar.p;

n = size(ps,1);

ps_errors = nan(n,length(ar.p));
chi2s = nan(1,n);
chi2sconstr = nan(1,n);
timing = nan(1,n);
exitflag = nan(1,n);
fun_evals = nan(1,n);

if(sum(ar.qFit==1)<6)
    figure(1)
    plotmatrix(ps(:,ar.qFit==1), 'x');
end

ar1 = ar;
t1 = tic;

parfor j=1:n
    thisworker = getCurrentWorker; % Worker object
    ar2 = ar1;
    ar2.p = ps(j,:);
    try
        ar2 = arChi2(ar2, sensis, ar2.p(ar2.qFit==1));
        timing(j) = ar2.stop/1e6;
        chi2s(j) = ar2.chi2fit;
        chi2sconstr(j) = ar2.chi2constr;
        exitflag(j) = 1;
        if(~silent) 
            fprintf('feval #%i (%s): objective function %g\n', j, ...
                thisworker.Name, ar2.chi2fit);
        end
    catch exception
        timing(j) = ar2.stop/1e6;
        ps_errors(j,:) = ar2.p;
        if(~silent) 
            fprintf('feval #%i (%s): %s\n', j, ...
                thisworker.Name, exception.message);
        end
        exitflag(j) = -1;
    end
end
toc(t1);

ar.ps = ps;
ar.ps_errors = ps_errors;
ar.chi2s = chi2s;
ar.chi2sconstr = chi2sconstr;
ar.timing = timing;
ar.exitflag = exitflag;
ar.fun_evals = fun_evals;

if(~silent) 
    fprintf('median feval time: %fsec, %i/%i errors\n', median(ar.timing(~isnan(ar.timing) & ...
        ar.exitflag>=0)), ...
        sum(ar.exitflag==-1),n);
end

ar.p = pReset;
try %#ok<TRYNC>
    arChi2(false, ar.p(ar.qFit==1));
end

if(~silent)
    if(sum(ar.qFit==1)<6)
        figure(2)
        plotmatrix(ar.ps(:,ar.qFit==1), 'x');
    end
    
    arPlotChi2s
end
