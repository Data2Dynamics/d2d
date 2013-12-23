% chi2 sequence
%
% arChi2s(ps, sensis)
%
% ps:        parameter values      

function arChi2s(ps, sensis, silent)

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

if(~silent) 
    arWaitbar(0); 
end
ar1 = ar;
parfor j=1:n
    ar2 = ar1;
    if(~silent) 
        arWaitbar(n-j+1, n);
    end
    ar2.p = ps(j,:);
    try
        tic;
        ar2 = arChi2(ar2, sensis, []); 
        timing(j) = ar2.stop;
        chi2s(j) = ar2.chi2fit;
        chi2sconstr(j) = ar2.chi2constr;
        exitflag(j) = 1;
    catch exception
        timing(j) = ar2.stop;
        ps_errors(j,:) = ar2.p;
        if(~silent) 
            fprintf('feval #%i: %s\n', j, exception.message);
        end
        exitflag(j) = -1;
    end
end

ar.ps = ps;
ar.ps_errors = ps_errors;
ar.chi2s = chi2s;
ar.chi2sconstr = chi2sconstr;
ar.timing = timing;
ar.exitflag = exitflag;
ar.fun_evals = fun_evals;

if(~silent) 
    fprintf('mean feval time: %fsec, %i/%i errors\n', mean(ar.timing(~isnan(ar.timing))), ...
        sum(ar.exitflag==-1),n);
    arWaitbar(-1);
end

ar.p = pReset;
try %#ok<TRYNC>
    arChi2(false, []);
end

if(~silent)
    if(sum(ar.qFit==1)<6)
        figure(2)
        plotmatrix(ar.ps(:,ar.qFit==1), 'x');
    end
    
    arPlotChi2s
end
