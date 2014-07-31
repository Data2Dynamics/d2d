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

ar.ps = ps;
ar.ps_errors = nan(n,length(ar.p));
ar.chi2s = nan(1,n);
ar.chi2sconstr = nan(1,n);
ar.timing = nan(1,n);
ar.exitflag = nan(1,n);
ar.fun_evals = nan(1,n);

if(sum(ar.qFit==1)<6)
    figure(1)
    plotmatrix(ps(:,ar.qFit==1), 'x');
end

if(~silent) 
    arWaitbar(0); 
end
t1 = tic;
for j=1:n
    if(~silent) 
        arWaitbar(j, n);
    end
    ar.p = ps(j,:);
    try
        t2 = tic;
        arChi2(sensis, []); 
        ar.timing(j) = ar.stop/1e6;
        if(ar.stop==0)
            ar.timing(j) = toc(t2);
        end
        ar.chi2s(j) = ar.chi2fit;
        ar.chi2sconstr(j) = ar.chi2constr;
        ar.exitflag(j) = 1;
    catch exception
        ar.timing(j) = ar.stop/1e6;
        if(ar.stop==0)
            ar.timing(j) = toc(t2);
        end
        ar.ps_errors(j,:) = ar.p;
        if(~silent) 
            fprintf('feval #%i: %s\n', j, exception.message);
        end
        ar.exitflag(j) = -1;
    end
end
toc(t1);

if(~silent) 
    fprintf('median feval time: %fsec, %i/%i errors\n', median(ar.timing(~isnan(ar.timing) & ...
        ar.exitflag>=0)), ...
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
