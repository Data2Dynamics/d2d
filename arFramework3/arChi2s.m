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
ar.ps_errors = [];
ar.chi2s = nan(1,n);
ar.timing = nan(1,n);
ar.exitflag = nan(1,n);
ar.fun_evals = [];

if(sum(ar.qFit==1)<6)
    figure(1)
    plotmatrix(ps(:,ar.qFit==1), 'x');
end

arWaitbar(0);
for j=1:n
    arWaitbar(j, n);
    ar.p = ps(j,:);
    tic;
    try
        arChi2(sensis);
        ar.timing(j) = toc;
        ar.chi2s(j) = ar.chi2fit;
        ar.exitflag(j) = 1;
    catch exception
        ar.timing(j) = toc;
        ar.ps_errors(end+1,:) = ar.p;
        fprintf('feval #%i: %s\n', j, exception.message);
        ar.exitflag(j) = -1;
    end
end
fprintf('mean feval time: %fsec, %i/%i errors\n', mean(ar.timing(~isnan(ar.timing))), ...
    sum(ar.exitflag==-1),n);
arWaitbar(-1);

ar.p = pReset;
arChi2(false);

if(~silent)
    if(sum(ar.qFit==1)<6)
        figure(2)
        plotmatrix(ar.ps(:,ar.qFit==1), 'x');
    end
    
    arPlotChi2s
end
