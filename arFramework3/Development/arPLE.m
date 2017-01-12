function arPLE(jk)

global ar

if(~exist('jk','var') || isempty(jk))
    jk = find(ar.qFit==1);
end
if(length(jk)>1)
    for j=1:length(jk)
        arPLE(jk(j));
    end
    return;
end

n = length(ar.scan.chi2s{jk});

if(~isfield(ar, 'ple'))
    ar.ple.chi2s = {};
    ar.ple.errors = {};
    ar.ple.ps = {};
    ar.ple.run = zeros(size(ar.p));
    ar.ple.alpha = 0.05;
    ar.ple.ndof = 1;
end

ar.ple.chi2s{jk} = nan(1,n);
ar.ple.errors{jk} = nan(1,n);
ar.ple.ps{jk} = nan(n,length(ar.p));
ar.ple.run(jk) = 1;

pReset = ar.p;
qFitReset = ar.qFit(jk);

ar.qFit(jk) = 0;

fprintf('PLE for #%i %s...', jk, ar.pLabel{jk});
arWaitbar(0);
tic;
for j=1:n
    arWaitbar(j,n, sprintf('PLE for %s', strrep(ar.pLabel{jk},'_','\_')));
    
    ar.p = pReset;
    ar.p(jk) = ar.scan.ps{jk}(j);
    try
        arFit(true);
        
        ar.ple.ps{jk}(j,:) = ar.p;
        ar.ple.chi2s{jk}(j) = ar.chi2fit;
        ar.ple.errors{jk}(j) = ar.fit.exitflag;
    catch errorid
        disp(errorid.message);
        ar.ple.errors{jk}(j) = -5;
    end
end
arWaitbar(-1);
fprintf('(%i errors, %i fit issues, %s elapse time)\n', ...
    sum(ar.ple.errors{jk}<0), ...
    sum(ar.ple.errors{jk}>1 | ar.ple.errors{jk}==0), ...
    secToHMS(toc));

ar.p = pReset;
ar.qFit(jk) = qFitReset;
arCalcMerit(false);
