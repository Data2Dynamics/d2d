function arPLE_direct(jk, n)

global ar

if(~exist('jk','var') || isempty(jk))
    jk = find(ar.qFit==1);
end
if(length(jk)>1)
    for j=1:length(jk)
        arPLE_direct(jk(j));
    end
    return;
end
if(~exist('n','var'))
    n = 50;
end

if(~isfield(ar, 'ple'))
    ar.ple.chi2s = {};
    ar.ple.errors = {};
    ar.ple.ps = {};
    ar.ple.run = zeros(size(ar.p));
    ar.ple.alpha = 0.05;
    ar.ple.ndof = 1;
end

arChi2(true);

chi2Reset = ar.chi2fit;
pReset = ar.p;
ar.ple.pStart = ar.p;
qFitReset = ar.qFit(jk);

fprintf('PLE #%i for %s...\n', jk, ar.pLabel{jk});
arWaitbar(0);
tic;

ar.qFit(jk) = 0;

[chi2sup, psup, errorsup] = ple_task(jk, n, 1, chi2Reset, pReset);
[chi2sdown, psdown, errorsdown] = ple_task(jk, n, -1, chi2Reset, pReset);

ar.ple.chi2s{jk} = [fliplr(chi2sdown) chi2Reset chi2sup];
ar.ple.ps{jk} = [flipud(psdown); pReset; psup];
ar.ple.errors{jk} = [fliplr(errorsdown) nan errorsup];
ar.ple.run(jk) = 1;

arWaitbar(-1);
fprintf('PLE #%i %i errors, %i fit issues, %s elapse time\n', jk, ...
    sum(ar.ple.errors{jk}<0), ...
    sum(ar.ple.errors{jk}==0), ...
    secToHMS(toc));

ar.p = pReset;
ar.qFit(jk) = qFitReset;
arChi2(false);



function [chi2s, ps, errors] = ple_task(jk, n, direction, chi2Reset, pReset)

global ar

chi2s = nan(1,n);
ps = nan(n,length(pReset));
errors = nan(1,n);

dchi2 = chi2inv(1-ar.ple.alpha, ar.ple.ndof);
rel_increase = 0.5;

dp = 0.1;
dp_min = 1e-6;
ar.p = pReset;

chi2Last = chi2Reset;
pLast = pReset;

for j=1:n
    if(direction>0)
        arWaitbar(j,2*n, sprintf('PLE up for %s', strrep(ar.pLabel{jk},'_','\_')));
    else
        arWaitbar(j+n,2*n, sprintf('PLE down for %s', strrep(ar.pLabel{jk},'_','\_')));
    end
   
    arChi2(true);
    pTrial = findTrial(jk, pLast, dp, direction);
    
    try
        ar.p = pTrial;
        arChi2(false);
        
        while(ar.chi2fit - chi2Last > rel_increase * dchi2)
            dp = dp/2;
            
            if(dp < dp_min)
                error('minimum step size reached'); 
            end
            
            pTrial = findTrial(jk, pLast, dp, direction);
            
            ar.p = pTrial;
            arChi2(false);
        end
        dp
        
        arFit(true);
        
        ps(j,:) = ar.p;
        chi2s(j) = ar.chi2fit;
        errors(j) = ar.fit.exitflag;
    catch errorid
        fprintf('PLE #%i ERROR %s\n', jk, errorid.message); 
        errors(j) = -5;
        break;
    end
    
    if(ar.chi2fit - chi2Reset > 2*dchi2)
        fprintf('PLE #%i reached confidence limit\n', jk); 
        break;
    end
    if(pTrial(jk) >= ar.ub(jk))
        fprintf('PLE #%i reached upper parameter bound\n', jk);
        break;
    elseif(pTrial(jk) <= ar.lb(jk))
        fprintf('PLE #%i reached lower parameter bound\n', jk);
        break;
    end
    
    pLast = ar.p;
    chi2Last = ar.chi2fit;
    
    dp = dp * 2;
end

ar.p = pReset;



function pTrial = findTrial(jk, pLast, dp, direction)

global ar

pTrial = pLast;
pTrial(jk) = pLast(jk) + dp*direction;

if(pTrial(jk) > ar.ub(jk))
    pTrial(jk) = ar.ub(jk);
elseif(pTrial(jk) < ar.lb(jk))
    pTrial(jk) = ar.lb(jk);
end

return

beta = transpose(-ar.res*ar.sres);
alpha = ar.sres'*ar.sres;

qFit = ar.qFit==1;
alpha2 = alpha(qFit,qFit);
beta2 = beta(qFit) + alpha(qFit,jk)*(dp*direction);

lb = zeros(size(beta2)) - dp;
ub = zeros(size(beta2)) + dp;

opts = optimset('Display','off');
dpTrial = quadprog(alpha2, beta2, [], [], [], [], lb, ub, [], opts);
% dpTrial = quadprog(alpha2, -beta2, [], [], [], [], lb, ub);

pTrial(qFit) = pTrial(qFit) + dpTrial';
