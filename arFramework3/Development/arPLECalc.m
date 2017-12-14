% profile likelilhood calculation
%
%   arPLE_direct(jk, n)
%
%   jk: parameter index or indices
%   n:  number of ple steps up and down

function varargout = arPLECalc(varargin)

global ar

if(nargin==0)
    qglobalar = true;
    jk = find(ar.qFit==1);
    n = 50;
else
    if(isstruct(varargin{1}))
        qglobalar = false;
        ar = varargin{1};
        if(nargin==1)
            jk = find(ar.qFit==1);
            n = 50;
        elseif(nargin==2)
            jk = varargin{2};
            n = 50;
        elseif(nargin==3)
            jk = varargin{2};
            n = varargin{3};
        end
    else
        qglobalar = true;
        if(nargin==1)
            jk = varargin{1};
            n = 50;
        elseif(nargin==2)
            jk = varargin{1};
            n = varargin{2};
        end
    end
end

if(length(jk)>1)
    for j=1:length(jk)
        if(qglobalar)
            arPLECalc(jk(j));
        else
            ar = arPLECalc(ar, jk(j), n);
        end
    end
    return;
end

if(~isfield(ar, 'ple'))
    ar.ple.chi2s = {};
    ar.ple.errors = {};
    ar.ple.ps = {};
    ar.ple.run = zeros(size(ar.p));
    ar.ple.chi2Reset = nan(size(ar.p));
    ar.ple.alpha = 0.05;
    ar.ple.ndof = 1;
end
if(~isfield(ar.ple, 'p_labels'))
    ar.ple.p_labels = ar.pLabel;
end
if(~isfield(ar.ple, 'alpha_level'))
    ar.ple.alpha_level = ar.ple.alpha;
end

% save original parameters
pReset = ar.p;
ar.ple.pStart = ar.p;
ar.ple.p = ar.ple.pStart;
qFitReset = ar.qFit(jk);

fprintf('PLE #%i for %s...\n', jk, ar.pLabel{jk});
arWaitbar(0);
tic;

dp_min = 1e-6;
dp_max = 1;

% fix parameter of interest
ar.qFit(jk) = 0;

% calculate initial chi^2
ar = arCalcMerit(ar, true, ar.p(ar.qFit==1));
chi2Reset = arGetMerit(true);
ar.ple.chi2Reset(jk) = chi2Reset;

[chi2sup, psup, errorsup] = ple_task(ar, jk, n, 1, chi2Reset, pReset, dp_min, dp_max);
[chi2sdown, psdown, errorsdown] = ple_task(ar, jk, n, -1, chi2Reset, pReset, dp_min, dp_max);

ar.ple.chi2s{jk} = [fliplr(chi2sdown) chi2Reset chi2sup];
ar.ple.ps{jk} = [flipud(psdown); pReset(:)'; psup];
ar.ple.errors{jk} = [fliplr(errorsdown) nan errorsup];
ar.ple.run(jk) = 1;

arWaitbar(-1);
fprintf('PLE #%i %i errors, %i fit issues, %s elapse time\n', jk, ...
    sum(ar.ple.errors{jk}<0), ...
    sum(ar.ple.errors{jk}==0), ...
    secToHMS(toc));

ar.p = pReset;
ar.qFit(jk) = qFitReset;
ar = arCalcMerit(ar, true, ar.p(ar.qFit==1));

if(nargout>0 && ~qglobalar)
    varargout{1} = ar;
else
    varargout = {};
end



function [chi2s, ps, errors] = ple_task(ar, jk, n, direction, chi2Reset, pReset, dp_min, dp_max)

chi2s = nan(1,n);
ps = nan(n,length(pReset));
errors = nan(1,n);

dchi2 = chi2inv(1-ar.ple.alpha, ar.ple.ndof);
rel_increase = 0.5;

dp = 0.1;
ar.p = pReset;

chi2Last = chi2Reset;
pLast = pReset;

for j=1:n
    if(direction>0)
        arWaitbar(j,2*n, sprintf('PLE up for %s', strrep(ar.pLabel{jk},'_','\_')));
    else
        arWaitbar(j+n,2*n, sprintf('PLE down for %s', strrep(ar.pLabel{jk},'_','\_')));
    end
   
    ar = arCalcMerit(ar, true, ar.p(ar.qFit==1));
    pTrial = findTrial(jk, pLast, dp, direction);
    
    try
        ar.p = pTrial;
        ar = arCalcMerit(ar, true, ar.p(ar.qFit==1));
        while(arGetMerit(true) - chi2Last > rel_increase * dchi2)
            dp = dp/2;
            
            if(dp < dp_min)
                error('minimum step size reached'); 
            end
            
            pTrial = findTrial(jk, pLast, dp, direction);
            
            ar.p = pTrial;
            ar = arCalcMerit(ar, true, ar.p(ar.qFit==1));
        end
        
        ar = arFit(ar, true);
        
        ps(j,:) = ar.p;
        chi2s(j) = arGetMerit(true);
        errors(j) = ar.fit.exitflag;
    catch errorid
        fprintf('PLE #%i ERROR %s\n', jk, errorid.message); 
        errors(j) = -5;
        break;
    end
    
    if(arGetMerit(true) - chi2Reset > 2*dchi2)
        if(direction>0)
            fprintf('PLE #%i reached upper confidence limit\n', jk);
        else
            fprintf('PLE #%i reached lower confidence limit\n', jk);
        end
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
    chi2Last = arGetMerit(true);
    
    if(dp * 2 < dp_max)
        dp = dp * 2;
    end
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
