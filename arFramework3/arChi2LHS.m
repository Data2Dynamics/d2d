% chi2 sequence using latin hyper cube sampling
%
% arChi2LHS(n, sensis)
%
% n:        number of runs      [10]

function arChi2LHS(n, sensis, silent)

global ar

if(~exist('n','var'))
    n = 10;
end

if(~exist('sensis','var'))
    sensis = false;
end


if(~exist('silent','var'))
    silent = false;
end

ps = ones(n,1) * ar.p;
psrand = lhsdesign(n,sum(ar.qFit==1));

psrand = psrand .* (ones(n,1)*(ar.ub(ar.qFit==1) - ar.lb(ar.qFit==1)));
psrand = psrand + (ones(n,1)*ar.lb(ar.qFit==1));

ps(:,ar.qFit==1) = psrand;

ps(1,:) = ar.p;

arChi2s(ps, sensis, silent);
