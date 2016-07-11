% Initialization of L1 scan
% jks    relative parameters to be investigated by L1 regularization

function l1Init(jks, means, lbs, ubs)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('jks','var') || isempty(jks))
    arPrint
    fprintf('\n')
    error('please specify L1 parameters from above')
end

if(~exist('means','var') || isempty(means))
    means = zeros(size(jks));
end

if(~exist('lbs','var') || isempty(lbs))
    lbs = ones(size(jks))*(-5);
end

if(~exist('ubs','var') || isempty(ubs))
    ubs = ones(size(jks))*5;
end

ar.mean(jks) = means;
ar.lb(jks) = lbs;
ar.ub(jks) = ubs;

ar.type(jks) = 3;
ar.qFit(jks) = 1;
ar.std(jks) = Inf;

try
    arPrint(jks)
end
    
ar.L1ps = [];
ar.L1chi2s = [];
ar.L1chi2fits = [];

ar.L1ps_unpen = [];
ar.L1chi2s_unpen = [];