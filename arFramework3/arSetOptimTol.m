% Set optimization tolerances
% MATLAB optimization routines terminate when the relative change of the
% log-likelihood is less than a given threshold. For likelhood-ratio tests
% however only absolute differences are meaningful. As the likelihood
% scales with the number of fitted data points the standard tolerances are
% not recommended.
% This functions adjusts the tolerances accordingly and should be called
% once before fitting.
%
% arSetOptimTol(tolFun,tolX)
%
% tolFun    [1e-6]
% tolX      [1e-6]

function arSetOptimTol(tolFun,tolX)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('tolFun','var'))
    tolFun = 1e-6;
end
if(~exist('tolX','var'))
    tolX = 1e-6;
end

tmpreserr = 0;

for jm=1:length(ar.model)
    arSimu(ar, 0, ~isfield(ar.model(jm), 'data'))
    
    for jd = 1:length(ar.model(jm).data)
        tmpreserr = tmpreserr + sum(sum(ar.model(jm).data(jd).reserr.^2));
    end
end

oom = floor(log10(tolFun / tmpreserr)); % order of magnitude
ar.config.optim.TolFun = round((tolFun / tmpreserr) / 10^oom) * 10^oom;

fprintf('Successfully set ar.config.optim.TolFun = %g\n', ar.config.optim.TolFun);

oom = floor(log10(tolX / sqrt(sum(ar.qFit==1)))); % order of magnitude
ar.config.optim.TolX = round((tolX / sqrt(sum(ar.qFit==1))) / 10^oom) * 10^oom;

fprintf('Successfully set ar.config.optim.TolX   = %g\n', ar.config.optim.TolX);
