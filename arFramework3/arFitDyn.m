% Fit only dynamical model parameters
%
% arFitDyn(silent)

function arFitDyn(silent)
global ar

if(nargin==0)
    silent = false;
end

qFitReset = ar.qFit + 0;

ar.qFit(ar.qFit==1 & ar.qDynamic==0) = 0;
arFit(silent);

ar.qFit = qFitReset;
