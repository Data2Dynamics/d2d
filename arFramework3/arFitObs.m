% Fit only observational model parameters
%
% arFitObs(silent)

function arFitObs(silent)
global ar

if(nargin==0)
    silent = false;
end

qFitReset = ar.qFit + 0;

ar.qFit(ar.qFit==1 & ar.qDynamic==1) = 0;
arFit(silent);

ar.qFit = qFitReset;