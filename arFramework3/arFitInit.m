% Fit only initial values parameters
%
% arFitInit(silent)

function arFitInit(silent)
global ar

if(nargin==0)
    silent = false;
end

qFitReset = ar.qFit + 0;

ar.qFit(ar.qFit==1 & ar.qInitial==0) = 0;
arFit(silent);

ar.qFit = qFitReset;
