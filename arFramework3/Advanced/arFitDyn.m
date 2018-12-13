% arFitDyn(silent)
%
% Fit only dynamic model parameters
%
%   silent - boolean, suppress output [false]

function arFitDyn(silent)
global ar

if(nargin==0)
    silent = false;
end

qFitReset = ar.qFit + 0;

ar.qFit(ar.qFit==1 & ar.qDynamic==0) = 0;
try	
	arFit(silent);
catch err
    ar.qFit = qFitReset;
    error(err.message)
end

ar.qFit = qFitReset;
