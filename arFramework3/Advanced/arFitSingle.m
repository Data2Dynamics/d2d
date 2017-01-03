% Fit single parameter
%
% arFitSingle(j, silent)

function arFitSingle(j, silent)
global ar

if(nargin==0)
    silent = false;
end

qFitReset = ar.qFit + 0;

ar.qFit(ar.qFit==1) = 0;
ar.qFit(j) = 1;
try	
	arFit(silent);
catch err
    ar.qFit = qFitReset;
    error(err.message)
end

ar.qFit = qFitReset;
