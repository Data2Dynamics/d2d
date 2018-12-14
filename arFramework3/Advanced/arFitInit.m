% arFitInit(silent)
%
% Fit only initial value parameters (ar.qInitial=1 && ar.qFit==1)
%
%   silent - boolean to suppress output [false]

function arFitInit(silent)
global ar

if(nargin==0)
    silent = false;
end

qFitReset = ar.qFit + 0;

ar.qFit(ar.qFit==1 & ar.qInitial==0) = 0;
try	
	arFit(silent);
catch err
    ar.qFit = qFitReset;
    error(err.message)
end
	
ar.qFit = qFitReset;
