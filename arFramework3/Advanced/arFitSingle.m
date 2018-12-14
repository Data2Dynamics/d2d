% arFitSingle(j, [silent])
% 
% Fit single parameter
%
%   silent  [false]  argument passed to arFit which might suppress output
%                    at the command line
% 
% Example:
% 
% ar.p(:) = -1; % set all parameters to -1
% arFitSingle(3)
% arPrint
% 
% See also arFitSome


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
