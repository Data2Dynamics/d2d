% arFitSome(ips, [silent])
% 
% Fitting only a specific subset of model parameters
% 
%   ips     parameter indices which are fitted
%   silent  [false]  argument passed to arFit which might suppress output
%                    at the command line
% 
% The standard flag ar.qFit used to specify which parameters are fitted are
% evaluated on top. This means only parameters with ar.qFit==1 and in ips
% are fitted.
% 
% Example:
% 
% ar.p(:) = -1; % set all parameters to -1
% arFitSome(3:4)
% arPrint
% 
% See also arFitSingle

function arFitSome(ips, silent)
global ar

if(nargin<2)
    silent = false;
end

qFitReset = ar.qFit + 0;

qDoFit = ismember(1:length(ar.p),ips);
ar.qFit(ar.qFit==1 & ~qDoFit) = 0;
try	
	arFit(silent);
catch err
    ar.qFit = qFitReset;
    error(err.message)
end

ar.qFit = qFitReset;
