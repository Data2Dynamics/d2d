% arFitObs([silent])
% 
% Fit only observational model parameters.
% 
%   silent      [false]  logical passed to arFit
% 
% Observation parameters are defined in this function as non-dynamic
% parameters which means that it might also comparise parameters of error
% models. 
% 
% This function is called by arTuner if the respective buttons are
% pressed.
%
% Example:
% ar.p(:) = -1; % set all parameters to -1
% arFitObs
% arPrint
% 
% see also arFit, arTuner

function arFitObs(silent)
global ar

if(nargin==0)
    silent = false;
end

qFitReset = ar.qFit + 0;

ar.qFit(ar.qFit==1 & ar.qDynamic==1) = 0;
% arCalcMerit(true,[]);  % I think, this is not requires => don't do it because of runtime.
try	
	arFit(silent);
catch err
    ar.qFit = qFitReset;
    rethrow(err)
end

ar.qFit = qFitReset;